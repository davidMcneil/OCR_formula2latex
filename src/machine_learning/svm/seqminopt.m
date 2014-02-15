function [alphas offset] = seqminopt(data, targetLabels, boxConstraints, ...
    kernelFunc, smoOptions)
%SEQMINOPT Sequential Minimal Optimization for SVM functions
%
%   [ALPHA OFFSET] = SEQMINOPT(TRAINING,Y,BOXC,KERNELFUNC,SMOOPTIONS)
%   trains SVM using Sequential Minimal Optimization. TRAINING is a numeric
%   matrix of predictor data. Rows of TRAINING correspond to observations;
%   columns correspond to features. Y is a column vector of labels. Each
%   values of Y is +1 or -1 representing which group the corresponding data
%   point belongs to. Two groups must be present in Y. BOXC is a column
%   vector of numbers, each of which must be larger than the tolerance with
%   which the KKT conditions are checked. KERNELFUNC is the kernel function
%   used. SMOOPTIONS is optional struct with a number of options. ALPHA is
%   a column vector containing the optimized Lagrange multipliers, one for
%   each data point. Each value in ALPHA is guaranteed to be within [0,C_i]
%   where C_i is the boxconstraint for the corresponding data. OFFSET is
%   the bias of the hyperplane.
%
% References:
% [1] J.C. Platt: A Fast Algorithm for Training  Support Vector
%     Machines,  http://research.microsoft.com/users/jplatt/smo.html
% [2] R.-E. Fan, P.-H. Chen, and C.-J. Lin. Working Set Selection Using
%     Second Order Information for Training SVM. Journal of Machine
%     Learning Research, 6(2005), 1889-1918.
% [3] L. Bottou and C.-J. Lin. Support Vector Machine Solvers. 2006,
%     available from http://www.csie.ntu.edu.tw/~cjlin/papers.html

% This function implements a sequential minimization algorithm to solve the
% optimization problem arising in support vector machine learning. The
% basic idea [1] is to pick a working set of just two alphas and solve the
% optimization for these two alphas while keeping all other variables
% fixed. This is done in the function 'updateAlphas' below. The working set
% selection is done according to the 'maximum gain' method [2] and is
% contained in the main loop below. The third element of the algorithm
% consists of a caching strategy for the kernel matrix for medium sized and
% large problems (when the kernel matrix is too large to fit into memory)
% in order to keep the number of kernel evaluations minimal. A further
% property is the ability of the algorithm to stop even if not all KKT
% conditions are satisfied. This can be specified in the options struct
% which contains a field with the fraction of points that are allowed to
% violate KKT. This is useful for big, noisy problems, where the algorithm
% sometimes spends a lot of time fluctuating very near the optimal
% solution. Note that this algorithm doesn't implement shrinking which is a
% heuristic method useful for speeding up the algorithm for large scale
% problems.
%

%   Copyright 2004-2012 The MathWorks, Inc.


% create standard smo options struct if not passed by user:
narginchk(4,5);
if nargin == 4
    smoOptions = svmsmoset;
end
switch(smoOptions.Display)
    case 'iter'
        smoOptions.verbose = 2;
    case 'final'
        smoOptions.verbose = 1;
    otherwise
        smoOptions.verbose = 0;
end
% Check consistency of input data
if length(size(data)) ~= 2
    error(message('stats:svmtrain:InvalidTraining'));
end

nPoints = size(data, 1);

% targetLabels: a column vector of +1 and -1 of length 'nPoints'
if ~isequal(size(targetLabels), [nPoints 1])
    error(message('stats:svmtrain:InvalidYSize'));
end

if ~all(boxConstraints >= smoOptions.TolKKT)
    error(message('stats:svmtrain:InvalidInput'));
end

% pass over to implementation function:
[alphas offset] = seqminoptImpl(data, targetLabels, ...
    boxConstraints, kernelFunc,  smoOptions);

% Check post conditions: we must have 0<=alpha<= boxconstraints ...
for i = 1:length(alphas)
    if ~( (alphas(i) >= -eps) && (alphas(i) <= boxConstraints(i) + eps))
        error(message('stats:svmtrain:BoxConstraintViolation', ...
                      i, sprintf('%g',alphas(i)), sprintf('%g',boxConstraints(i))));
    end
end
% linear equality constraint: \sum y_i*alpha_i=0
linEqConstraint = targetLabels' * alphas;
if abs(linEqConstraint) > smoOptions.TolKKT
    error(message('stats:svmtrain:ViolatedPostcondition'));
end
end % end of svmtrain

function [alphas offset] = seqminoptImpl(data, targetLabels, ...
    boxConstraints, kernelFunc, smoOptions)
%SEQMINOPTIMPLMOD: Implementation of Sequential Minimal Optimization

% -------------------------------------
% INITIALIZATION
% -------------------------------------

% Working copy of tolerance by which the KKT conditions and stopping
% criteria are checked
tolKKT = smoOptions.TolKKT;

% Tolerance by which support vectors are identified (compare svmtrain function):
svTol = sqrt(eps);

% Number of data points:
nPoints = length(targetLabels);

% Number of accepted KKT violations
acceptedKKTviolations = smoOptions.KKTViolationLevel * nPoints;

% Kernel cache: see createKernelCache function for more documentation
fullKernel = [];  % initialized by call to 'createKernelCache'.
kernelCache = createKernelCache();
kerDiag = kernelCache.getKernelDiag();

% Initialize alphas and gradient of objective function:
alphas = zeros(nPoints, 1);
objGrad = ones(nPoints, 1);
offset = NaN;  % Note that [1] uses a different sign convention here.

% These quantities are useful for the working set selection, see equations
% (7) and (11) of [3].
Avec = zeros(size(targetLabels));
Avec(targetLabels==-1) = - boxConstraints(targetLabels==-1);
Bvec = zeros(size(targetLabels));
Bvec(targetLabels==1) = boxConstraints(targetLabels==1);
upMask = targetLabels .* alphas < (Bvec - svTol);
downMask = targetLabels .* alphas > (Avec + svTol);
idxHelper = 1:nPoints;

% -------------------------------------
% MAIN LOOP
% -------------------------------------

% The main loop finds a pair of alphas as working set using the maximum
% gain method, see WSS3 in [2].

for itCount = 1 : smoOptions.MaxIter

    % find the first alpha:
    [val1 idx1] = max(targetLabels(upMask) .* objGrad(upMask));
    tmp = idxHelper(upMask);
    idx1 = tmp(idx1);

    % find the second alpha according to maximum gain criterion:
    [val2 idx2] = getMaxGain(val1, idx1);

    if isempty(idx2)
         % Max gain method failed, so we select the 'maximum violating pair'
        % in this case:
        [val2 idx2] = min(targetLabels(downMask) .* objGrad(downMask));
        tmp = idxHelper(downMask);
        idx2 = tmp(idx2);
    end  % end if isempty(idx2)
    
    if val1 - val2 <= tolKKT
        offset = (val1 + val2) / 2;  %% see eq. (11) and (14) of [3]
        if smoOptions.verbose
            disp(getString(message('stats:svmtrain:SMOFinished')));
            reportStatus(itCount, val1 - val2);
        end
        return;
    end
   % We now have the working set. Next do the analytical solution:
    updateAlphas(idx1, idx2);
        
    % Regularly check if number of KKT violations is below acceptable limit
    % and print progress in verbose mode:
    if mod(itCount, 500) == 0
        val1 = max(targetLabels(upMask) .* objGrad(upMask));
        val2 = min(targetLabels(downMask) .* objGrad(downMask));
        offset =  (val1 + val2) / 2;  % needed in checkKKT function !
        kktViolationCount = sum(~checkKKT());
        % check number of KKT violations, if KKTViolationLevel > 0 has been set.
        if (acceptedKKTviolations > 0) &&  ...
                (kktViolationCount <= acceptedKKTviolations)
            if smoOptions.verbose
                disp(getString(message('stats:svmtrain:SMOFinished')));
                reportStatus(itCount, val1 - val2);
            end
            return;
        end
        if smoOptions.verbose == 2
            reportStatus(itCount, val1 - val2);
        end
    end
end % end of itCount loop

% Error exit because we never reached the right exit conditions
error(message('stats:svmtrain:NoConvergence'))
   

% -------------------------------------
% HELPER FUNCTIONS
% -------------------------------------

    function [val2 idx2] = getMaxGain(val1, idx1)
        % Given a first index this function finds a second index (which defines
        % the working set together with the first index) based on the 'maximum
        % gain method'. This function either returns a valid index idx2 or an
        % empty index idx2, if it fails to find a second index.
        % The function never returns idx2=NaN.

        % Create the mask of usable down-indices. Note that this definition
        % makes sure that idx2 ~= idx1 in the end.
        mask = downMask & (targetLabels .* objGrad < val1);
        gainNumerator = (targetLabels(mask) .* objGrad(mask) - val1) .* ...
            (targetLabels(mask) .* objGrad(mask) - val1);
        if isempty(fullKernel)
            kerCol1 = kernelCache.getColumn(idx1);
            gainDenominator = -4 * kerCol1(mask) + 2 * kerDiag(mask) + 2 * ...
                kerDiag(idx1);
        else
            gainDenominator = -4 * fullKernel(mask, idx1) + 2 * kerDiag(mask) + 2 * ...
                kerDiag(idx1);
        end
        % Note that the gainNumerator array is strictly positive.
        % So we cannot have a '0/0' situation which would possibly result in
        % idx2 = NaN.
        [~, idx2] = max(gainNumerator ./ gainDenominator);
        tmp = idxHelper(mask);
        idx2 = tmp(idx2);
        val2 = targetLabels(idx2)' * objGrad(idx2);
    end  % end of getMaxGain

    function updateAlphas(i, j)
        % This function calculates new values for alpha_i and alpha_j
        % according to eqs. (16) and (18) of [1].
        % j corresponds to index 2 and i corresponds to index 1.

        % Get relevant kernel matrix elements: eq. (15) of [1].
        if isempty(fullKernel)
            [Kii Kjj Kij] = kernelCache.getElements(i,j);
            eta = Kii + Kjj - 2 * Kij;
        else
            eta = fullKernel(i,i) + fullKernel(j,j) - 2 * fullKernel(i,j);
        end

        % Calculate clip limits: eq. (13) and (14) of [1].
        if targetLabels(i) * targetLabels(j) == 1
            Low = max(0, alphas(j) + alphas(i) - boxConstraints(i));
            High = min(boxConstraints(j), alphas(j) + alphas(i));
        else
            Low = max(0, alphas(j) - alphas(i));
            High = min(boxConstraints(j), boxConstraints(i) + alphas(j) - alphas(i));
        end

        if eta > eps
            % Calculate new values for alpha(i) and alpha(j) (but don't store
            % them yet in the global alpha array). This corresponds to
            % finding the right orientation of the separating plane,
            % eq. (16) of [1].
            lambda = - targetLabels(i) * objGrad(i) + targetLabels(j) * objGrad(j);
            alpha_j = alphas(j) + targetLabels(j) / eta * lambda;
            % Clip alpha: eq. (17) of [1].
            if alpha_j < Low
                alpha_j = Low;
            elseif alpha_j > High
                alpha_j = High;
            end
        else
            % The case 'eta < eps' should not happen too often (only for duplicate
            % data points and illegal kernels)! In this case we do
            % the ugly calculation of eq. (19) in [1]
            [psi_l, psi_h Low High] = evalPsiAtEnd(i,j, Low, High);
            if psi_l < (psi_h - eps)
                alpha_j  = Low;
            elseif psi_l > (psi_h + eps)
                alpha_j = High;
            else % no progress :-(
                alpha_j = alphas(j);
            end
        end

        % New value for alpha(i): eq. (18) of [1]
        alpha_i = alphas(i) + targetLabels(j) * targetLabels(i) * ...
            (alphas(j) - alpha_j);

        % We do have to make sure that the new alpha definitely
        % sits in its box.
        if alpha_i < eps
            alpha_i = 0;
        elseif alpha_i > (boxConstraints(i) - eps)
            alpha_i = boxConstraints(i);
        end

        % Update gradient of objective function:
        if isempty(fullKernel)
            objGrad = objGrad - ...
                (kernelCache.getColumn(i) .* targetLabels) * ...
                (alpha_i - alphas(i)) * targetLabels(i) - ...
                (kernelCache.getColumn(j) .* targetLabels) * ...
                (alpha_j - alphas(j)) * targetLabels(j);
        else
            objGrad = objGrad - ...
                (fullKernel(:,i) .* targetLabels) * ...
                (alpha_i - alphas(i)) * targetLabels(i) - ...
                (fullKernel(:,j) .* targetLabels) * ...
                (alpha_j - alphas(j)) * targetLabels(j);
        end

        % Update up and down masks:
        upMask(i) = targetLabels(i) * alpha_i < (Bvec(i) - svTol);
        downMask(i) = targetLabels(i) * alpha_i > (Avec(i) + svTol);
        upMask(j) = targetLabels(j) * alpha_j < (Bvec(j) - svTol);
        downMask(j) = targetLabels(j) * alpha_j > (Avec(j) + svTol);

        % Finally, update global alpha array
        alphas(i) = alpha_i;
        alphas(j) = alpha_j;
    end % end of updateAlphas

    function [psi_l, psi_h Low High] = evalPsiAtEnd(i, j, Low, High)
        % This function evaluates the objective function at the ends of
        % the feasible region. This is necessary in the hopefully rare
        % case of eta < 0. It looks a bit ugly but simply implements eq. (19) of [1].

        [Kii Kjj Kij] = kernelCache.getElements(i,j);

        s = targetLabels(i) * targetLabels(j);
        fi = - objGrad(i) - alphas(i) * Kii - s * alphas(j) * Kij;
        fj = - objGrad(j) - alphas(j) * Kjj - s * alphas(i) * Kij;
        Li = alphas(i) + s * (alphas(j) - Low);
        Hi = alphas(i) + s * (alphas(j) - High);
        psi_l = Li * fi + Low * fj + Li * Li * Kii / 2 + Low * Low * Kjj / 2 + ...
            s * Low * Li * Kij;
        psi_h = Hi * fi + High * fj + Hi * Hi * Kii / 2 + High * High * Kjj / 2 + ...
            s * High * Hi * Kij;
    end

    function flags = checkKKT()
        % This function checks which alphas satisfy the KKT conditions.
        % The return value is a logical array of the same size as alphas
        % and indicates which alphas satisfy the KKT conditions.

        amount = - objGrad + targetLabels * offset;
        flags = false(size(amount));

        % SV: check y_i * u_i == 1
        freeSVmask = (alphas > svTol) & (alphas < (boxConstraints - svTol));
        flags(freeSVmask) = (abs(amount(freeSVmask)) < smoOptions.TolKKT);

        % Lower bound alphas: correctly classified and not an SV:
        % check u_i * y_i >= 1
        mask = alphas < svTol;
        flags(mask) = amount(mask) > -smoOptions.TolKKT;

        % Upper bound alphas: margin violators:
        % check u_i * y_i <= 1
        mask = (boxConstraints - alphas) < svTol;
        flags(mask) = amount(mask) <= smoOptions.TolKKT;
    end

    function reportStatus(numIter, stopCrit)
        % This function reports the current status of the algorithm and is called
        % only in verbose mode.

        % Calculate value of objective function
        objFunc = sum(alphas) + 1/2 * (objGrad'-1) * alphas;

        disp(getString(message('stats:svmtrain:SMOStatus')));
        fprintf('%s',getString(message('stats:svmtrain:NumberOfIterations', numIter)));
        fprintf('%s',getString(message('stats:svmtrain:ValueOfStoppingCriterion', sprintf('%f',stopCrit))));
        fprintf('%s',getString(message('stats:svmtrain:ValueOfObjectiveFunction', sprintf('%f',objFunc))));
        if isempty(fullKernel)
            fprintf('%s\n',getString(message('stats:svmtrain:NumberOfCachedKernelColumns', ...
                length(kernelCache.getStoredColIdxs()),kernelCache.getCacheSize() )));              
        end
        disp('---------------------------');
    end

% -----------------------------------
% Kernel cache closure
% -----------------------------------
    function kernelCache = createKernelCache()
        % This function implements a caching mechanism for the kernel
        % matrix as a closure.
        % If the number of data points is smaller than smoOptions.kernelCacheLimit
        % then the full kernel is stored in the variable 'fullKernel' and everything
        % is very simple. The 'fullKernel' variable is function-global and not
        % hidden in this closure in order to have the fastest possible element access.
        % Also the diagonal of the kernel matrix is stored (redundantly if also
        % the full kernel matrix is stored) in the function-global variable
        % 'kerDiag' in order to have fastest possible access to subsets of the
        % kernel diagonal.
        %
        % If the number of data points is larger than smoOptions.kernelCacheLimit
        % the function-global variable 'fullKernel' is left empty and only a subset
        % of the kernel matrix is stored inside this closure. In this case element
        % access is done through the member functions of this kernel. However,
        % also in this case the function-global variable 'kerDiag' contains the
        % full kernel diagonal.

        % Internally, a subset of the kernel columns is stored in 'subKernel'.
        % The idea is that 'subKernel' stores the most recently used kernel columns.
        % If a kernel column not contained in 'subKernel' is requested, the
        % column of 'subKernel' that has not been used for the longest time is
        % replaced by the newly requested column. To achieve this an array of most
        % recently requested kernel indices is maintained in the variable
        % 'leastIndices'.

        % Private fields:
        % -----------------------
        % If the full kernel is not stored in the 'fullKernel' variable we here
        % keep its diagonal. It is also stored in the function-global variable 'kerDiag'.
        kernelDiag = [];

        % Subset of kernel columns, i.e. the number of rows equals the number of
        % rows of the original kernel, but it contains less columns.
        subKernel = [];

        % An array of size nPoints: an element contains 0 if the column with
        % this index is not contained in subKernel. Otherwise it contains the
        % index of that column in subKernel.
        subKernelIndices = [];

        % An array of size equal to the number of columns in subKernel,
        % maintaining the order of most recently used kernel columns.
        % The last element contains the last recently used column in subKernel.
        % and the first element contains the most recently used.
        % The whole stuff gets filled as kernel columns are requested.
        % Before, leastIndices contains 0 elements.
        leastIndices = [];

        % Constructor section
        % ------------------------------
        % fill the full kernel matrix if possible
        if nPoints <= smoOptions.KernelCacheLimit
            try
                fullKernel = kernelFunc(data, data);
            catch ME
                error(message('stats:svmtrain:KernelFunctionError', ME.message));
            end
        else
            try 
            kernelDiag = arrayfun(@(idx) kernelFunc(data(idx, :), data(idx, :)), ...
                (1:nPoints)');
            catch ME
                error(message('stats:svmtrain:KernelFunctionError', ME.message));
            end
                
            s = max(1, floor(smoOptions.KernelCacheLimit^2 / nPoints));
            leastIndices = zeros(1,s);
            subKernel = zeros(nPoints, s);
            subKernelIndices = zeros(1, nPoints);
        end

        % Private Methods
        % ------------------------
        function subKernelIndex = loadColumn(idx)
            % This function should must only be called if a kernel column not
            % contained in the cache is requested. Then this function calculates
            % the requested column, inserts it into the 'subKernel' matrix and
            % updates the 'subKernelIndices' and 'leastIndices' variables.

            % Determine position where the new kernel column in subKernel
            % should be stored:
            deleteIndex = leastIndices(end);
            if deleteIndex == 0
                % subKernel not yet completely filled
                subKernelIndex = 1 + max(subKernelIndices);
            else
                subKernelIndex = subKernelIndices(deleteIndex);
                subKernelIndices(deleteIndex) = 0;
            end

            % This is where we put the new column
            subKernelIndices(idx) = subKernelIndex;
            % Calculate new kernel column and put it into column with
            % index 'kernelColIdx'.
            try
            subKernel(:,subKernelIndex) = kernelFunc(data, data(idx, :));
            catch ME
            error(message('stats:svmtrain:KernelFunctionError', ME.message));
           
            end
            % Finally, add idx to leastIndices array
            leastIndices(2:end) = leastIndices(1:end-1);
            leastIndices(1) = idx;
        end

        % Public Methods
        % ------------------------
        function [Kii Kjj Kij] = getElements(i, j)
            % This function provides some often used kernel elements.

            % Case 1: the full kernel matrix is stored:
            if ~isempty(fullKernel)
                Kii = fullKernel(i,i);
                Kjj = fullKernel(j,j);
                Kij = fullKernel(i,j);
            else
                Kii = kernelDiag(i);
                Kjj = kernelDiag(j);
                if subKernelIndices(i) > 0
                    Kij = subKernel(j, subKernelIndices(i));
                elseif subKernelIndices(j) > 0
                    Kij = subKernel(i, subKernelIndices(j));
                else
                    Kij = subKernel(i, loadColumn(j));
                end
            end
        end  % end of getElementsGeneral

        function kerCol = getColumn(colIdx)
            % This function returns the kernel column with column index colIdx.

            if ~isempty(fullKernel)
                kerCol = fullKernel(:, colIdx);
            else
                if subKernelIndices(colIdx) == 0
                    kerCol = subKernel(:, loadColumn(colIdx));
                else
                    kerCol = subKernel(:, subKernelIndices(colIdx));
                end
            end
        end

        function ret = getKernelDiag()
            % return the diagonal of the kernel
            if ~isempty(fullKernel)
                ret = diag(fullKernel);
            else
                ret = kernelDiag;
            end
        end

        function ret = getStoredColIdxs()
            % returns the indices of the columns contained in the cache.
            % Just needed for reporting and debugging.
            if ~isempty(fullKernel)
                ret = 1:size(fullKernel,2);
            else
                ret = find(subKernelIndices > 0);
            end
        end

        function ret = getCacheSize()
            % returns the number of kernel columns that could possibly fit into the cache.
            % It is just needed for reporting.
            if ~isempty(fullKernel)
                ret = 1:size(fullKernel,2);
            else
                ret = length(leastIndices);
            end
        end

        % Define return values
        %--------------------------
        kernelCache.getElements = @getElements;
        kernelCache.getColumn = @getColumn;
        kernelCache.getKernelDiag = @getKernelDiag;
        kernelCache.getStoredColIdxs = @getStoredColIdxs;
        kernelCache.getCacheSize = @getCacheSize;
    end  % end of createKernelCache

end  % end of seqMinOptImpl


