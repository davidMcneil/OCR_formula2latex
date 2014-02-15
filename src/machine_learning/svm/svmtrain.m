function [svm_struct, svIndex] = svmtrain(training, groupnames, varargin)
%SVMTRAIN Train a support vector machine classifier
%   SVMSTRUCT = SVMTRAIN(TRAINING, Y) trains a support vector machine (SVM)
%   classifier on data taken from two groups. TRAINING is a numeric matrix
%   of predictor data. Rows of TRAINING correspond to observations; columns
%   correspond to features. Y is a column vector that contains the known
%   class labels for TRAINING. Y is a grouping variable, i.e., it can be a
%   categorical, numeric, or logical vector; a cell vector of strings; or a
%   character matrix with each row representing a class label (see help for
%   groupingvariable). Each element of Y specifies the group the
%   corresponding row of TRAINING belongs to. TRAINING and Y must have the
%   same number of rows. SVMSTRUCT contains information about the trained
%   classifier, including the support vectors, that is used by SVMCLASSIFY
%   for classification. SVMTRAIN treats NaNs, empty strings or 'undefined'
%   values as missing values and ignores the corresponding rows in
%   TRAINING and Y.
%
%   SVMSTRUCT = SVMTRAIN(TRAINING, Y, 'PARAM1',val1, 'PARAM2',val2, ...)
%   specifies one or more of the following name/value pairs:
%
%      Name                Value
%      'kernel_function'  A string or a function handle specifying the
%                         kernel function used to represent the dot
%                         product in a new space. The value can be one of
%                         the following:
%                         'linear'     - Linear kernel or dot product
%                                        (default). In this case, SVMTRAIN
%                                        finds the optimal separating plane
%                                        in the original space.
%                         'quadratic'  - Quadratic kernel
%                         'polynomial' - Polynomial kernel with default
%                                        order 3. To specify another order,
%                                        use the 'polyorder' argument.
%                         'rbf'        - Gaussian Radial Basis Function
%                                        with default scaling factor 1. To
%                                        specify another scaling factor,
%                                        use the 'rbf_sigma' argument.
%                         'mlp'        - Multilayer Perceptron kernel (MLP)
%                                        with default weight 1 and default
%                                        bias -1. To specify another weight
%                                        or bias, use the 'mlp_params'
%                                        argument.
%                         function     - A kernel function specified using
%                                        @(for example @KFUN), or an
%                                        anonymous function. A kernel
%                                        function must be of the form
%
%                                        function K = KFUN(U, V)
%
%                                        The returned value, K, is a matrix
%                                        of size M-by-N, where M and N are
%                                        the number of rows in U and V
%                                        respectively.
%
%   'rbf_sigma'           A positive number specifying the scaling factor
%                         in the Gaussian radial basis function kernel.
%                         Default is 1.
%
%   'polyorder'           A positive integer specifying the order of the
%                         polynomial kernel. Default is 3.
%
%   'mlp_params'          A vector [P1 P2] specifying the parameters of MLP
%                         kernel.  The MLP kernel takes the form:
%                         K = tanh(P1*U*V' + P2),
%                         where P1 > 0 and P2 < 0. Default is [1,-1].
%
%   'method'              A string specifying the method used to find the
%                         separating hyperplane. Choices are:
%                         'SMO' - Sequential Minimal Optimization (SMO)
%                                 method (default). It implements the L1
%                                 soft-margin SVM classifier.
%                         'QP'  - Quadratic programming (requires an
%                                 Optimization Toolbox license). It
%                                 implements the L2 soft-margin SVM
%                                 classifier. Method 'QP' doesn't scale
%                                 well for TRAINING with large number of
%                                 observations.
%                         'LS'  - Least-squares method. It implements the
%                                 L2 soft-margin SVM classifier.
%
%   'options'             Options structure created using either STATSET or
%                         OPTIMSET.
%                         * When you set 'method' to 'SMO' (default),
%                           create the options structure using STATSET.
%                           Applicable options:
%                           'Display'  Level of display output.  Choices
%                                    are 'off' (the default), 'iter', and
%                                    'final'. Value 'iter' reports every
%                                    500 iterations.
%                           'MaxIter'  A positive integer specifying the
%                                    maximum number of iterations allowed.
%                                    Default is 15000 for method 'SMO'.
%                         * When you set method to 'QP', create the
%                           options structure using OPTIMSET. For details
%                           of applicable options choices, see QUADPROG
%                           options. SVM uses a convex quadratic program,
%                           so you can choose the 'interior-point-convex'
%                           algorithm in QUADPROG.
%
%  'tolkkt'              A positive scalar that specifies the tolerance
%                        with which the Karush-Kuhn-Tucker (KKT) conditions
%                        are checked for method 'SMO'. Default is
%                        1.0000e-003.
%
%  'kktviolationlevel'   A scalar specifying the fraction of observations
%                        that are allowed to violate the KKT conditions for
%                        method 'SMO'. Setting this value to be positive
%                        helps the algorithm to converge faster if it is
%                        fluctuating near a good solution. Default is 0.
%
%  'kernelcachelimit'    A positive scalar S specifying the size of the
%                        kernel matrix cache for method 'SMO'. The
%                        algorithm keeps a matrix with up to S * S
%                        double-precision numbers in memory. Default is
%                        5000. When the number of points in TRAINING
%                        exceeds S, the SMO method slows down. It's
%                        recommended to set S as large as your system
%                        permits.
%
%  'boxconstraint'       The box constraint C for the soft margin. C can be
%                        a positive numeric scalar or a vector of positive
%                        numbers with the number of elements equal to the
%                        number of rows in TRAINING.
%                        Default is 1.
%                        * If C is a scalar, it is automatically rescaled
%                          by N/(2*N1) for the observations of group one,
%                          and by N/(2*N2) for the observations of group
%                          two, where N1 is the number of observations in
%                          group one, N2 is the number of observations in
%                          group two. The rescaling is done to take into
%                          account unbalanced groups, i.e., when N1 and N2
%                          are different.
%                        * If C is a vector, then each element of C
%                          specifies the box constraint for the
%                          corresponding observation.
%
%   'autoscale'          A logical value specifying whether or not to
%                        shift and scale the data points before training.
%                        When the value is true, the columns of TRAINING
%                        are shifted and scaled to have zero mean unit
%                        variance. Default is true.
%
%   'showplot'           A logical value specifying whether or not to show
%                        a plot. When the value is true, SVMTRAIN creates a
%                        plot of the grouped data and the separating line
%                        for the classifier, when using data with 2
%                        features (columns). Default is false.
%
%   SVMSTRUCT is a structure having the following properties:
%
%   SupportVectors       Matrix of data points with each row corresponding
%                        to a support vector. 
%                        Note: when 'autoscale' is false, this field
%                        contains original support vectors in TRAINING.
%                        When 'autoscale' is true, this field contains
%                        shifted and scaled vectors from TRAINING.
%   Alpha                Vector of Lagrange multipliers for the support
%                        vectors. The sign is positive for support vectors
%                        belonging to the first group and negative for
%                        support vectors belonging to the second group.
%   Bias                 Intercept of the hyperplane that separates
%                        the two groups.
%                        Note: when 'autoscale' is false, this field
%                        corresponds to the original data points in
%                        TRAINING. When 'autoscale' is true, this field
%                        corresponds to shifted and scaled data points.
%   KernelFunction       The function handle of kernel function used.
%   KernelFunctionArgs   Cell array containing the additional arguments
%                        for the kernel function.
%   GroupNames           A column vector that contains the known
%                        class labels for TRAINING. Y is a grouping
%                        variable (see help for groupingvariable).
%   SupportVectorIndices A column vector indicating the indices of support
%                        vectors.
%   ScaleData            This field contains information about auto-scale.
%                        When 'autoscale' is false, it is empty. When
%                        'autoscale' is set to true, it is a structure
%                        containing two fields:
%                        shift       - A row vector containing the negative
%                                      of the mean across all observations
%                                      in TRAINING.
%                        scaleFactor - A row vector whose value is
%                                      1./STD(TRAINING).
%   FigureHandles        A vector of figure handles created by SVMTRAIN
%                        when 'showplot' argument is TRUE.
%
%   Example:
%       % Load the data and select features for classification
%       load fisheriris
%       X = [meas(:,1), meas(:,2)];
%       % Extract the Setosa class
%       Y = nominal(ismember(species,'setosa'));
%       % Randomly partitions observations into a training set and a test
%       % set using stratified holdout
%       P = cvpartition(Y,'Holdout',0.20);
%       % Use a linear support vector machine classifier
%       svmStruct = svmtrain(X(P.training,:),Y(P.training),'showplot',true);
%       C = svmclassify(svmStruct,X(P.test,:),'showplot',true);
%       errRate = sum(Y(P.test)~= C)/P.TestSize  %mis-classification rate
%       conMat = confusionmat(Y(P.test),C) % the confusion matrix
%
%   See also SVMCLASSIFY, NAIVEBAYES, CLASSREGTREE, CLASSIFY, TREEBAGGER,
%            GROUPINGVARIABLE

%   Copyright 2004-2012 The MathWorks, Inc.


%   References:
%
%     [1] Cristianini, N., Shawe-Taylor, J An Introduction to Support
%         Vector Machines, Cambridge University Press, Cambridge, UK. 2000.
%         http://www.support-vector.net
%     [2] Kecman, V, Learning and Soft Computing,
%         MIT Press, Cambridge, MA. 2001.
%     [3] Suykens, J.A.K., Van Gestel, T., De Brabanter, J., De Moor, B.,
%         Vandewalle, J., Least Squares Support Vector Machines,
%         World Scientific, Singapore, 2002.
%     [4] J.C. Platt: A Fast Algorithm for Training  Support Vector
%         Machines,  Advances in Kernel Methods - Support Vector Learning,
%         MIT Press, 1998.
%     [5] J.C. Platt: Fast Training of Support Vector Machines using
%         Sequential Minimal Optimization Microsoft Research Technical
%         Report MSR-TR-98-14, 1998.
%     [6] http://www.kernel-machines.org/papers/tr-30-1998.ps.gz
%
%   SVMTRAIN(...,'KFUNARGS',ARGS) allows you to pass additional
%   arguments to kernel functions.

narginchk(2, Inf);

% check group is a vector or a char array
if ~isvector(groupnames) && ~ischar(groupnames)
    error(message('stats:svmtrain:GroupNotVector'));
end
% make sure that the data are correctly oriented.
if size(groupnames,1) == 1
    groupnames = groupnames';
end

if ~isnumeric(training) || ~ismatrix(training) 
    error(message('stats:svmtrain:TrainingBadType'));
end

% grp2idx sorts a numeric grouping var ascending, and a string grouping
% var by order of first occurrence
[groupIndex, groupString] = grp2idx(groupnames);

% make sure data is the right size
if size(training,1) ~= size(groupIndex,1)
    if size(training,2) == size(groupIndex,1)
        training = training';
    else
        error(message('stats:svmtrain:DataGroupSizeMismatch'))
    end
end

if isempty(training)
    error(message('stats:svmtrain:NoData'))
end

nans = isnan(groupIndex) | any(isnan(training),2);
if any(nans)
    training(nans,:) = [];
    groupIndex(nans) = [];
end
if isempty(training)
    error(message('stats:svmtrain:NoData'))
end

ngroups = length(unique(groupIndex));
nPoints = length(groupIndex);

if ngroups > 2
    error(message('stats:svmtrain:TooManyGroups', ngroups))
end
if length(groupString) > ngroups
    warning(message('stats:svmtrain:EmptyGroups'));
        
end
% convert to groupIndex from 2 to -1.
groupIndex = 1 - (2* (groupIndex-1));

pnames = {'kernel_function','method','showplot', 'polyorder','mlp_params',...
    'boxconstraint','rbf_sigma','autoscale', 'options',...
    'tolkkt','kktviolationlevel','kernelcachelimit'...
    'kfunargs', 'quadprog_opts','smo_opts'};
dflts =  { 'linear',         [],      false,      [],         [],   ....
    1,              [],         true ,        [] ,    ....
    [],      [],                 [],...
    {} ,          []  ,           []};
[kfun,optimMethod, plotflag, polyOrder, mlpParams, boxC,  rbf_sigma, ...
    autoScale, opts, tolkkt, kktvl,kerCL, kfunargs, qpOptsInput, ...
    smoOptsInput] = internal.stats.parseArgs(pnames, dflts, varargin{:});

usePoly = false;
useMLP = false;
useSigma = false;
%parse kernel functions
if ischar(kfun)
    okfuns = {'linear','quadratic', 'radial','rbf','polynomial','mlp'};
    [~,i] = internal.stats.getParamVal(kfun,okfuns,'kernel_function');
    switch i
        case 1
            kfun = @linear_kernel;
        case 2
            kfun = @quadratic_kernel;
        case {3,4}
            kfun = @rbf_kernel;
            useSigma = true;
        case 5
            kfun = @poly_kernel;
            usePoly = true;
        case 6
            kfun = @mlp_kernel;
            useMLP = true;
    end
elseif ~isa(kfun,  'function_handle')
    error(message('stats:svmtrain:BadKernelFunction'));
end

%parse optimization method
optimList ={'QP','SMO','LS'};
i = 2; % set to 'SMO'

if ~isempty(optimMethod)
    [~,i] = internal.stats.getParamVal(optimMethod,optimList,'Method');
    if i==1 &&  ( ~license('test', 'optimization_toolbox') ...
            || isempty(which('quadprog')))
        warning(message('stats:svmtrain:NoOptim'));
        i = 2;
    end
end

if i == 2 && ngroups==1
    error(message('stats:svmtrain:InvalidY'));
end
optimMethod = optimList{i};

% The large scale solver cannot handle this type of problem, so turn it off.
% qp_opts = optimset('LargeScale','Off','display','off');
% We can use the 'interior-point-convex' option 
qp_opts = optimset('Algorithm','active-set','display','off');
smo_opts = statset('Display','off','MaxIter',15000);
%parse opts. opts will override 'quadprog_opt' and 'smo_opt' argument
if ~isempty(opts)
    qp_opts = optimset(qp_opts,opts);
    smo_opts = statset(smo_opts,opts);
else
    % only consider undocumented 'quadprog_opts' arguments
    % when 'opts' is empty; Otherwise, ignore 'quadprog_opts'
    if ~isempty(qpOptsInput)
        if isstruct(qpOptsInput)
            qp_opts = optimset(qp_opts,qpOptsInput);
        elseif iscell(qpOptsInput)
            qp_opts = optimset(qp_opts,qpOptsInput{:});
        else
            error(message('stats:svmtrain:BadQuadprogOpts'));
        end
    end
end

if ~isempty(smoOptsInput) && isempty(tolkkt) && isempty(kktvl) ...
        && isempty(kerCL) && isempty(opts)
    %back-compatibility.
    smo_opts = svmsmoset(smoOptsInput);
else
    if isempty(tolkkt)
        tolkkt = 1e-3;
    end
    if isempty(kerCL)
        kerCL = 5000;
    end
    if isempty(kktvl)
        kktvl = 0;
    end
    smo_opts = svmsmoset(smo_opts,'tolkkt',tolkkt,'KernelCacheLimit',kerCL,....
        'KKTViolationLevel',kktvl);
end

if ~isscalar(smo_opts.TolKKT) || ~isnumeric(smo_opts.TolKKT) || smo_opts.TolKKT <= 0
    error(message('stats:svmtrain:badTolKKT'));
end

if ~isscalar(smo_opts.KKTViolationLevel) || ~isnumeric(smo_opts.KKTViolationLevel)...
        || smo_opts.KKTViolationLevel < 0 || smo_opts.KKTViolationLevel > 1
    error(message('stats:svmtrain:badKKTVL'));
end

if  ~isscalar(smo_opts.KernelCacheLimit) || ~isnumeric(smo_opts.KernelCacheLimit)...
        ||smo_opts.KernelCacheLimit < 0
    error(message('stats:svmtrain:badKerCL'));
end

%parse plot flag
plotflag = opttf(plotflag,'showplot');
if plotflag && size(training,2) ~=2
    plotflag = false;
    warning(message('stats:svmtrain:OnlyPlot2D'));
end

if ~isempty(kfunargs) &&  ~iscell(kfunargs)
    kfunargs = {kfunargs};
end

%polyOrder
if ~isempty(polyOrder)
    
    %setPoly = true;
    if ~usePoly
        warning(message('stats:svmtrain:PolyOrderNotPolyKernel'));
    else
        kfunargs = {polyOrder};
    end
end

% mlpparams
if ~isempty(mlpParams)
    if ~isnumeric(mlpParams) || numel(mlpParams)~=2
        error(message('stats:svmtrain:BadMLPParams'));
    end
    if mlpParams(1) <= 0
        error(message('stats:svmtrain:MLPWeightNotPositive'))
    end
    if mlpParams(2) >= 0
        warning(message('stats:svmtrain:MLPBiasNotNegative'))
    end
    if ~useMLP
        warning(message('stats:svmtrain:MLPParamNotMLPKernel'));
    else
        kfunargs = {mlpParams(1), mlpParams(2)};
    end
end

%rbf_sigma
if ~isempty(rbf_sigma)
    if useSigma
        kfunargs = {rbf_sigma};
    else
        warning(message('stats:svmtrain:RBFParamNotRBFKernel'))
    end
end

% box constraint: it can be a positive numeric scalar or a numeric vector
% of the same length as the number of data points
if isscalar(boxC) && isnumeric(boxC) && boxC > 0
    % scalar input: adjust to group size and transform into vector
    % set default value of box constraint
    boxconstraint = ones(nPoints, 1); 
    n1 = length(find(groupIndex==1));
    n2 = length(find(groupIndex==-1));
    c1 = 0.5 * boxC * nPoints / n1;
    c2 = 0.5 * boxC * nPoints / n2;
    boxconstraint(groupIndex==1) = c1;
    boxconstraint(groupIndex==-1) = c2;
elseif isvector(boxC) && isnumeric(boxC) && all(boxC > 0) && (length(boxC) == nPoints)
    % vector input
    boxconstraint = boxC;
else
    error(message('stats:svmtrain:InvalidBoxConstraint'));
end
% If boxconstraint == Inf then convergence will not
% happen so fix the value to 1/sqrt(eps).
boxconstraint = min(boxconstraint,repmat(1/sqrt(eps(class(boxconstraint))),...
    size(boxconstraint)));

autoScale = opttf(autoScale,'autoscale');

% plot the data if requested
if plotflag
    [hAxis,hLines] = svmplotdata(training,groupIndex);
    legend(hLines,cellstr(groupString));
end

% autoscale data if required,
scaleData = [];
if autoScale
    scaleData.shift = - mean(training);
    stdVals = std(training);
    scaleData.scaleFactor = 1./stdVals;
    % leave zero-variance data unscaled:
    scaleData.scaleFactor(~isfinite(scaleData.scaleFactor)) = 1;
    
    % shift and scale columns of data matrix:
    for c = 1:size(training, 2)
        training(:,c) = scaleData.scaleFactor(c) * ...
            (training(:,c) +  scaleData.shift(c));
    end
end

if strcmpi(optimMethod, 'SMO')
    % if we have a kernel that takes extra arguments we must define a new
    % kernel function handle to be passed to seqminopt
    if ~isempty(kfunargs)
        tmp_kfun = @(x,y) feval(kfun, x,y, kfunargs{:});
    else
        tmp_kfun = kfun;
    end
    
    [alpha, bias] = seqminopt(training, groupIndex, ...
        boxconstraint, tmp_kfun, smo_opts);
    
    svIndex = find(alpha > sqrt(eps));
    sv = training(svIndex,:);
    alphaHat = groupIndex(svIndex).*alpha(svIndex);
    
else % QP and LS both need the kernel matrix:
    
    % calculate kernel function and add additional term required
    % for two-norm soft margin
    try
        kx = feval(kfun,training,training,kfunargs{:});
        % ensure function is symmetric
        kx = (kx+kx')/2 + diag(1./boxconstraint);
    catch ME
        m = message('stats:svmtrain:KernelFunctionError',func2str(kfun));
        throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
    end
    
    % create Hessian
    H =((groupIndex * groupIndex').*kx);
    
    if strcmpi(optimMethod, 'QP')
        if strncmpi(qp_opts.Algorithm,'inte',4)
            X0 = [];
        else
            X0= ones(nPoints,1);
        end
        [alpha, ~, exitflag, output] = quadprog(H,-ones(nPoints,1),[],[],...
            groupIndex',0,zeros(nPoints,1), Inf *ones(nPoints,1),...
            X0, qp_opts);
        
        if exitflag <= 0
            error(message('stats:svmtrain:UnsolvableOptimization', output.message));
        end
        
        % The support vectors are the non-zeros of alpha.
        % We could also use the zero values of the Lagrangian (fifth output of
        % quadprog) though the method below seems to be good enough.
        svIndex = find(alpha > sqrt(eps));
        sv = training(svIndex,:);
        
        % calculate the parameters of the separating line from the support
        % vectors.
        alphaHat = groupIndex(svIndex).*alpha(svIndex);
        
        % Calculate the bias by applying the indicator function to the support
        % vector with largest alpha.
        [~,maxPos] = max(alpha);
        bias = groupIndex(maxPos) - sum(alphaHat.*kx(svIndex,maxPos));
        % an alternative method is to average the values over all support vectors
        % bias = mean(groupIndex(sv)' - sum(alphaHat(:,ones(1,numSVs)).*kx(sv,sv)));
        
        % An alternative way to calculate support vectors is to look for zeros of
        % the Lagrangian (fifth output from QUADPROG).
        %
        % [alpha,fval,output,exitflag,t] = quadprog(H,-ones(nPoints,1),[],[],...
        %             groupIndex',0,zeros(nPoints,1),inf *ones(nPoints,1),zeros(nPoints,1),opts);
        %
        % sv = t.lower < sqrt(eps) & t.upper < sqrt(eps);
    else  % Least-Squares
        % now build up compound matrix for solver
        A = [0 groupIndex';groupIndex,H];
        b = [0;ones(size(groupIndex))];
        x = A\b;
        
        % calculate the parameters of the separating line from the support
        % vectors.
        sv = training;
        bias = x(1);
        alphaHat = groupIndex.*x(2:end);
        svIndex = (1:nPoints)';
    end
end
svm_struct.SupportVectors = sv;
svm_struct.Alpha = alphaHat;
svm_struct.Bias = bias;
svm_struct.KernelFunction = kfun;
svm_struct.KernelFunctionArgs = kfunargs;
svm_struct.GroupNames = groupnames;
svm_struct.SupportVectorIndices = svIndex;
svm_struct.ScaleData = scaleData;
svm_struct.FigureHandles = [];
if plotflag
    hSV = svmplotsvs(hAxis,hLines,groupString,svm_struct);
    svm_struct.FigureHandles = {hAxis,hLines,hSV};
end
