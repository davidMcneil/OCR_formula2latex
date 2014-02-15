function options = svmsmoset(varargin)
%SVMSMOSET Obsolete function.
%   SVMSMOSET is obsolete. Use SVMTRAIN to create and change SMO OPTIONS.
%   
%   OPTIONS = SVMSMOSET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an
%   options structure OPTIONS in which the named properties have the
%   specified values. Any unspecified properties have default values. It is
%   sufficient to type only the leading characters that uniquely identify
%   the property. Case is ignored for property names.
%
%   OPTIONS = SVMSMOSET(OLDOPTS,'NAME1',VALUE1,...) alters an existing
%   options structure OLDOPTS.
%
%   OPTIONS = SVMSMOSET(OLDOPTS,NEWOPTS) combines an existing options
%   structure OLDOPTS with a new options structure NEWOPTS. Any new
%   properties overwrite corresponding old properties.
%
%   SVMSMOSET with no input arguments displays all property names and their
%   possible values.
%
%   SVMSMOSET has the following properties:
%
%   TolKKT
%   Tolerance with which the Karush-Kuhn-Tucker (KKT) conditions are
%   checked. Default value is 1e-3.
%
%   MaxIter
%   Maximum number of iterations of main loop. If this number is exceeded
%   before the algorithm converges then the algorithm stops and gives an
%   error. Default value is 15000.
%
%   Display
%   Controls the level of information about the optimization iterations
%   that is displayed as the algorithm runs. The value can be 'off', which
%   displays nothing, 'iter', which reports every 500 iterations, and
%   'final', which reports when the algorithm finishes. Default value is
%   'off'.
%
%   KKTViolationLevel
%   This number specifies the fraction of alphas that are allowed to
%   violate the KKT conditions. Setting this to a value greater than 0 will
%   help the algorithm to converge if it is fluctuating near a good
%   solution. Default value is 0.
%
%   KernelCacheLimit
%   This number specifies the size of the kernel matrix cache. The
%   algorithm keeps a matrix with up to KernelCacheLimit * KernelCacheLimit
%   double numbers in memory. Default value is 5000.
%
%   Examples:
%
%       opts = svmsmoset('Display','final','MaxIter',20000,...
%                                      'KernelCacheLimit',1000);
%       alt_opts = svmsmoset(opts,'Display','iter','KKTViolationLevel',.05);
%
% See also SVMCLASSIFY, SVMTRAIN.


%   References:
%
%     [1] Cristianini, N., Shawe-Taylor, J An Introduction to Support
%         Vector Machines, Cambridge University Press, Cambridge, UK. 2000.
%         http://www.support-vector.net
%     [2] J.C. Platt: A Fast Algorithm for Training  Support Vector
%         Machines,  http://research.microsoft.com/users/jplatt/smo.html
%     [3] R.-E. Fan, P.-H. Chen, and C.-J. Lin. Working Set Selection Using
%         Second Order Information for Training SVM. Journal of Machine
%         Learning Research, 6(2005), 1889-1918.
%     [4] L. Bottou and C.-J. Lin. Support Vector Machine Solvers. 2006,
%         available from http://www.csie.ntu.edu.tw/~cjlin/papers.html

%   Copyright 2006-2012 The MathWorks, Inc.


%   
%   MaxNonBoundsIter -- may get added at a later date. Currently hardcoded
%   Maximum number of iterations of the loop which tries to make the set of
%   non-bound alphas (true support vectors) consistent. If this number is
%   exceeded the algorithm continues with loop over the full set of alphas.
%   Tuning this number can speed up the algorithm. Default value is 25.

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    fprintf('            Display: [ off | iter | final ]\n');
    fprintf('             TolKKT: [ positive scalar ]\n');
    fprintf('            MaxIter: [ positive scalar ]\n');
    fprintf('   KernelCacheLimit: [ positive scalar ]\n');
    fprintf('  KKTViolationLevel: [ positive scalar]\n');
    fprintf('\n');
    return;
end

% Create a struct of all the fields with all values set to
Options = {...
    'Display', 'off';
    'TolKKT', 1e-3;
    'MaxIter', 15000;
    'KKTViolationLevel', 0;
    'KernelCacheLimit', 5000;};

Names = Options(:,1);
Defaults = Options(:,2);

m = size(Names,1);

% Combine all leading options structures o1, o2, ... in odeset(o1,o2,...).
for j = 1:m
    options.(Names{j}) = Defaults{j};
end
% work through the inputs until we find a parameter name. Handle options
% structures as we go.
i = 1;
while i <= nargin
    arg = varargin{i};
    if ischar(arg)                         % arg is an option name
        break;
    end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error(message('stats:svmtrain:NoPropNameOrStruct', i));
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),Names{j}))
                val = arg.(Names{j});
            else
                val = [];
            end
            if ~isempty(val)
                options.(Names{j}) = val;
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error(message('stats:svmtrain:ArgNameValueMismatch'));
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};

    if ~expectval
        if ~ischar(arg)
            error(message('stats:svmtrain:NoPropName', i));
        end
        k = find(strncmpi(arg, Names,numel(arg)));
        if isempty(k)
            error(message('stats:svmtrain:UnknownParameterName', arg));
        elseif length(k)>1
            error(message('stats:svmtrain:AmbiguousParameterName', arg));
        end
        expectval = 1;                      % we expect a value next

    else
        options.(Names{k}) = arg;
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error(message('stats:svmtrain:NoValueForProp', arg));
end

%check tolkkt
