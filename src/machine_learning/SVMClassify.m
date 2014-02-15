function [result] = SVMClassify(net, FMT)
    % Classify a testing set of data
    addpath([fileparts(mfilename('fullpath')) '/svm']);
    N = size(FMT, 1);
    NC = numel(net); % Number of SVM classes
    
    result = ones(N,1) * -1;
    %classify test cases
    res = zeros(NC, 1); % Results for one row
    for j=1:N
        for k=1:NC
            [o, res(k)] = svmclassify(net(k),FMT(j,:));
        end
%         res
        result(j) = find(res == max(res)) - 1; % Scale back to zero 
    end
end

