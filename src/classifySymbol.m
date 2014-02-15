function [ class ] = classifySymbol(img)
    % Use the SVM to classify the given image and return the class
    % img- img or filename to img
    addpath([pwd '/machine_learning'], [pwd '/features']);
    load('MinMax');
    load('net');
    % Check parameter
    if ischar(img)
       img = imread(img);
    end
    if ~islogical(img)
        img = im2bw(img);
    end
    imtool(img)
    
    fm = getFeatureMatrix(img);
    fm = normalize(fm, MinMax);
    fm
    class = SVMClassify(net, fm);
end

function [FM] = normalize(FM, MinMax)
    % Normalize a single feature vector
    MinMax = MinMax';
    maxmat = MinMax(2,:);
    minmat = MinMax(1,:);
    FM = (FM - minmat)./(maxmat - minmat);
    % Scale to boundries
    FM(FM < 0) = 0;
    FM(FM > 1) = 1; 
end
