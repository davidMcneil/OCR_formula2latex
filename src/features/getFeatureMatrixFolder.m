function [FM, IM] = getFeatureMatrixFolder(filepath)
    keyset = {'zero', 'one', 'two', 'three', 'four', 'five', 'six', ...
              'seven', 'eight', 'nine', 'addition', 'subtraction', ...
              'multiplication', 'division', 'equals'};
    valueset = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
    lookup = containers.Map(keyset, valueset);

    
    files = dir(filepath);
    files = files(3:end);
    N = numel(files); % Number of images
    NF = 27; % Number feature vectors
    FM = zeros(numel(files), NF); % Feature matrix
    IM = ones(numel(files), 1) .* -1; % Identity Matrix
    
    for k=lookup.keys
        k = char(k);
        if numel(findstr(filepath, k)) > 0
            IM = IM .* -lookup(k);
            break;
        end
    end
    
    for i=1:N
        FM(i, :) = getFeatureMatrix([filepath, '/', files(i).name]);
    end
end