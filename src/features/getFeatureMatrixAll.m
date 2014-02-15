function [FM, IM, MinMax] = getFeatureMatrixAll(filepath)
    % Get the feature matrix of a folder with labeled subfolders
    subdirs = {'zero', 'one', 'two', 'three', 'four', 'five', 'six', ...
              'seven', 'eight', 'nine', 'addition', 'subtraction', ...
              'multiplication', 'division', 'equals'};
    FM = [];
    IM = [];
    for s=subdirs
       F = [filepath '/' char(s)];
       disp(F);
       [fm, im] = getFeatureMatrixFolder(F);
       FM = [FM; fm];
       IM = [IM; im];
    end
    [FM, MinMax] = normalizeDataSet(FM);
end

function [FM, IM] = getFeatureMatrixFolder(filepath)
    % Get the feature matrix of all subdirs in folder
    keyset = {'zero', 'one', 'two', 'three', 'four', 'five', 'six', ...
              'seven', 'eight', 'nine', 'addition', 'subtraction', ...
              'multiplication', 'division', 'equals'};
    valueset = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
    lookup = containers.Map(keyset, valueset);

    
    files = dir(filepath);
    files = files(3:end);
    N = numel(files); % Number of images
    NF = 26; % Number feature vectors
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
        F = [filepath, '/', files(i).name];
%         disp(F);
        FM(i, :) = getFeatureMatrix(F);
    end
end

function [MFM, MinMax] = normalizeDataSet(MFM)
    % Normalize a set of  feature matricies
    % Normalize using (data - min)/(max - min)
    % Store max and min of each column in a separate matrix
    [nr, nc] = size(MFM);
    maxc = max(MFM);
    minc = min(MFM);
    maxmat = repmat(maxc, nr, 1);
    minmat = repmat(minc, nr, 1);
    MinMax = [minc', maxc'];
    MFM = (MFM - minmat)./(maxmat - minmat);
end
