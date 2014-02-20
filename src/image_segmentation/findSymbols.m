function [symbols] = findSymbols(img, V)
    % Extract the symbols from an image
    if ischar(img)
       img = imread(img);
    end
    if nargin < 2
        V = 0;
    end 

    % Convert to binary image
    if islogical(img);
        BW = img;
    else
        BW = toBinary(img, V);
    end
    
    % Extract symbols
    [BB, N] = findBoundingBoxes(BW);
    symbols = cell(N, 1); % Used to store image symbols
    for i=1:N
        symbols{i} = processSymbol(imcrop(BW, BB(i,:))); % Crop and save images
        show(symbols{i}, ['Symbol: ' int2str(i)], V);
%         imwrite(symbols{i}, ['results/', 'new_' int2str(i), '.png']);
    end    
end

function [BB, N] = findBoundingBoxes(bw)
    % find bounding boxes of symbols within a bw image
    [L, N] = bwlabel(bw); % L-labled regions, N-num regions
    S = regionprops(L); % S-stats on regions
    C = [S.Centroid]; % Centroids
    C = reshape([C(:,2:2:end) C(:,1:2:end)], N, 2); % Centroids
    D = squareform(pdist(C)); % Distances between centroids
    D = D - tril(D); % 0 lower half of matrix, duplicate data
    % TODO: Make thresh more robust
    thresh = numel(bw) * 0.0001; % Threshold for single symbol
    [r, c] = find(D < thresh & D ~= 0);
    nc = numel(r); % Number of changes to be made
    for i=1:nc
        n = max(r(i), c(i)); % New
        o = min(r(i), c(i)); % Old
        L(L == o) = n;
    end
    S = regionprops(L); % Get the new region props
    BB = [S.BoundingBox];
    BB = reshape([BB(:,1:4:end) BB(:,2:4:end) BB(:,3:4:end) BB(:,4:4:end)], N, 4);
    area = [S.Area]; 
    % Remove all bounding boxes that have 0 area
    idx = find(area > 0); 
    BB = BB(idx, :);
    % Change the number of regions
    N = N - nc;
end

function show(img, cap, V)
    if V
        figure('Name', cap);
        imshow(img);
    end
end