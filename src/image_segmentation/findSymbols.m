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
%         imtool(symbols{i});
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

function [image] = processSymbol(img)
    [r, c] = size(img);
    if r > c
        padding = zeros(r, round((r - c)/2));
        img = [padding img padding];
    end
    if r < c
        padding = zeros(round((c - r)/2), c);
        img = [padding; img; padding];
    end
    side_padding = zeros(size(img, 1), 50);
    t_b_padding = zeros(50, size(img, 2) + 100);
    img = [side_padding img side_padding];
    img = [t_b_padding; img; t_b_padding];
%     Process all symbols to a standard size to extract data from
    img_size = size(img,1)*size(img,2);
    element = strel('square',3);

    while length(find(img==1))/img_size < 0.06
        img = imdilate(img,element);
    end
    image = img;
    props = regionprops(img);
    img = imcrop(img, props.BoundingBox);
    [r, c] = size(img);
    if r > c
        padding = zeros(r, round((r - c)/2));
        img = [padding img padding];
    end
    if r < c
        padding = zeros(round((c - r)/2), c);
        img = [padding; img; padding];
    end
    image = logical(image);
    image = logical(imresize(img,[50,50]));
end

function show(img, cap, V)
    if V
        figure('Name', cap);
        imshow(img, 'border', 'tight');
    end
end