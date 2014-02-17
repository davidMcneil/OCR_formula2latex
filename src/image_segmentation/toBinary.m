function [BW] = toBinary(img)
    % Convert an image to a binary black and white image
    V = 0; 
    show(img, V);
    
    % Convert to gray scale and apply contrast equalization
    grayimg = adapthisteq(rgb2gray(img));
    [r, c] = size(grayimg);
    show(img, V);
    
    % Use matlabs conversion to BW
    BW = ~im2bw(img, graythresh(img)); 
    show(BW, 1);
    
    % Use kmeans algorithm to extract 2 colors
    K = BGvsEQ(grayimg, 2);
    show(K, 1);
    
    % Find the edges of the image
    E = edge(grayimg, 'sobel');
    thresh = ceil(max(r, c) * .02); % Threshold to throw out an edge because it is too small  
    E = bwareaopen(E, thresh); % Remove small edges
    show(E, 1);
    
    % Calculate average width of symbol
    D = bwdist(E); % Get distance to edges
    show(uint8(D), 1);
    M = imregionalmax(D) .* D; % Get peaks of local maximimums
    show(M, 1);
    M(M==0) = NaN;
    W = nanmedian(nanmedian(M)) * 2; % Width of average symbol
    D = D < W; % Fill in the symbols that are within the width
    show(D, V);
    D = imerode(D, ones(W)); % Erode the extra 
    show(D, 1);
    
    % Combine (K & D & BW) inorder to flter out noise
    BW = K & D & BW;
    
    % Filter out symbols based on average size
    A = regionprops(bwlabel(BW), 'Area'); % Get areas of all regions
    thresh = ceil(mean(mean([A.Area])) * .2); % Threshhold to remove 1/10 of mean
    BW = bwareaopen(BW, thresh);
%     show(BW, V);
end

function [out] = BGvsEQ(grayimg, K)
    % Extract the EQuation from the BackGround
    % Find "K" clusters in "grayimg" and return "out" as 
    % resulting segmentation with regions labeled 0-(K-1) 
    % starting with lighter to darker regions
    % Note: none deterministic, may not find all "K"
    grayimg = double(grayimg);
    [r, c, d] = size(grayimg);
    grayimg = reshape(grayimg, r*c, d); % Resize for use by kmeans
    idx = kmeans(grayimg, K, ...  % Compute K means
        'start', 'uniform', ...
        'emptyaction', 'singleton');
    % Resize back to image format, -1 to start regions at 0
    out = reshape(idx, r, c) - 1;
    idxD = find(out == K - 2); % Current dark region
    idxL = find(out == K - 1); % Current light region
    if (grayimg(idxD(1)) < grayimg(idxL(1))) % If the dark is actually lighter
        out = ~out; % Flip the 0's and 1's
    end
end

function show(img, V)
    if V
        imtool(img);
    end
end