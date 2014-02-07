function [BW] = toBinary(img)
    % Convert an image to a binary black and white image
    V = 0; 
    addpath('edge_linking'); % Add edge linking functions
    show(img, V);
    
    % Convert to gray scale and apply contrast equalization
    grayimg = adapthisteq(rgb2gray(img));
    [r, c] = size(grayimg);
    show(img, V);
    
    % Use matlabs conversion to BW
    BW = ~im2bw(img, graythresh(img)); 
    show(BW, V);
    
    % Use kmeans algorithm to extract 2 colors
    K = BGvsEQ(grayimg, 2);
    show(K, V);
    
    % Combine (K & BW) inorder to flter out noise
    BW = K & BW;
    show(BW, V);
    
    % Find the edges of the image
    E = edge(grayimg, 'sobel');
    thresh = ceil(max(r, c) * .01); % Threshold to throw out an edge because it is too small  
    E = cleanedgelist(edgelink(E), thresh); % Link edges and remove short edges
    E = edgelist2image(E, [r, c]); % Convert back to image
    E = bwareaopen(E, thresh); % Remove small edges
    show(E, V);
    
    % Calculate average width of symbol
    D = bwdist(E); % Get distance to edges
    M = imregionalmax(D) .* D; % Get peaks of local maximimums
    M(M==0) = NaN;
    W = nanmedian(nanmedian(M)) * 2; % Width of average symbol
    D = D < W; % Fill in the symbols
    show(D, V)
    
    % Combine (D & BW) inorder to flter out noise
    BW = D & BW;
    show(BW, V);
    
    % Filter out symbols based on average size
    A = regionprops(bwlabel(BW), 'Area'); % Get areas of all regions
    thresh = ceil(mean(mean([A.Area])) * .2); % Threshhold to remove 1/10 of mean
    BW = bwareaopen(BW, thresh);
end

function show(img, V)
    if V
        imtool(img);
    end
end