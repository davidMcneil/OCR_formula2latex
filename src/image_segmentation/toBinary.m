function [BW] = toBinary(img)
    % Convert an image to a binary black and white image
    V = 0; 
    addpath('edge_linking'); % Add edge linking functions
    show(img, V);
    
    % Convert to gray scale and apply contrast equalization
    grayimg = adapthisteq(rgb2gray(img));
    show(img, V);
    
    % Use matlabs conversion to BW
    BW = ~im2bw(img, graythresh(img)); 
    show(BW, V);
    
    % Use kmeans algorithm to extract 2 colors
    K = BGvsEQ(grayimg);
    show(K, V);
    
    % Combine (K & BW) inorder to flter out noise
    BW = K & BW;
    show(BW, V);
    
    % Find the edges of the image
    E = edge(grayimg, 'sobel');
    show(E, V);
end

function show(img, V)
    if V
        imtool(img);
    end
end
    