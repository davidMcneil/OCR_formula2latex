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
    
%     imtool(uint8(bwdist(E)));
    % Connect edges
%     imtool(E);
%     [edgelist, labelededgeim] = edgelink(E);
%     imtool(labelededgeim);
% %     drawedgelist(edgelist, size(img), 1, 'rand', 2); axis off
%     
%     % Fit line segments to the edgelists
%     tol = 2;         % Line segments are fitted with maximum deviation from
% 		     % original edge of 2 pixels.
%     seglist = lineseg(edgelist, tol);
% 
%     % Draw the fitted line segments stored in seglist in figure window 3 with
%     % a linewidth of 2 and random colours
%     drawedgelist(seglist, size(img), 2, 'rand', 3); axis off
end

function show(img, V)
    if V
        imtool(img);
    end
end