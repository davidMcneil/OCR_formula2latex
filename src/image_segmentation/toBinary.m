function [BW] = toBinary(img)
    grayimg = adapthisteq(rgb2gray(img)); % Apply contrast equalization
%     imtool(grayimg);
    
    BW = ~im2bw(img, graythresh(img));
%     imtool(BW);
    
    C = BGvsEQ(grayimg);
% %     imtool(C);
    
    BW = C & BW;
    imtool(BW);
    
    E = findEdge(grayimg);
%     imtool(E);
end