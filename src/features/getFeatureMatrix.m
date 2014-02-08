function [F, L] = getFeatureMatrix(img, V)
   % Operates on a segment of a bw image
   % img- can be eiter a filename to a binary image or a binary image
   % V- verbose
   % F- Feature vector
   % L- Corresponding featrue names
   
   % Check parameters
   if nargin < 2
    V = 0;
   end 
   if ischar(img)
       img = imread(img);
   end
   
   S = regionprops(uint8(img), 'All');
   [r, c] = size(img);
   BBP = (2*r + 2*c); % Bounding Box Perimeter
   N = 23; % Number of features
   
   % Feature Vector
   [F(1), F(2)] = circularity(img);
   F(3) = skewness(skewness(img));
   F(4) = kurtosis(kurtosis(img));
   F(5) = elongation(img);
   F(6) = S.Eccentricity;
   F(7) = S.EquivDiameter/BBP;
   F(8) = S.EulerNumber;
   F(9) = S.Orientation;
   F(10) = S.Solidity;
   F(11) = S.Extent;
   F(12) = percentFill(img(1:fix(r/3), 1:fix(c/3)));
   F(13) = percentFill(img(fix(r/3):fix(2*r/3), 1:fix(c/3)));
   F(14) = percentFill(img(fix(2*r/3):r, 1:fix(c/3)));
   F(15) = percentFill(img(1:fix(r/3), fix(c/3):fix(2*c/3)));
   F(16) = percentFill(img(fix(r/3):fix(2*r/3), fix(c/3):fix(2*c/3)));
   F(17) = percentFill(img(fix(2*r/3):r, fix(c/3):fix(2*c/3)));
   F(18) = percentFill(img(1:fix(r/3), fix(2*c/3):c));
   F(19) = percentFill(img(fix(r/3):fix(2*r/3), fix(2*c/3):c));
   F(20) = percentFill(img(fix(2*r/3):r, fix(2*c/3):c));
   F(21) = S.Perimeter/BBP;
   F(22) = S.Centroid(2) / r;
   F(23) = S.Centroid(1) / c;
   
   % Feature Label Vector
   L = cell(N, 1);
   L{1} = 'Circularity by P^2/A';
   L{2} = 'Circularity by sigma/mu';
   L{3} = 'Skewness';
   L{4} = 'Kurtosis';
   L{5} = 'Elongation';
   L{6} = 'Eccentricity';
   L{7} = 'EquivDiameter/BoundingBoxPerimeter';
   L{8} = 'EulerNumber';
   L{9} = 'Orientation';
   L{10} = 'Solidity';
   L{11} = 'Extent';
   L{12} = 'Extent (TL) region';
   L{13} = 'Extent (TC) region';
   L{14} = 'Extent (TR) region';
   L{15} = 'Extent (CL) region';
   L{16} = 'Extent (CC) region';
   L{17} = 'Extent (CR) region';
   L{18} = 'Extent (BL) region';
   L{19} = 'Extent (BC) region';
   L{20} = 'Extent (BR) region';
   L{21} = 'Perimeter/Perimeter of BoundingBox';
   L{22} = '% Row Centroid';
   L{23} = '% Col Centroid';
   
   % Display outbput if verbose
   if V
     horzcat(L, strread(num2str(F),'%s'))
   end
end