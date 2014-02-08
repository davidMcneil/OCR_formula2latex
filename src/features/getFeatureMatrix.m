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
   [nr, nc] = size(img);
   [r, c] = find(img == 1);
   BBP = (2*nr + 2*nc); % Bounding Box Perimeter
   N = 23; % Number of features
   
   % Feature Vector
   [F(1), F(2)] = circularity(img);
   F(3) = elongation(img);
   F(4) = S.Orientation;
   F(5) = S.EulerNumber;
   F(6) = S.Solidity;
   F(7) = S.Perimeter/BBP;
   F(8) = S.EquivDiameter/BBP;
   F(9) = S.Centroid(2) / nr;
   F(10) = S.Centroid(1) / nc;
   F(11) = var(r);
   F(12) = var(c);
   F(13) = skewness(r);
   F(14) = skewness(c);
   F(15) = kurtosis(r);
   F(16) = kurtosis(c);
   F(17) = S.Eccentricity;
   F(18) = S.Extent;
   F(19) = percentFill(img(1:fix(nr/3), 1:fix(nc/3)));
   F(20) = percentFill(img(fix(nr/3):fix(2*nr/3), 1:fix(nc/3)));
   F(21) = percentFill(img(fix(2*nr/3):nr, 1:fix(nc/3)));
   F(22) = percentFill(img(1:fix(nr/3), fix(nc/3):fix(2*nc/3)));
   F(23) = percentFill(img(fix(nr/3):fix(2*nr/3), fix(nc/3):fix(2*nc/3)));
   F(24) = percentFill(img(fix(2*nr/3):nr, fix(nc/3):fix(2*nc/3)));
   F(25) = percentFill(img(1:fix(nr/3), fix(2*nc/3):nc));
   F(26) = percentFill(img(fix(nr/3):fix(2*nr/3), fix(2*nc/3):nc));
   F(27) = percentFill(img(fix(2*nr/3):nr, fix(2*nc/3):nc));

   
   % Feature Label Vector
   L = cell(N, 1);
   L{1} = 'Circularity (P^2/A)';
   L{2} = 'Circularity (std/mean)';
   L{3} = 'Elongation';
   L{4} = 'Orientation';
   L{5} = 'EulerNumber';
   L{6} = 'Solidity';
   L{7} = 'Perimeter/BondingBoxPerimeter';
   L{8} = 'EquivDiameter/BoundingBoxPerimeter';
   L{9} = '% Row Centroid';
   L{10} = '% Col Centroid';
   L{11} = 'Row variance';
   L{12} = 'Col variance';
   L{13} = 'Row Skewness';
   L{14} = 'Col Skewness';
   L{15} = 'Row Kurtosis';
   L{16} = 'Col Kurtosis';
   L{17} = 'Eccentricity';  
   L{18} = 'Extent';
   L{19} = 'Extent (TL) region';
   L{20} = 'Extent (TC) region';
   L{21} = 'Extent (TR) region';
   L{22} = 'Extent (CL) region';
   L{23} = 'Extent (CC) region';
   L{24} = 'Extent (CR) region';
   L{25} = 'Extent (BL) region';
   L{26} = 'Extent (BC) region';
   L{27} = 'Extent (BR) region';

   
   % Display outbput if verbose
   if V
     horzcat(L, strread(num2str(F),'%s'))
   end
end