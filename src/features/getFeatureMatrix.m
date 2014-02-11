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
   F(1) = circularity(img);
   F(2) = elongation(img);
   F(3) = S.Orientation;
   F(4) = S.EulerNumber;
   F(5) = S.Solidity;
   F(6) = S.Perimeter/BBP;
   F(7) = S.EquivDiameter/BBP;
   F(8) = S.Centroid(2) / nr;
   F(9) = S.Centroid(1) / nc;
   F(10) = std(r);
   F(11) = std(c);
   F(12) = skewness(r);
   F(13) = skewness(c);
   F(14) = kurtosis(r);
   F(15) = kurtosis(c);
   F(16) = S.Eccentricity;
   F(17) = S.Extent;
   F(18) = percentFill(img(1:fix(nr/3), 1:fix(nc/3)));
   F(19) = percentFill(img(fix(nr/3):fix(2*nr/3), 1:fix(nc/3)));
   F(20) = percentFill(img(fix(2*nr/3):nr, 1:fix(nc/3)));
   F(21) = percentFill(img(1:fix(nr/3), fix(nc/3):fix(2*nc/3)));
   F(22) = percentFill(img(fix(nr/3):fix(2*nr/3), fix(nc/3):fix(2*nc/3)));
   F(23) = percentFill(img(fix(2*nr/3):nr, fix(nc/3):fix(2*nc/3)));
   F(24) = percentFill(img(1:fix(nr/3), fix(2*nc/3):nc));
   F(25) = percentFill(img(fix(nr/3):fix(2*nr/3), fix(2*nc/3):nc));
   F(26) = percentFill(img(fix(2*nr/3):nr, fix(2*nc/3):nc));
  
   % Feature Label Vector
   L = cell(N, 1);
   L{1} = 'Circularity (P^2/A)';
%    L{2} = 'Circularity (std/mean)';
   L{2} = 'Elongation';
   L{3} = 'Orientation';
   L{4} = 'EulerNumber';
   L{5} = 'Solidity';
   L{6} = 'Perimeter/BBP';
   L{7} = 'EquivDiameter/BBP';
   L{8} = '% Row Centroid';
   L{9} = '% Col Centroid';
   L{10} = 'Row std';
   L{11} = 'Col std';
   L{12} = 'Row Skewness';
   L{13} = 'Col Skewness';
   L{14} = 'Row Kurtosis';
   L{15} = 'Col Kurtosis';
   L{16} = 'Eccentricity';  
   L{17} = 'Extent';
   L{18} = 'Extent (TL) region';
   L{19} = 'Extent (TC) region';
   L{20} = 'Extent (TR) region';
   L{21} = 'Extent (CL) region';
   L{22} = 'Extent (CC) region';
   L{23} = 'Extent (CR) region';
   L{24} = 'Extent (BL) region';
   L{25} = 'Extent (BC) region';
   L{26} = 'Extent (BR) region';

   
   % Display outbput if verbose
   if V
     horzcat(L, strread(num2str(F),'%s'))
   end
end