function [F, L] = getFeatureMatrix(img, V)
   % Operates on a segment of a bw image
   % img- can be either a filename to a binary image or a binary image
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
   N = 26; % Number of features
   
   % Feature Vector
   F(1) = circularity(img);
   F(2) = elongation(img);
   F(3) = S.Orientation;
   F(4) = S.EulerNumber;
   F(5) = S.Solidity;
   F(6) = S.Perimeter / BBP;
   F(7) = S.EquivDiameter / BBP;
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
   F(18) = percentFill(img(1:fix(nr/3), 1:fix(nc/3)), S.Area);
   F(19) = percentFill(img(fix(nr/3):fix(2*nr/3), 1:fix(nc/3)), S.Area);
   F(20) = percentFill(img(fix(2*nr/3):nr, 1:fix(nc/3)), S.Area);
   F(21) = percentFill(img(1:fix(nr/3), fix(nc/3):fix(2*nc/3)), S.Area);
   F(22) = percentFill(img(fix(nr/3):fix(2*nr/3), fix(nc/3):fix(2*nc/3)), S.Area);
   F(23) = percentFill(img(fix(2*nr/3):nr, fix(nc/3):fix(2*nc/3)), S.Area);
   F(24) = percentFill(img(1:fix(nr/3), fix(2*nc/3):nc), S.Area);
   F(25) = percentFill(img(fix(nr/3):fix(2*nr/3), fix(2*nc/3):nc), S.Area);
   F(26) = percentFill(img(fix(2*nr/3):nr, fix(2*nc/3):nc), S.Area);
  
   % Feature Label Vector
   L = cell(N, 1);
   L{1} = 'Circularity (P^2/A)';
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


function elong = elongation(img)
    [r, c] = find(img == 1);
    cMean = mean(c);
    cNorm = c - cMean;
    cNormSquared = (cNorm).^2;
    cSumSquared = sum(cNormSquared);
    s = size(r);
    N = s(1);
    cm(1,1) = cSumSquared / N;
    rMean = mean(r);
    rNorm = r - rMean;
    rNormSquared = (rNorm).^2;
    rSumSquared = sum(rNormSquared);
    cm(2,2) = rSumSquared / N;
    crNormSquared = cNorm .* rNorm;
    %crNormSquared = (crNorm).^2;
    crSumSquared = sum(crNormSquared);
    cm(2,1) = crSumSquared / N;
    cm(1,2) = cm(2,1);
    evals = eig(cm);
    if evals(1) > evals(2)
        elong = sqrt(evals(1)/evals(2));
    else
        elong = sqrt(evals(2)/evals(1));
    end
end

function c1 = circularity(img)
    % Calculate the circularity of the region using |P^2|/A
    [r, c] = find(img == 1);
    b = bwtraceboundary(img, [r(1) c(1)], 'N', 8);
    N = size(b, 1); % pixels on permiter
    bRows = b(:,1);
    bCols = b(:,2);    
    A = size(r, 1);
    P = 0;
    for idx = 1:(N-1)
        diffR = bRows(idx)-bRows(idx+1);
        diffC = bCols(idx)-bCols(idx+1);
        P = P + sqrt(diffR^2 + diffC^2);
    end
    c1 = (P*P)/A;
end

function per = percentFill(img, area)
   per = size(find(img == 1), 1) / area;
end