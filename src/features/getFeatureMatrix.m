function features = getFeatureMatrix(imgSeg)
   % Operates on a segment of a bw image
   % By offset:
   % Circularity by P^2/A
   % Circularity by sigma/mu
   % Row skew
   % Column skew
   % Row kurtosis
   % Column kurtosis
   % Elongation
   % Row size
   % Column size
   % Percent fill of region
   % Regionprops
   if ischar(imgSeg)
       imgSeg = imread(imgSeg);
   end
   stats = regionprops(imgSeg);
   stats
   features(1:2) = circularity(imgSeg);
   features(3:4) = skew(imgSeg);
   features(5:6) = kurtosis(imgSeg);
   features(7) = elongation(imgSeg);
%    features(8:9) = centroid(imgSeg);
   [r, c] = size(imgSeg);
   features(8:9) = [r, c];
   features(10) = percentFill(imgSeg(1:fix(r/3), 1:fix(c/3)));
   features(11) = percentFill(imgSeg(fix(r/3):fix(2*r/3), 1:fix(c/3)));
   features(12) = percentFill(imgSeg(fix(2*r/3):r, 1:fix(c/3)));
   features(13) = percentFill(imgSeg(1:fix(r/3), fix(c/3):fix(2*c/3)));
   features(14) = percentFill(imgSeg(fix(r/3):fix(2*r/3), fix(c/3):fix(2*c/3)));
   features(15) = percentFill(imgSeg(fix(2*r/3):r, fix(c/3):fix(2*c/3)));
   features(16) = percentFill(imgSeg(1:fix(r/3), fix(2*c/3):c));
   features(17) = percentFill(imgSeg(fix(r/3):fix(2*r/3), fix(2*c/3):c));
   features(18) = percentFill(imgSeg(fix(2*r/3):r, fix(2*c/3):c));
   features = [features];
end