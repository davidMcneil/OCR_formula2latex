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
   % Percent fill of region
   % Row size
   % Column size
   % Regionprops
   stats = regionprops(imgSeg);
   features(1:2) = circularity(imgSeg);
   features(3:4) = skew(imgSeg);
   features(5:6) = kurtosis(imgSeg);
   features(7) = elongation(imgSeg);
%    features(8:9) = centroid(imgSeg);
   features(8) = percentFill(imgSeg);
   features(9:10) = size(imgSeg);
   features = [features, stats];
end