function features = getFeatureMatrix(imgSeg)
   % Operates on a segment of a bw image %
   features(1:2) = circularity(imgSeg);
   features(3:4) = skew(imgSeg);
   features(5:6) = kurtosis(imgSeg);
   features(7) = elongation(imgSeg);
   features(8:9) = centroid(imgSeg);
   features(10) = percentFill(imgSeg);
   features(11:12) = size(imgSeg);
end