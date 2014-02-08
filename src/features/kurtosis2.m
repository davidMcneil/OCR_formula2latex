function [rKurtosis, cKurtosis] = kurtosis(img)
   [r, c] = find(img == 1);
   rNorm = r - mean(r);
   cNorm = c - mean(c);
   rNormSquared = (rNorm).^2;
   cNormSquared = (cNorm).^2;
   rNormFourthed = (rNorm).^4;
   cNormFourthed = (cNorm).^4;
   rNormSum = sum(rNormFourthed);
   cNormSum = sum(cNormFourthed);
   rNormSquaredSum = sum(rNormSquared);
   cNormSquaredSum = sum(cNormSquared);
   rKurtosis = (rNormSum / size(r, 1)) / (rNormSquaredSum/size(r, 1) ^ (2));
   cKurtosis = (cNormSum / size(c, 1)) / (cNormSquaredSum/size(c, 1) ^ (2)); 
end