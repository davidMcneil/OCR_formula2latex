function [cSkew, rSkew] = skew(img)
    [r, c] = find(img == 1);
    rSkew = skewness(r);
    cSkew = skewness(c);
%     rNorm = r - mean(r);
%     cNorm = c - mean(c);
%     rNormSquared = (rNorm).^2;
%     cNormSquared = (cNorm).^2;
%     rNormCubed = (rNorm).^3;
%     cNormCubed = (cNorm).^3;
%     
%     rSkew = mean(rNormCubed) / (mean(rNormSquared) ^ (3/2));
%     cSkew = mean(cNormCubed) / (mean(cNormSquared) ^ (3/2));
end