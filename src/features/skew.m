function [cSkew, rSkew] = skew(img)
    [r, c] = find(img == 1);
    rNorm = r - mean(r);
    cNorm = c - mean(c);
    rNormSquared = (rNorm).^2;
    cNormSquared = (cNorm).^2;
    rNormCubed = (rNorm).^3;
    cNormCubed = (cNorm).^3;
    rNormSum = sum(rNormCubed);
    cNormSum = sum(cNormCubed);
    rNormSquaredSum = sum(rNormSquared);
    cNormSquaredSum = sum(cNormSquared);
    rSkew = (rNormSum / size(r, 1)) / (rNormSquaredSum/size(r, 1) ^ (3/2));
    cSkew = (cNormSum / size(c, 1)) / (cNormSquaredSum/size(c, 1) ^ (3/2));
end