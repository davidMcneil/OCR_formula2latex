function [rMean, cMean] = centroid(img)
    [r, c] = find(img == 1);
    rMean = mean(r);
    cMean = mean(c);
end