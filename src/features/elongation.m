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