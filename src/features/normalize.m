function [FM, minmax] = normalize(FM)
    % Normalize a feature matrix
    % Normalize using (data - min)/(max - min)
    % Save max and min of each column in a separate matrix
    [nr, nc] = size(FM);
    maxc = max(FM);
    minc = min(FM);
    maxmat = repmat(maxc, nr, 1);
    minmat = repmat(minc, nr, 1);
    minmax = [minc', maxc'];
    FM = (FM - minmat)./(maxmat - minmat);
end

