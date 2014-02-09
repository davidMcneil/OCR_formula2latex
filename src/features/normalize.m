function [FM] = normalize(FM)
    % Normalize a feature matrix
    [nr, nc] = size(FM);
    N = [0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1];
    N = repmat(N, nr, 1);
    FM = FM .* N;
end

