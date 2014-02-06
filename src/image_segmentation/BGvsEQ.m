function [out] = BGvsEQ(grayimg)
    % Extract the EQuation from the BackGround
    % Find "K" clusters in "grayimg" and return "out" as 
    % resulting segmentation with regions labeled 0-(K-1) 
    % starting with lighter to darker regions
    % Note: none deterministic, may not find all "K"
    K = 2; % Always looking for background vs equation, 2 regions
    grayimg = double(grayimg);
    [r, c, d] = size(grayimg);
    grayimg = reshape(grayimg, r*c, d); % Resize for use by kmeans
    idx = kmeans(grayimg, K, ...  % Compute K means
        'start', 'uniform', ...
        'emptyaction', 'singleton');
    % Resize back to image format, -1 to start regions at 0
    out = reshape(idx, r, c) - 1;
    idxD = find(out == K - 2); % Current dark region
    idxL = find(out == K - 1); % Current light region
    if (grayimg(idxD(1)) < grayimg(idxL(1))) % If the dark is actually lighter
        out = ~out; % Flip the 0's and 1's
    end
end

