function [out] = BGvsEQ(Gimg)
    % Extract the EQuations from the BackGround
    % Find "K" clusters in "Gimg (G=grayscale)" and return "out" as 
    % resulting segmentation with regions labeled ~ 0-(K-1) 
    % starting with lighter to darker regions
    % Note: none deterministic, may not find all "K"
    K = 2; % Always looking for background vs equation, 2 regions
    Gimg = double(Gimg);
    [r, c, d] = size(Gimg);
    Gimg = reshape(Gimg, r*c, d); % Resize for use by kmeans
    idx = kmeans(Gimg, K); % Compute K means
    % Resize back to image format, -1 to start regions at 0
    out = reshape(idx, r, c) - 1;
    idxD = find(out == K - 2); % Current dark region
    idxL = find(out == K - 1); % Current light region
    if (Gimg(idxD(1)) < Gimg(idxL(1))) % If the dark is actually lighter
        out = ~out; % Flip the 0's and 1's
    end
end

