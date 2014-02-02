function [out] = findEdge(Gimg)
    % Find the edges of the "Gimg"
    horizMask = (1/8)*[-1, -2, -1; 0, 0, 0; 1, 2, 1];
    vertMask = (1/8)*[-1, 0, 1; -2, 0, 2; -1, 0, 1];
    sobelHoriz = filter2(horizMask, Gimg);
    sobelVert = filter2(vertMask, Gimg);
    out = uint8(sqrt(sobelHoriz .^ 2 + sobelVert .^ 2)); % Gradient
    c = 25; % Cutoff for edge
    out(find(out < c)) = 0;
    out(find(out >= c)) = 1;
    out = double(out);
end

