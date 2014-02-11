function [c1, c2] = circularity(img)
    [r, c] = find(img == 1);
    b = bwtraceboundary(img, [r(1) c(1)], 'N', 8);
    N = size(b, 1); % pixels on permiter
    bRows = b(:,1);
    bCols = b(:,2);  
    % Calculate the circularity of the region using |P^2|/A
    A = size(r, 1);
    P = 0;
    for idx = 1:(N-1)
        diffR = bRows(idx)-bRows(idx+1);
        diffC = bCols(idx)-bCols(idx+1);
        P = P + sqrt(diffR^2 + diffC^2);
    end
    c1 = (P*P)/A;
    
    % Calculate the circularity of the region using u/o
    % Difference in each row and column values from mean
%     diffR = bRows - mean(r);
%     diffC = bCols - mean(c);
%     u = sum(sqrt(diffR.^2 + diffC.^2))/N; % mean
%     v = sum((sqrt(diffR.^2 + diffC.^2) - u).^2)/N; % variance
%     o = sqrt(v); % standard deviation
%     c2 = u/o;
    c2 = 0;
end

