function [symbols] = processImage(filename)
    img = imread(filename);
%     imtool(img);

    bw = toBinary(img);
%     imtool(bw);
    
    % Extract symbols
    [L, N] = bwlabel(bw); % L-labled regions, N-num regions
    S = regionprops(bw); % S-stats on regions
    symbols = cell(N); % Used to store symbols in image
    for i=1:N
        symbols{i} = imcrop(L, S(i).BoundingBox);
%         imtool(symbols{i});
    end    
end

