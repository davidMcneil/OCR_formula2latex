function [symbols] = findSymbols(filename)
    V = 0;
    img = imread(filename);
    show(img, V);

    % Convert to binary image
    BW = toBinary(img);
    show(BW, V);
    
    % Extract symbols
    [BB, N] = findBoundingBoxes(BW);
    symbols = cell(N); % Used to store image symbols
    for i=1:N
        symbols{i} = imcrop(BW, BB(i,:)); % Crop and save images
        show(symbols{i}, V)
    end    
end

function show(img, V)
    if V
        imtool(img);
    end
end
    