function [symbols] = findSymbols(img)

    V = 0;
    if ischar(img)
       img = imread(img);
    end
    show(img, V);

    % Convert to binary image
    if islogical(img);
        BW = img;
    else
        BW = toBinary(img);
    end
    show(BW, V);
    
    % Extract symbols
    [BB, N] = findBoundingBoxes(BW);
    symbols = cell(N, 1); % Used to store image symbols
    for i=1:N
        symbols{i} = imcrop(BW, BB(i,:)); % Crop and save images
        show(symbols{i}, V);
%         imwrite(symbols{i}, ['results/', int2str(i), '.png']);
    end    
end

function show(img, V)
    if V
        imtool(img);
    end
end