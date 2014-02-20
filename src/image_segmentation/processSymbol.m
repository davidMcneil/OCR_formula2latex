function [image] = processSymbol(img)
    [r, c] = size(img);
    if r > c
        padding = zeros(r, (r - c)/2);
        img = [padding img padding];
    end
    if r < c
        padding = zeros((c - r)/2, c);
        img = [padding; img; padding];
    end
    side_padding = zeros(size(img, 1), 50);
    t_b_padding = zeros(50, size(img, 2) + 100);
    img = [side_padding img side_padding];
    img = [t_b_padding; img; t_b_padding];
    % Process all symbols to a standard size to extract data from
    img_size = size(img,1)*size(img,2);
    element = strel('square',3);

    while length(find(img==1))/img_size < 0.3
        img = imdilate(img,element);
    end
%     props = regionprops(img);
%     img = imcrop(img, props.BoundingBox);
    image = imresize(img,[50,50]);
end