function [image] = processSymbol(img)
    % Process all symbols to a standard size to extract data from
    img_size = size(img,1)*size(img,2);
    element = strel('square',3);

    while length(find(img==1))/img_size < 0.3
        img = imdilate(img,element);
    end
    image = imresize(img,[50,50]);
end