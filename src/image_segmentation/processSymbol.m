function [image] = processSymbol(img)

    img_size = size(img,1)*size(img,2);
    element = strel('square',3);

    while length(find(img==1))/img_size < 0.5

        img = imdilate(img,element);

    end

    image = imresize(img,[23,23]);
end