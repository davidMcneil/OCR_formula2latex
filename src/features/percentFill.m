function per = percentFill(img)
   per = size(find(img == 1), 1)/(size(img, 1) * size(img, 2));
end