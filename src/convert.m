function res = convert(img, V)
    addpath('image_segmentation');
    if nargin < 2
        V = 0;
    end 
    if ischar(img)
       img = imread(img);
    end
    if V 
        figure('Name', 'Oringinal Image');
        imshow(img, 'border', 'tight');
    end
    
    keySet = {0, 1, 2, 3, 4, 5,  6, 7, 8, 9, 10, 11, 12, 13, 14};
    valueSet = {'0', '1', '2', '3', '4', '5',  '6', '7', '8', '9', '+', '-', '*', '/', '='};
    L = containers.Map(keySet,valueSet);
    
    symbols = findSymbols(img, V);
    res = '';
    for i=1:numel(symbols)
        res = [res L(classifySymbol(symbols{i}))];
    end
end