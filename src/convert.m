function res = convert(img)
    addpath('image_segmentation');
    imtool(img);
    if ischar(img)
       img = imread(img);
    end
    
    keySet = {0, 1, 2, 3, 4, 5,  6, 7, 8, 9, 10, 11, 12, 13, 14};
    valueSet = {'0', '1', '2', '3', '4', '5',  '6', '7', '8', '9', '+', '-', '*', '/', '='};
    L = containers.Map(keySet,valueSet);
    
    symbols = findSymbols(img);
    
    res = '';
    for i=1:numel(symbols)
        imtool(symbols{i});
        res = [res L(classifySymbol(symbols{i}))];
    end
end