function res = convert(img)
    addpath('image_segmentation');
    imtool(img);
    if ischar(img)
       img = imread(img);
    end
    
    symbols = findSymbols(img);
    
    res = '';
    for i=1:numel(symbols)
        imtool(symbols{i});
        res = [res ', ' int2str(classifySymbol(symbols{i}))];
    end
    
end