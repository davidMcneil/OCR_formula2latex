function [fm, im] = getFeatureMatrixFolder(filepath)
    filepath
    keyset = {'zero', 'one', 'two', 'three', 'four', 'five', 'six', ...
                'seven', 'eight', 'nine', 'addition', 'subtraction', ...
                'multiplication', 'division', 'equals'};
    valueset = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
    lookup = containers.Map(keyset, valueset)
    files = dir(filepath);
    fm = zeros(numel(files));
    for i=1:numel(keyset)
        key = keyset(i)
%         str = findstr(filepath, key); % This doesn't work for some reason
        str = findstr(filepath, 'four');
        if str
            k = 'four';
            k
        end
    end
    
    num = lookup(k);
    for i=3:numel(files)
        files(i)
        img = imread([filepath, '\', files(i).name]);
        fm(i, :) = getFeatureMatrix(img); % neither does this
        im(i, 1) = num;
    end
end