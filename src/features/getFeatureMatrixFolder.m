function [fm, im] = getFeatureMatrixFolder(filepath)
    lookup = ['zero', 'one', 'two', 'three', 'four', 'five', 'six', ...
        'seven', 'eight', 'nine', 'addition', 'subtraction', ...
        'multiplication', 'division', 'equals'; ...
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
    files = dir(filepath);
    fm = zeros(numel(files));
    k = findstr(filepath, lookup(1, :));
    num = lookup(2, k);
    for i=1:numel(files)
        fm(i,:) = getFeatureMatrix(files(i));
        im(i,1) = num;
    end
end