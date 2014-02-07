function fm = getFeatureMatrixFolder(filepath)
    files = dir(filepath);
    fm = zeros(numel(files));
    for i=1:numel(files)
        fm(i,:) = getFeatureMatrix(files(i));
    end
end