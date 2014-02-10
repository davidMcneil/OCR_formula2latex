function [FM, IM] = getFeatureMatrixAll(filepath)
    subdirs = {'zero', 'one', 'two', 'three', 'four', 'five', 'six', ...
              'seven', 'eight', 'nine', 'addition', 'subtraction', ...
              'multiplication', 'division', 'equals'};
    FM = [];
    IM = [];
    for s=subdirs
       F = [filepath '/' char(s)];
       disp(F);
       [fm, im] = getFeatureMatrixFolder(F);
       FM = [FM; fm];
       IM = [IM; im];
    end
    FM = normalize(FM);
%     save('testingFeatures.mat', 'FM', 'IM');
end