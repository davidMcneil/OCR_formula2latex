function [] = getTraingAndTestingData()
    % Utility function to retrive data, saved in mat file
    path = '../../images';
    disp('TRAINING');
    [FM, IM, MinMax] = getFeatureMatrixAll([path '/training']);
    save('training.mat', 'FM', 'IM');
    save('MinMax.mat', 'MinMax');
    
    disp('TESTING');
    [FM, IM, MinMax] = getFeatureMatrixAll([path '/testing']);
    save('testing.mat', 'FM', 'IM');
end
