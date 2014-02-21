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
    
%     disp('MNISTTRAINING');
%     [FM, IM, MinMax] = getFeatureMatrixAll([path '/MNISTtraining']);
%     save('MNISTtraining.mat', 'FM', 'IM');
%     save('MNISTMinMax.mat', 'MinMax');
%     
%     disp('MNISTTESTING');
%     [FM, IM, MinMax] = getFeatureMatrixAll([path '/MNISTtesting']);
%     save('MNISTtesting.mat', 'FM', 'IM');
end
