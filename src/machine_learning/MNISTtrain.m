function [] = train()
    load('MNISTtraining')
    S = 1;
    FM = FM(1:S:end, :);
    IM = IM(1:S:end, :);
    net = SVMTrain(FM, IM);
    save('MNISTnet.mat', 'net');
end