function [] = train()
    load('training')
    S = 1;
    FM = FM(1:S:end, :);
    IM = IM(1:S:end, :);
    net = SVMTrain(FM, IM);
    save('net.mat', 'net');
end

