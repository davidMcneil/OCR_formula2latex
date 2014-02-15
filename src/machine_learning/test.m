function [] = test()
    load('testing');
    load('net');
    S = 5;
    FM = FM(1:S:end, :);
    IM = IM(1:S:end, :);
    result = SVMClassify(net, FM);
    acc = statistics(result, IM);
    fprintf('ACC: %.4f\n', acc);
end

function [acc] = statistics(results, IM)
    acc = numel(find((results - IM) == 0))/numel(IM);
end