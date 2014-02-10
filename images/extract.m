clear;
addpath('../src/image_segmentation');
W = 28;
H = 28;

% load('10kTrain.mat');
% fea = full(fea);
% path = 'training/';
% 
% for i=1:numel(fea)
%     switch gnd(i)
%     case 0
%         dir = 'zero/';
%     case 1
%         dir = 'one/';
%     case 2
%         dir = 'two/';
%     case 3
%         dir = 'three/';
%     case 4
%         dir = 'four/';
%     case 5
%         dir = 'five/';
%     case 6
%         dir = 'six/';
%     case 7
%         dir = 'seven/';
%     case 8
%         dir = 'eight/';
%     case 9
%         dir = 'nine/'; 
%     end
%     syms = findSymbols(im2bw(reshape(fea(i,:),[H,W])'));
%     % Find the most predominent (largest) symbol and save it
%     [maxsize, maxidx]= max(cellfun(@numel,syms));
%     imwrite(syms{maxidx}, [path dir int2str(i) '.png']);
% end

load('10kTest.mat');
fea = full(fea);
path = 'testing/';

for i=1:numel(fea)
    switch gnd(i)
    case 0
        dir = 'zero/';
    case 1
        dir = 'one/';
    case 2
        dir = 'two/';
    case 3
        dir = 'three/';
    case 4
        dir = 'four/';
    case 5
        dir = 'five/';
    case 6
        dir = 'six/';
    case 7
        dir = 'seven/';
    case 8
        dir = 'eight/';
    case 9
        dir = 'nine/'; 
    end
    syms = findSymbols(im2bw(reshape(fea(i,:),[H,W])'));
    % Find the most predominent (largest) symbol and save it
    [maxsize, maxidx]= max(cellfun(@numel,syms));
    imwrite(syms{maxidx}, [path dir int2str(i) '.png']);
end
