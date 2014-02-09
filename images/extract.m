clear;
load('2k2k.mat');
addpath('../src/image_segmentation');

W = 28;
H = 28;
numPerLine = 10;

path = 'training/';

for i=1:numel(trainIdx)
    
    switch gnd(trainIdx(i))
    case 0
        dir = 'zero/';
    case 1
        dir = 'one/';cd
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
    syms = findSymbols(im2bw(reshape(fea(trainIdx(i),:),[H,W])'));
    % Find the most predominent (largest) symbol and save it
    [maxsize, maxidx]= max(cellfun(@numel,syms));
    imwrite(syms{maxidx}, [path dir int2str(i) '.png']);
end
