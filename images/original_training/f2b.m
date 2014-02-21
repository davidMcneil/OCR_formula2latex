subdirs = {'zero', 'one', 'two', 'three', 'four', 'five', 'six', ...
              'seven', 'eight', 'nine', 'addition', 'subtraction', ...
              'multiplication', 'division', 'equals'};
for s=subdirs
filepath = char(s);
files = dir(filepath);
files = files(3:end);
N = numel(files); % Number of images

for i=1:N
    img = imread([filepath '/' files(i).name]);
    img = im2bw(img);
% imtool(img);
    imwrite(img, [filepath '/' files(i).name]);
end
end