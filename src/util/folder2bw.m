filepath = '../machine_learning/test';
files = dir(filepath);
files = files(3:end);
N = numel(files); % Number of images

for i=1:N
    img = imread([filepath '/' files(i).name]);
    img = im2bw(img);
%     imtool(img);
    imwrite(img, [filepath '/' files(i).name]);
end