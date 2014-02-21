function [] = getTrainingAndTesting()
    addpath('../../src/image_segmentation');
    subdirs = {'zero', 'one', 'two', 'three', 'four', 'five', 'six', ...
              'seven', 'eight', 'nine', 'addition', 'subtraction', ...
              'multiplication', 'division', 'equals'};

    C = 0;
    for s=subdirs % All symbols
       D = char(s);
       files = dir(D); files = files(3:end);
       N = numel(files); % Number of images
       mkdir(['training/' D]);
       mkdir(['testing/' D]);
       for i=1:N % All files in symbols
            n = files(i).name;
            F = [D, '/', n];
            I = imreadAutoRot(F);
%             imtool(I);
            symb = findSymbols(I);
            for y=1:numel(symb)
%                 C
                if mod(C, 5) == 0
                    imwrite(symb{y}, ['testing/' D '/' int2str(C) '.png']);
                else
                    imwrite(symb{y}, ['training/' D '/' int2str(C) '.png']);
                end
                C = C + 1;
            end
        end
    end
end

function im = imreadAutoRot(filename)

    % This is a replacement of imread in Matlab to handle auto-rotation in JPEG.
    % Matlab seems not able to handle automatic rotation of image in imread 
    % (at least until R2012a version).
    % Therefore, I wrote this file to automatically rotate the image into the 
    % write direction based on EXIF orientation. 
    % I have tested this function on iPhone 5 with iOS 6.
    % Jianxiong Xiao
    % Reference: JPEG image format at http://www.impulseadventure.com/photo/exif-orientation.html

    im = imread(filename);

    try
        info = imfinfo(filename);
        switch info.Orientation
            case 1

            case 2
                im = im(:,end:-1:1,:);
            case 3
                im = imrotate(im,180);
                im = im(:,end:-1:1,:);
            case 4
                im = im(:,end:-1:1,:);

            case 6
                im = imrotate(im,-90);            
            case 5
                im = im(:,end:-1:1,:);
                im = imrotate(im,-90);            

            case 8
                im = imrotate(im,90);            
            case 7
                im = imrotate(im,90);            
                im = im(:,end:-1:1,:);
        end
    catch
    end
end