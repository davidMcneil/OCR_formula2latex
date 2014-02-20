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
            imtool(F);
            symb = findSymbols(F);
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

