load('result.mat');
load('testingFeatures.mat');

n = size(result,1);
right = 0;

for i=1:n
   if IMT(i)+1==result(i)
       right = right+1;
   end
end

disp(right);