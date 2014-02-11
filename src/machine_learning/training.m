
clear all
close all
load('feature.mat')

TrainingSet = FM; 
TestSet= FMTest ; 
GroupTrain= IM; 

u=unique(GroupTrain);
numClasses=length(u);
result = zeros(length(TestSet(:,1)),1);
options = statset('maxiter',inf);

%build models
for k=1:numClasses
    %Vectorized statement that binarizes Group
    %where 1 is the current class and 0 is all other classes
    G1vAll=(GroupTrain==u(k));
    net(k) = svmtrain(TrainingSet,G1vAll,'Kernel_Function','rbf','options',options);
    disp(k);
end

%classify test cases
for j=1:size(TestSet,1)
    for k=1:numClasses
        if(svmclassify(net(k),TestSet(j,:))) 
            break;
        end
    end
    result(j) = k;
end

disp('multi class problem'); 
disp(result);
plotboundary(models(1), [0,20], [0,20]);