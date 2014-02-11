
clear all
close all
load('testingFeatures.mat');
load('trainingFeatures.mat');

TrainingSet = FM; 
TestSet= FMT ; 
GroupTrain= IM; 

u=unique(GroupTrain);
numClasses=length(u);
result = zeros(length(TestSet(:,1)),1);
options = statset('maxiter',inf);

%build models
%for k=1:numClasses
%    G1vAll=(GroupTrain==u(k));
%    net(k) = svmtrain(TrainingSet,G1vAll,'Kernel_Function','rbf','options',options);
%    disp(k);
%end
load('net.mat');
TestSet(find(isnan(TestSet)))=0;
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
