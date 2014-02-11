function net = TrainMultiSVM( FM, IM )

TrainingSet = FM; 
GroupTrain= IM; 
u=unique(GroupTrain);
numClasses=length(u);
options = statset('maxiter',inf);

for k=1:numClasses
    G1vAll=(GroupTrain==u(k));
    net(k) = svmtrain(TrainingSet,G1vAll,'Kernel_Function','rbf','options',options);
end

end

