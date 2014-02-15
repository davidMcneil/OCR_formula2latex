function net = SVMTrain( FM, IM )
    % Train a multiclass SVM
    % Return a multidimensional net object with multiple nets
    addpath('svm');
    TrainingSet = FM; 
    GroupTrain= IM; 
    u=unique(GroupTrain);
    numClasses=length(u);

    for k=1:numClasses
        disp(['Starting Class: ' int2str(u(k))])
        G1vAll=(GroupTrain==u(k));
        net(k) = svmtrain(TrainingSet,G1vAll,'Kernel_Function','rbf');
    end

end

