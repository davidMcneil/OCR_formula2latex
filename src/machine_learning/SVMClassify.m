function [ result, rate ] = SVMClassify( net, FMT, IMT, numClasses )

result = zeros(size(IMT,1),1);
TestSet= FMT ; 
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

numData = size(result,1);
right = 0;

for i=1:numData
   if IMT(i)+1==result(i)
       right = right+1;
   end
end

rate = right/numData;

end

