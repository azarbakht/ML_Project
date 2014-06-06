
function [FaultNumber, SelectedDataIndices] = activeLearning (data,kMeansLabelCluster1Data,  kmedoidsLabelCluster1Data, nystromCluster1Data)


SelectedDataIndices = [];
SelectedLabels = [];
%1- find the center of the clusters, All Clustering algorithms should have the same vote for the
%selected data point

%2- get the label of the selected data point 
SelectedLabels = [];
%3- 









FaultNumber(1)=1;


for i = 2 : size(SelectedDataIndices , 2)

    FaultNumber(i) = size(unique(SelectedLabels(1:i,1) ),1 );
    
end

end 