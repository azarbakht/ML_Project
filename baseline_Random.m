function [FaultNumber, SelectedDataIndices] = baseline_Random (data, datalabel, numberOfRuns)

SelectedDataIndices = randperm(size(data,1), numberOfRuns);
SelectedLabels = datalabel(SelectedDataIndices,:);
FaultNumber (1)= 1;
for i = 2 : size(SelectedDataIndices , 2)
    FaultNumber(i) = size(unique(SelectedLabels(1:i,1) ),1 );
end

end