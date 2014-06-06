% main fucntion
clc
clear all
close all

% 0- read the dataset

%  datafcov = dsread('../fuzzer_data/d.fcov.ds');
%  dataTokenLabel = dsread('../fuzzer_data/fuzzer.fmin.data'); % this dataset contains Bag of words representation of TOKENS for each test case and its labels
%  datafcovLabel= dataTokenLabel(:,end);
%  save dataSet

load dataSet

K=100;

% Base line; random 
[FaultRand , SelectedDataIndicesRand] = baseline_Random (datafcov, datafcovLabel, K);

% Spectral clustering 
[nystromLabel, ~,~,~] = nystrom_no_orth(datafcov, 300, 50, K);


% Kmedoid clustering
kmedoidsLabel = kmedoids(datafcov', K);
% clustering K-means
kmeansLabel = Kmeans_clustering_algorithm(datafcov,K);
kmedoidsLabel = kmedoidsLabel';
% Hierarchical clustering algorithm
X = datafcov;
Y = pdist(X);
% Y2 = squareform(Y);

% hierarchical centroid
% Z = linkage(Y, 'centroid');

% hierarchical centroid
% Zcentroid = linkage(Y, 'centroid');

% hierarchical complete link
Z = linkage(Y, 'complete');

% % hierarchical single link
% Zsingle = linkage(Y, 'single');
% 
% % hierarchical average link
% Zaverage = linkage(Y, 'average');
% 
% % hierarchical weighted link
% Zweighted = linkage(Y, 'weighted');
% 
% % hierarchical weighted link
% Zward = linkage(Y, 'ward');

hierarchicalLabel = cluster(Z,'maxclust', K);
% hierarchicalLabelZcentroid = cluster(Zcentroid,'maxclust', K);
% hierarchicalLabelZcomplete = cluster(Zcomplete,'maxclust', K);
% hierarchicalLabelZsingle = cluster(Zsingle,'maxclust', K);
% hierarchicalLabelZaverage = cluster(Zaverage,'maxclust', K);
% hierarchicalLabelZweighted = cluster(Zweighted,'maxclust', K);
% hierarchicalLabelZward = cluster(Zward,'maxclust', K);


% [h,nodes] = dendrogram(Z,0);

% Transformation

% Labels = zeros(size(kmeansLabel,1),4);
Labels = [kmeansLabel kmedoidsLabel nystromLabel hierarchicalLabel datafcovLabel];

transformedKmeans = zeros(K,size(kmeansLabel,1));
transformedKmedoid = zeros(K,size(kmedoidsLabel,1));
transformedSpectral = zeros(K,size(nystromLabel,1));
transformedHierarchical = zeros(K,size(hierarchicalLabel,1));

for i = 1 : K,
    transformedKmeans(i,1:size(find(Labels(:,1)==i),1)) = (find(Labels(:,1)==i))';
    transformedKmedoid(i,1:size(find(Labels(:,2)==i),1)) = (find(Labels(:,2)==i))';
    transformedSpectral(i,1:size(find(Labels(:,3)==i),1)) = (find(Labels(:,3)==i))';
    transformedHierarchical(i,1:size(find(Labels(:,4)==i),1)) = (find(Labels(:,4)==i))';
end

transformedDatas(1).val = transformedKmeans;
transformedDatas(2).val = transformedKmedoid;
transformedDatas(3).val = transformedSpectral;
transformedDatas(4).val = transformedHierarchical;

colors={'r','g','b','c', 'k'};

hold off;
figure(10);
for counterClusteringMethod=1:4
    FaultNumber=[];
    SelectedDataIndices=[];
    for i=1:K,
        
        vector = (transformedDatas(counterClusteringMethod).val(i, find(transformedDatas(counterClusteringMethod).val(i,:)>0)));
        if (size(vector,2)>0)
            SelectedDataIndices(i) = vector(ceil(size(vector,2)/2));
        end
    end
    SelectedLabels = datafcovLabel(SelectedDataIndices(find(SelectedDataIndices>0)),:);
    FaultNumber (1)= 1;
    for j = 2 : size(SelectedLabels , 1)
        FaultNumber(j) = size(unique(SelectedLabels(1:j,1) ),1 );
    end
    
    hold on;
    plot((1:1:size(SelectedLabels , 1)), FaultNumber, colors{counterClusteringMethod});
    hold on;
    transformedDatas(counterClusteringMethod).FaultNumber=FaultNumber;
    transformedDatas(counterClusteringMethod).SelectedLabels=SelectedLabels;

    transformedDatas(counterClusteringMethod).SelectedDataIndices=SelectedDataIndices;

    
end

hold on;
plot((1:1:size(SelectedDataIndicesRand , 2)), FaultRand, colors{counterClusteringMethod+1});
legend ('Kmeans','Kmedoid','Spectral','Hierarchical', 'Baseline');

saveas(10, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'ClusteringAlgorithmsComparison'], 'epsc2');
saveas(10, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'ClusteringAlgorithmsComparison'], 'fig');
saveas(10, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'ClusteringAlgorithmsComparison'], 'png');

allBugLabels=unique([ unique(transformedDatas(1,1).SelectedLabels); ...
    unique(transformedDatas(1,2).SelectedLabels);...
     unique(transformedDatas(1,3).SelectedLabels);...
      unique(transformedDatas(1,4).SelectedLabels)]);

allBugLabels34=unique([ unique(transformedDatas(1,3).SelectedLabels);...
      unique(transformedDatas(1,4).SelectedLabels)]);

allBugLabels24=unique([ unique(transformedDatas(1,2).SelectedLabels);...
      unique(transformedDatas(1,4).SelectedLabels)]);

allBugLabels23=unique([ unique(transformedDatas(1,2).SelectedLabels);...
      unique(transformedDatas(1,3).SelectedLabels)]);

allBugLabels13=unique([ unique(transformedDatas(1,1).SelectedLabels);...
      unique(transformedDatas(1,3).SelectedLabels)]);

allBugLabels14=unique([ unique(transformedDatas(1,1).SelectedLabels);...
      unique(transformedDatas(1,4).SelectedLabels)]);

allBugLabels12=unique([ unique(transformedDatas(1,1).SelectedLabels);...
      unique(transformedDatas(1,2).SelectedLabels)]);

allBugLabels234=unique([unique(transformedDatas(1,2).SelectedLabels);...
     unique(transformedDatas(1,3).SelectedLabels);...
      unique(transformedDatas(1,4).SelectedLabels)]);

allBugLabels134=unique([unique(transformedDatas(1,1).SelectedLabels);...
     unique(transformedDatas(1,3).SelectedLabels);...
      unique(transformedDatas(1,4).SelectedLabels)]);

  
  
% save('K_51.mat');



% find similar clusters

% for i = 1 : K
%     maxIntersect = -1;
%     selectedJ = -1;
%     selectedk = -1;
%     
%     for j = 1 : K
% 
%         for k = 1 : K
%             
%             if maxIntersect <= size(intersect(intersect(transformedKmeans(i,:), transformedKmedoid(j,:)), transformedSpectral(k,:)),2),
%                 selectedJ = j;
%                 selectedk = k;
%                 maxIntersect = size(intersect(intersect(transformedKmeans(i,:), transformedKmedoid(j,:)), transformedSpectral(k,:)),2);
%             end % end if
%             
%             
%         end % end k
%         
%     end % end j
%     
%     transformedKmedoidNew(i,:) = transformedKmedoid(selectedJ,:);
%     transformedSpectralNew(i,:) = transformedSpectral(selectedk,:);
%     
% end %i

% select the best k number of clusters

% 2-2 code of the active learning


% save('K_28');
