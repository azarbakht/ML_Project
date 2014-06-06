% K-means clustering algorithm
% Authors: Amir Azarbakht and Mandana Hamidi
% azarbaam@eecs.oregonstate.edu
% 2014-05-30

function dataLabel = Kmeans_clustering_algorithm(data, K)

numberOfTrainingData = size(data,1);
numberOfFeatureDimensions = size(data,2);

randCenterIndices = randperm(numberOfTrainingData, K);
clusterCenters = data(randCenterIndices,:);

dataLabel = zeros(numberOfTrainingData, 1);
dataDistanceFromClusterCenter = Inf(numberOfTrainingData, K);

% Assign each data point to a cluster: label each data point with a cluster
% number

changeInLabelsFlag = 1;

while changeInLabelsFlag == 1,
    
    changeInLabelsFlag = 0;
    
    for i = 1:numberOfTrainingData,
        for j = 1:K,
            dataDistanceFromClusterCenter(i,j) = EuclideanDistance(data(i,:),clusterCenters(j,:));
        end
        % check to see if any data point's label is updated
        [~, newLabel] = min(dataDistanceFromClusterCenter(i,:));
        if  dataLabel(i) ~= newLabel,
            changeInLabelsFlag = 1;
        end
        % update label maybe?
        [~, dataLabel(i)] = min(dataDistanceFromClusterCenter(i,:));
        
    end
    
    
    % Re-estimate cluster centers
    
    Mu = Inf(K, numberOfFeatureDimensions);
    for i = 1:K,
        Mu(i,:) = mean(data(find(dataLabel == i),:));
        minDiff = EuclideanDistance(data(1,:), Mu(i,:));
        
        for j = 2:numberOfTrainingData
            if minDiff > EuclideanDistance(data(j,:), Mu(i,:)),
                minDiff = EuclideanDistance(data(j,:), Mu(i,:));
                closestDataPointToMuIndex = j;
            end
        end
        clusterCenters(i,:) = Mu(i,:);
    end
end

withinClusterSum = zeros(K,1);

% clustering done
% calculate withing cluster sum of distances
for i = 1:K,
    clusterI = data(find(dataLabel == i),:);
    for j = 1:size(clusterI,1),
        withinClusterSum(i) = withinClusterSum(i) + (EuclideanDistance(Mu(i,:), clusterI(j,:)))^2;
    end
end
% 
% title({'K-Means clustering with K = 5, for 200 times'});
% legend('Original Data Points', 'Cluster Centers');
% xlabel('X');
% ylabel('Y');

% save as PNG and EPS
% saveas(1, 'fig_Kmeans_with_random_initialization', 'epsc2');
% saveas(1, 'fig_Kmeans_with_random_initialization', 'png');
% saveas(1, 'fig_Kmeans_with_random_initialization', 'fig');
% saveas(1, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'Kmeans_with_random_initialization'], 'epsc2');
% saveas(1, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'Kmeans_with_random_initialization'], 'fig');
% saveas(1, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'Kmeans_with_random_initialization'], 'png');
% 
% % save the environment variables
save (['K-means_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS')]);
% 

end


