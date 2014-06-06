% K-means clustering algorithm
% Authors: Amir Azarbakht and Mandana Hamidi
% azarbaam@eecs.oregonstate.edu
% 2014-05-30
% 
% clc;
% clear all;
% close all;
 startTime = tic;
% 

load dataSet
minK = 63;
maxK = 100;
numberOfRuns = 10;

minKwithinClusterSum = Inf(maxK,1);
run_K_WithinClusterSumDistance = [];


for K = minK : maxK,
    numberOfTrainingData = size(datafcov,1);
    numberOfFeatureDimensions = size(datafcov,2);
    
    for index = 1 : numberOfRuns,
        randCenterIndices = randperm(numberOfTrainingData, K);
        clusterCenters = datafcov(randCenterIndices,:);
        
        datafcovLabel = zeros(numberOfTrainingData, 1);
        data1DistanceFromClusterCenter = Inf(numberOfTrainingData, K);
        
        % Assign each data point to a cluster: label each data point with a cluster
        % number
        
        changeInLabelsFlag = 1;
        
        while changeInLabelsFlag == 1,
            
            changeInLabelsFlag = 0;
            
            for i = 1:numberOfTrainingData,
                for j = 1:K,
                    data1DistanceFromClusterCenter(i,j) = EuclideanDistance(datafcov(i,:),clusterCenters(j,:));
                end
                % check to see if any data point's label is updated
                [~, newLabel] = min(data1DistanceFromClusterCenter(i,:));
                if  datafcovLabel(i) ~= newLabel,
                    changeInLabelsFlag = 1;
                end
                % update label maybe?
                [~, datafcovLabel(i)] = min(data1DistanceFromClusterCenter(i,:));
                
            end
            
            
            % Re-estimate cluster centers
            Mu = Inf(K, numberOfFeatureDimensions);
            
            for i = 1:K,
                Mu(i,:) = mean(datafcov(find(datafcovLabel == i),:));
                minDiff = EuclideanDistance(datafcov(1,:), Mu(i,:));
                for j = 2:numberOfTrainingData
                    if minDiff > EuclideanDistance(datafcov(j,:), Mu(i,:)),
                        minDiff = EuclideanDistance(datafcov(j,:), Mu(i,:));
                        closestDataPointToMuIndex = j;
                    end
                end
                
                clusterCenters(i,:) = Mu(i,:);
            end
        end
        withinClusterSum = zeros(K, 1);
        
        % clustering done
        % calculate withing cluster sum of distances
        for i = 1:K,
            
            clusterI = datafcov(find(datafcovLabel == i),:);
            for j = 1:size(clusterI,1),
                withinClusterSum(i) = withinClusterSum(i) + (EuclideanDistance(Mu(i,:), clusterI(j,:)))^2;
            end
            
        end
        
        if min(withinClusterSum) == 0
            test = 0;
        end
        
        MinMaxMeanSD(index, 1) = sum(withinClusterSum);        
        
    end % for # of runs index
    
    kMinMaxMeanSD(K, 1) = min(MinMaxMeanSD);
    save (['K=' int2str(K) '_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS')]);

end % for K # of clusters index


title({'K-Means clustering with K = 2:15, for 10 times'});
legend('Original Data Points', 'Cluster Centers');
xlabel('X');
ylabel('Y');


% % plot the metrics
figure(2);
plot((minK:1:maxK), kMinMaxMeanSD(minK:1:maxK), 'b');
title({'Min Within-Cluster Sum of Squared Distances for K-Means clustering with K = 5, for 200 times'});
legend('Minimum Within-Cluster Sum of Squared Distances');
xlabel('# of runs');
ylabel('Within-Cluster Sum of Squared Distances');
 
saveas(2, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'Within_cluster_sum_of_squared_distances__K_2_15'], 'epsc2');
saveas(2, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'Within_cluster_sum_of_squared_distances__K_2_15'], 'fig');
saveas(2, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'Within_cluster_sum_of_squared_distances__K_2_15'], 'png');


% save the environment variables
save ([datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS')]);

runningTime = toc(startTime);
runningTime
