%% Neal Ramseier 12.19.25 
%This code inputs the X and Y coordinates from ThunderSTORM reconstruction and clusters the data using DBSCAN, then 
%calculates the  nearest neighbor distance, maximum diameter, and area of the clusters. The values are then exported 
% in .csv file format. 
%% Close all and Import Data 
clear variables; clc; %close all; close all hidden;
disp("Previous Data Cleared"); 
if not(isfolder('MATLABData'))
    mkdir('MATLABData')
end
folder = string(pwd)+'\MATLABData\';
[fname,directory] = uigetfile('*.csv','Please choose a .csv file'); %the user selects a file
tic;
fileName = fullfile(directory, fname); %takes directory/filename.txt and stores to fullname var
Data = readmatrix(fileName);%reads in file to Data var
RawXData = Data(:,1); %separate X data
RawYData = Data(:,2); %separate Y data
RawXYData = [RawXData,RawYData]; 
Xmin = min(RawXData); Xmax = max(RawXData); Ymin = min(RawYData); Ymax = max(RawYData); 
figure(); %create figure of unclustered data
pr = scatter(RawXData, RawYData,10,'filled'); 
axis equal; xlim([Xmin-1000 Xmax+1000]); ylim([Ymin-1000 Ymax+1000]); 
disp("Data Loaded");

%% Cluser using DBSCAN
%This section of code clusters the data using DBSCAN. The clustered data is then denoised, then organized into a struct for grouping. 

%Change these two values below to tweak the clustering of the data.
minPTS = 17; %Set to the minimum nuber of points you want in a single cluster
Epsilon = 20; %Set the search radius to find those minimum number of points

Clusters = struct('XYG',{},'Area',{},'MaxDiameter',{},'NND',{},'NND_ID',{}); clusterLabel = [];
disp("Begin Clustering");
ClusterID = dbscan(RawXYData,Epsilon,minPTS); %perform the DBSCAN clustering of the X and Y data
StructData=[RawXData,RawYData,ClusterID]; %set up for removing noise points
[ClusterID2,val]=find(StructData==-1); %find noise points  
StructData(ClusterID2,:)=[]; % remove noise points
X = StructData(:,1); Y = StructData(:,2); ID = StructData(:,3);
for n = 1:max(ClusterID)
    temp = ID(StructData(:,3)==n);
    for iii =  1:length(temp)
        clusterLabel(iii,1) = n;
    end
    holdData = horzcat(X(StructData(:,3)==n), Y(StructData(:,3)==n),clusterLabel); %find clusters of 1:X and sort into XY 
    Clusters(n).XYG = holdData;
    clear clusterLabel
end
figure()
pl = gscatter(X,Y,ID,[],[],10);
TheLegend = legend({' '}); 
set(TheLegend,'visible','off'); 
xlabel("X"); ylabel("Y"); title("Clusters");
axis equal; xlim([Xmin-1000 Xmax+1000]); ylim([Ymin-1000 Ymax+1000]);
return
%% Attribute Calculation of Clusters
%This section of code calculates the max diameter of each cluster. 
disp("Calculating MD & Area");
%Area
nameCluster = 'Cluster'; holdArea = []; holdCluster = []; numberOfClusters = n; holdIDX = []; holdNND = []; MDHold = [];
for i =1:numberOfClusters 
    holdTheCluster = Clusters(i).XYG;
    X = holdTheCluster(:,1);
    Y = holdTheCluster(:,2);
    C = holdTheCluster(:,3);
    if length(X) < 4 %leave out clusters with less than 4 points
        continue
    end
    K = boundary(X,Y); %compute boundary
    area = polyarea(X(K),Y(K)); %Area of the boundary previously calculated.
    holdArea = [holdArea;area]; %add to matrix keeping track of area
    holdCluster = [holdCluster;i];%add to matrix keeping track of cluster number
    f =strcat(nameCluster,num2str(i),'Area'); 
    Clusters(i).Area = area; 

    %Maximum Diameter of Cluster
    xydistances = pdist2([X,Y],[X,Y]);
    maxDiameter = max(xydistances(:));
    MDHold = [MDHold;maxDiameter];
    f =strcat(nameCluster,num2str(i),'MaxDiameter'); 
    Clusters(i).MaxDiameter = maxDiameter; 
end
%% New NND Calculation
   %This section of code calculates a new NND. Rather than measuring the
   %distance of the center of a cluster to the nearest cluster center, this
   %code will calculte the distance of all points in the cluster to all
   %points of the nearest cluster. 
disp("Calculating NND");
for k = 1:numberOfClusters %loop through the clusters
    ClusterK = Clusters(k).XYG; %assign the cluster X and Y to varaibles
    Xk = ClusterK(:,1);
    Yk = ClusterK(:,2);
    Ck = [Xk,Yk]; %combine X and Y
    minDist = 3000000; %set minDist to something really high
    idx = -1; %initialize the index of the cluser number
    if length(X) < 4 %leave out clusters with less than 4 points
     continue
    end
    for j = 1:numberOfClusters %nested for loop to loop through clusters to compare with cluster from first for loop
       if j == k, continue; end %if the cluster being examined is the same as the one in this loop, skip it. 
       ClusterJ = Clusters(j).XYG; %assign the cluster X and Y to varaibles
       Xj = ClusterJ(:,1);
       Yj = ClusterJ(:,2);
       Cj = [Xj,Yj]; %combine X and Y
       dist = pdist2(Ck, Cj); %calculate the distances of all points between the two clusters
       d = min(dist(:)); %find the minimum value from the output of the pdist2 function
       if d < minDist %if statement to keep replacing lower numbers with the minimum 
           minDist = d;
           idx = j;
       end
    end
    holdNND = [holdNND;minDist]; %add the NND to the hold matrix
    Clusters(k).NND = minDist; %add to struct
    Clusters(k).NND_ID = idx; %add to struct
end

%display average Area and MD
averageMD = sum(MDHold)/length(MDHold);
disp(['average Max Diamater: ',num2str(averageMD), ' nm']);

averageArea = sum(holdArea)/length(holdArea);
disp(['average Area: ',num2str(averageArea), ' nm^2']);
averageNND = sum(holdNND)/length(holdNND); %calculate the average
disp(['Average NND: ', num2str(averageNND), ' nm']);
toc; %Time
