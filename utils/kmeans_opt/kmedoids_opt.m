function [IDX,C,SUMD,D,K,d,Var,PC]=kmedoids_opt(X,varargin)
%%% [IDX,C,SUMD,K]=kmeans_opt(X,varargin) returns the output of the k-medoids
%%% algorithm with the optimal number of clusters, as determined by the ELBOW
%%% method. this function treats NaNs as missing data, and ignores any rows of X that
%%% contain NaNs.
%%%
%%% [IDX]=kmeans_opt(X) returns the cluster membership for each datapoint in
%%% vector X.
%%%
%%% [IDX]=kmeans_opt(X,MAX) returns the cluster membership for each datapoint in
%%% vector X. The Elbow method will be tried from 1 to MAX number of
%%% clusters (default: square root of the number of samples)
%%% [IDX]=kmeans_opt(X,MAX,CUTOFF) returns the cluster membership for each datapoint in
%%% vector X. The Elbow method will be tried from 1 to MAX number of
%%% clusters and will choose the number which explains a fraction CUTOFF of
%%% the variance (default: 0.95)
%%% [IDX]=kmeans_opt(X,MAX,CUTOFF,REPEATS) returns the cluster membership for each datapoint in
%%% vector X. The Elbow method will be tried from 1 to MAX number of
%%% clusters and will choose the number which explains a fraction CUTOFF of
%%% the variance, taking the best of REPEATS runs of k-means (default: 3).
%%% [IDX,C]=kmeans_opt(X,varargin) returns in addition, the location of the
%%% centroids of each cluster.
%%% [IDX,C,SUMD]=kmeans_opt(X,varargin) returns in addition, the sum of
%%% point-to-cluster-centroid distances.
%%% [IDX,C,SUMD,K]=kmeans_opt(X,varargin) returns in addition, the number of
%%% clusters.

%%% sebastien.delandtsheer@uni.lu
%%% sebdelandtsheer@gmail.com
%%% Thomas.sauter@uni.lu


[m,~]=size(X); %getting the number of samples

if nargin>1, ToTest=cell2mat(varargin(1)); else, ToTest=ceil(sqrt(m)); end
if nargin>2, Cutoff=cell2mat(varargin(2)); else, Cutoff=0.95; end
if nargin>3, Repeats=cell2mat(varargin(3)); else, Repeats=3; end
if nargin>4, Distance=varargin{4}; else, Distance='sqeuclidean'; end

% %unit-normalize
% MIN=min(X'); MAX=max(X'); 
% X=(X'-MIN)./(MAX-MIN);X=X';
% %check normalization

d=zeros(ToTest,1); %initialize the results matrix
for c=1:ToTest %for each sample
    [~,~,dist]=kmedoids(X,c,'Distance',Distance); %compute the sum of intra-cluster distances
    tmp=sum(dist); %best so far
    
    for cc=2:Repeats %repeat the algo
        [~,~,dist]=kmedoids(X,c,'Distance',Distance);
        tmp=min(sum(dist),tmp);
    end
    d(c,1)=tmp; %collect the best so far in the results vecor
end

Var=d(1:end-1)-d(2:end); %calculate %variance explained
PC=cumsum(Var)/(d(1)-d(end));

[r,~]=find(PC>Cutoff); %find the best index
K=1+r(1,1); %get the optimal number of clusters
[IDX,C,SUMD,D]=kmedoids(X,K,'Distance',Distance); %now rerun one last time with the optimal number of clusters

%C=C.*(MAX'-MIN')+MIN';

end