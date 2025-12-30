% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);

X=rand(2,5);
n=size(X,2);
k=2;
M=diag([1,1,0.5,0.5]);
M=diag([1,1]);

%find the clusters with 
[clusterId,clusterCenter,clusterDia]=kMeans(X,M,k);
totalDia1=sum(clusterDia);
%Plot clusters
figure
gscatter(X(1,:),X(2,:),clusterId)
hold on;
for j=1:k
    plot(clusterCenter(1,j)+clusterDia(j)*circCos,clusterCenter(2,j)+clusterDia(j)*circSin)
end

%from MATLAB's kmeans function
clusterId=kmeans(X',k);

%Plot clusters
figure
gscatter(X(1,:),X(2,:),clusterId)


%%
%Manual check
Ba=de2bi(1:2^(n-1)-1,n);
%Ba=Ba';
for jj=1:length(Ba)
B=Ba(jj,:)';
oneVec=ones(5,1);
Xc(:,2)=X*B/(oneVec'*B);
Xc(:,1)=X*(oneVec-B)/(oneVec'*(oneVec-B));
Xc_arr{jj}=Xc;
dist=[];
for j=1:2
    ind=find(B==j-1);
    for i=1:length(ind)
        dist{j}(i)=norm(X(:,ind(i))-Xc(:,j));
    end
    clusterDia(j,jj)=max(dist{j});
end
totalDia0(jj)=sum(clusterDia(:,jj));
end
[totalDia3,indMinDia3]=min(totalDia0);

figure
gscatter(X(1,:),X(2,:),Ba(indMinDia3,:))
hold on
for j=1:2
    plot(Xc_arr{indMinDia3}(1,j)+clusterDia(j,indMinDia3)*circCos,Xc_arr{indMinDia3}(2,j)+clusterDia(j,indMinDia3)*circSin)
end
