function [clusterId,Xc,clusterDia]=kMeans(X,M,k)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This function calculates  k cluster of a given data points X=[x_1, x_2, x_3,... x_n] with distance
%metric M such that distance between x_i and x_j is d(x_i,x_j)=sqrt((x_i-x_j)'*M*(x_i-x_j))

n=size(X,2);
indc0=randi([1,n],k,1);
Xc0=X(:,indc0);
err=10;
clusterId=[];
while err>1e-10
    for i=1:n
        for j=1:k
            dist(i,j)=sqrt((X(:,i)-Xc0(:,j))'*M*(X(:,i)-Xc0(:,j)));
        end 
        [~,clusterId(i)]=min(dist(i,:));
    end
    %Find new cluster centers
    err=0;
    Xc=[];
    for j=1:k
        Xj=X(:,clusterId==j);
        Xc(:,j)=sum(Xj,2)/size(Xj,2);
        err=err+sqrt((Xc0(:,j)-Xc(:,j))'*M*(Xc0(:,j)-Xc(:,j)));
        clusterDia(j)=max(dist(clusterId==j,j));
    end
    Xc0=Xc;
end

end