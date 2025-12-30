function PlotShortestPath2(tanG,Path,pathColor,tanG_prime,figNumber)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

SetPlotDefaults;
fontSize=18;
gammaMax=0.99999999;

path=Path.nodes;
GammaOk=tanG.GammaOk;

NO=tanG.NO;
nVO=tanG.nVO;
xOk=tanG.xOk;
yOk=tanG.yOk;
dxOk_dg=tanG.dxOk_dg;
dyOk_dg=tanG.dyOk_dg;
GammaOk=tanG.GammaOk;
PeriSOk=tanG.PeriSOk;
PeriO=tanG.PeriO;
G=tanG_prime.G;
G_pathType=tanG_prime.G_pathType;
rVO=tanG_prime.rVO;
rVO2=tanG_prime.rVO2;
gTO_all=tanG_prime.gTO_all;
rTO_all=tanG_prime.rTO_all;
obsId_all=tanG_prime.obsId_all;

%generate the path between the nodes on the shortest path
r_path=[];
Np=50;
countV0=1;
rV_path(:,1)=rTO_all(:,path(1));
for i=1:length(path)-1
    v1=path(i);
    v2=path(i+1);
    r_path0=[];
    if G_pathType(v1,v2)==3  %Tangent segments
        r_path0=[rTO_all(:,v1),rTO_all(:,v2)];
        k=obsId_all(v2);
    else  %Boundary segments
        g1=gTO_all(v1);
        g2=gTO_all(v2);
        k=obsId_all(v2);
        if G_pathType(v1,v2)==1            
            if g1<=g2
                gamma=g1:(g2-g1)/Np:g2;
                ind1=find(g1<GammaOk{k} & GammaOk{k}<g2);
            else
                gamma=g2:(1-g2)*2/Np:1;
                gamma=[gamma,0:2*g1/Np:g1];
                
                ind1=[find(g2<GammaOk{k} & GammaOk{k}<gammaMax) find(0<=GammaOk{k} & GammaOk{k}<g1)];
            end            
        elseif G_pathType(v1,v2)==2
            if g1>=g2
                gamma=g1:(g2-g1)/Np:g2;
                ind1=find(g2<GammaOk{k} & GammaOk{k}<g1);
            else
                gamma=g1:-2*g1/Np:0;
                gamma=[gamma,1:(g2-1)*2/Np:g2];
                ind1=[find(g2<GammaOk{k} & GammaOk{k}<gammaMax); find(0<=GammaOk{k} & GammaOk{k}<g1)];
            end
            ind1=flipud(ind1);
        end
        if ~isempty(ind1)
        rV_path(:,countV0+1:countV0+length(ind1))=rVO2{k}(:,ind1);
        rV_path_type(:,countV0+1:countV0+length(ind1))=1;  %intersection of straight line and circular ars, (outer vertices)
        rV_obsId(countV0+1:countV0+length(ind1))=k;
        end
        countV0=countV0+length(ind1);             
        
        g0=zeros(NO,1);
        for j=1:length(gamma)
            g0(k)=gamma(j);
            r_path0(:,j)=[xOk{k}(g0); yOk{k}(g0)];
        end
    end
     %Add the second vertex on the segment
        countV0=countV0+1;
        rV_path(:,countV0)=rTO_all(:,v2);
        rV_path_type(:,countV0)=0; %Tangent points
        rV_obsId(countV0)=k;
    r_path=[r_path,r_path0];
end
NrV_path=length(rV_path(1,:));
rV_obsId(NrV_path)=0;  %replace the obsId of the last vertex (final point) by 0
%Remove the intermediate vertices which are redundent  (only keep first, last and intermediate outer vertices on each obstacle in the path)
obsind=find(rV_obsId==rV_obsId(2));
rV_path2=rV_path(:,1);
while obsind(1)~=1 && obsind(1)~=NrV_path
    ind1=[obsind(1) obsind(find(rV_path_type(obsind)==1)) obsind(end)];
    rV_path2=[rV_path2,rV_path(:,ind1)];
    %find the indices for the next obstacle
    obsind=find(rV_obsId==rV_obsId(obsind(end)+1));
end
rV_path2=[rV_path2 rV_path(:,end)];
NrV_path2=length(rV_path2(1,:));
figure(figNumber)
hold on
if(0)
%all the vertices
plot(rV_path(1,:),rV_path(2,:),'mo')
%the outer vertices
ind2=find(rV_path_type==1);
plot(rV_path(1,ind2),rV_path(2,ind2),'go')
end

%minimum required vertices
plot(rV_path2(1,:),rV_path2(2,:),'o','color',pathColor)
%the path
plot(r_path(1,:),r_path(2,:),'color',pathColor,'linewidth',1.5)

%Add text corresponding to the speeds at intersection points
%initial and final
if (0)
text(rV_path2(1,1)-1,rV_path2(2,1),['$v_{i}$'],'fontsize',fontSize);
text(rV_path2(1,NrV_path2)-1,rV_path2(2,NrV_path2),['$v_{f}$'],'fontsize',fontSize)
%rest
for i=2:NrV_path2-1
    text(rV_path2(1,i)-1,rV_path2(2,i),['$v_{',num2str(i-1),'}$'],'fontsize',fontSize)
end
end
end
