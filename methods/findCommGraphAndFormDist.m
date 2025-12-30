function [W, Rij_tilde]=findCommGraphAndFormDist(N,shapeId,R0)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This function finds communication graph for given number of agents, shape
%and reference distances between agents based on the given parameter R0

%shapeId=0: straight line (R0 is inter-agent distance)
%shapeId=1: regular polygon formation with no center agent (R0 is radius of
%the circumcircle)
%shapeId=2: regular polygon formation with one center agent (R0 is radius of
%the circumcircle)
%shapeId=3: regular open semicircular formation with no center agent(R0 is radius of the
%circle)
%shapeId=4: regular open semicircular formation with one center agent(R0 is radius of the
%circle)
%shapeId=5: two tier regular polygon with no center agent (R0 is radius of
%the circumcircle, seconnd polygon with 2*R0 radius)
%shapeId=6:  Kite like shape (arrow like)

if shapeId==0   %straight line formation
    Rij0=R0;
    Rij1=(N-1)*R0;
    %Rij_bar=Rij1*ones(N);
    if N>1
        W=ones(N);
        Rij_tilde=zeros(N);
        for i=1:N
            W(i,i)=0;
            for ii=1:N
                Rij_tilde(i,ii)=abs(i-ii)*R0;
            end             
        end
    else
        W=0;
    end
    
elseif shapeId==1
    Rij_tilde=zeros(N);
    if N>1
        W=zeros(N);
        %NC=N-1;  %One center agent, others on the boundary
        NC=N;  % all at the periphery, No center agent
        Rij0=R0*sqrt(2*(1-cos(2*pi/NC)));
        Rij1=Rij0*sqrt(2*(1-cos(pi-2*pi/NC)));
        temp=[1:NC,1:NC];
        for i=1:NC
            W(i,temp([i+1:i+2,i+NC-1]))=ones(1,3);
            Rij_tilde(i,temp([i+1,i+NC-1]))=Rij0;
            Rij_tilde(i,temp(i+2))=Rij1;
            Rij_tilde(i,i)=0;           
        end
    else
        W=0;
    end
elseif shapeId==2
     
    Rij_tilde=zeros(N);
    if N>1
        W=zeros(N);
        NC=N-1;  %One center agent, others on the boundary
        Rij0=R0*sqrt(2*(1-cos(2*pi/NC)));
        Rij1=Rij0*sqrt(2*(1-cos(pi-2*pi/NC)));
        %NC=N;  % all at the periphery, No center agent
        temp=[1:NC,1:NC];
        for i=1:NC
            W(i,temp([i+1:i+2,i+NC-1]))=ones(1,3);
            Rij_tilde(i,temp([i+1,i+NC-1]))=Rij0;
            Rij_tilde(i,temp(i+2))=Rij1;
            Rij_tilde(i,i)=0;           
        end
        W(N,1:N-1)=ones(1,N-1);
        W(1:N-1,N)=ones(N-1,1);
        Rij_tilde(N,1:N-1)=R0;
        Rij_tilde(1:N-1,N)=R0;
    else
        W=0;
    end
    
elseif shapeId==3
    Rij_tilde=zeros(N);
    if N>1
        W=zeros(N);
        NC=N;  %All on the boundary
        Rij0=R0*sqrt(2*(1-cos(pi/(NC-1)))); %vertices on semicircular arc
        Rij1=Rij0*sqrt(2*(1-cos(pi-pi/(NC-1))));
        %NC=N;  % all at the periphery, No center agent
        temp=[1:NC,1:NC];
        for i=1:NC
            W(i,temp([i+1:i+2,i+NC-1]))=ones(1,3);
            Rij_tilde(i,temp([i+1,i+NC-1]))=Rij0;
            Rij_tilde(i,temp(i+2))=Rij1;
            Rij_tilde(i,i)=0;           
        end        
        %No connection between the first and the last defender 
        W(1,N)=0;
        W(N,1)=0;
        Rij_tilde(1,N)=0;
        Rij_tilde(N,1)=0;
    else
        W=0;
    end
    
elseif shapeId==4
     Rij_tilde=zeros(N);
    if N>1
        W=zeros(N);
        NC=N-1;  %One center agent, others on the boundary
        Rij0=R0*sqrt(2*(1-cos(pi/(NC-1)))); %vertices on semicircular arc
        Rij1=Rij0*sqrt(2*(1-cos(pi-pi/(NC-1))));
        %NC=N;  % all at the periphery, No center agent
        temp=[1:NC,1:NC];
        for i=1:NC
            W(i,temp([i+1:i+2,i+NC-1]))=ones(1,3);
            Rij_tilde(i,temp([i+1,i+NC-1]))=Rij0;
            Rij_tilde(i,temp(i+2))=Rij1;
            Rij_tilde(i,i)=0;           
        end
        W(N,1:N-1)=ones(1,N-1);
        W(1:N-1,N)=ones(N-1,1);
        Rij_tilde(N,1:N-1)=R0;
        Rij_tilde(1:N-1,N)=R0;
        
        %No connection between the first and the last defender (N^th agent is virtaul)
        W(1,N-1)=0;
        W(N-1,1)=0;
        Rij_tilde(1,N-1)=0;
        Rij_tilde(N-1,1)=0;
    else
        W=0;
    end
    
elseif shapeId==5
    W=ones(N);
    for i=1:N
      W(i,i)=0;
      for ii=1:N
      Rij_tilde(i,ii)=R0*abs(i-ii);
      end
    end
    
elseif shapeId==6
    W=ones(N);
    for i=1:N
      W(i,i)=0;
      if i==1
          Rij_tilde(1,2)=R0;
          Rij_tilde(1,3)=R0;
          for ii=4:N
             Rij_tilde(1,ii) = sqrt(2)*R0+(ii-4)*R0;
          end
      elseif i==2
          Rij_tilde(2,1)=R0;
          Rij_tilde(2,3)=sqrt(2)*R0;
          for ii=4:N
             Rij_tilde(2,ii) =sqrt(R0^2+(R0*(ii-4))^2-2*R0*(R0*(ii-4))*cos(3*pi/4));
          end
      elseif i==3
          Rij_tilde(3,1)=R0;
          Rij_tilde(3,2)=sqrt(2)*R0;
          for ii=4:N
             Rij_tilde(3,ii) =sqrt(R0^2+(R0*(ii-4))^2-2*R0*(R0*(ii-4))*cos(3*pi/4));
          end
      else
          Rij_tilde(i,1) = Rij_tilde(1,i);
          Rij_tilde(i,2) = Rij_tilde(2,i);
          Rij_tilde(i,3) = Rij_tilde(3,i);
           for ii=5:N
             Rij_tilde(i,ii) = (ii-i)*R0;
          end
      end
    end
else
    disp('No valid shape specified.\n')
end

end