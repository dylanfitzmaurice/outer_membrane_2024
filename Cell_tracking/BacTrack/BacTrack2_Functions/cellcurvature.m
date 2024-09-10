%K=cellcurvature(X,Y,neighbor)
%This program calculates the curvature of a cell contour and averages over
%the neighboring points.
%
%INPUT
%X,Y: x and y coordinates of cell boundaries.
%neighbor: number of points to average curvature over.

function K=cellcurvature(X,Y,neighbor)

         [sX,temp]=size(X);

         dX=diff(X);
         dY=diff(Y);
         
         dS=sqrt(dX.^2+dY.^2);
         S=[0 cumsum(dS)'];
         DS=S(end)/(sX-1);
         
         NNx=zeros(sX-1,1);
         NNy=zeros(sX-1,1);
         
         NNx=dY./DS;
         NNy=-dX./DS;
         NN3=[NNx,NNy,zeros(sX-1,1)];
         
         Nlen=sqrt(NNx.^2+NNy.^2);
         
         Ndot=NNx(1:end-1).*NNx(2:end)+NNy(1:end-1).*NNy(2:end);
         Ncross=cross(NN3(1:end-1,:),NN3(2:end,:));
         dphi=acos(Ndot./(Nlen(1:end-1).*Nlen(2:end))).*sign(Ncross(:,3));
         
         K=zeros(sX,1);
         K(2:end-1)=dphi/DS;
         
         K=smooth(K,neighbor,'moving');