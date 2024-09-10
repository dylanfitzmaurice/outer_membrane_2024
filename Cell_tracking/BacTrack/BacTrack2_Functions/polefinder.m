function [X,Y,pole]=polefinder(X,Y);

         K=cellcurvature(X,Y,1);
         lK=length(K);

         sigma=15;
         b=3;
         win=50;
         testx=[-win/2:win/2];
         test=b*exp(-testx.^2/2/(sigma/2).^2);
         cK=conv([K;K;K],test,'same');
         
%          dS=sqrt((X(2:end)-X(1:end-1)).^2+(Y(2:end)-Y(1:end-1)).^2);
%          S=[0 cumsum(dS)']*0.08;
%          figure,plot(S,K),hold on,
%          figure,plot(S,cK(lK+1:2*lK))
%          %figure, hold on, axis equal, axis off,plot(X,Y)
%          %pause
         
         [kmax,imax]=extrema(cK);
         tkmax=imax>lK & imax<=2*lK;
         imax=imax(tkmax);
         
         if length(imax)>1
            poles=imax(1:2)';
         else
            poles(1)=mod(imax(1),lK)+1;
            poles(2)=mod(poles(1)+round(lK/2),lK)+1;
         end
         mpole=min(poles);
         poles=poles-mpole+1;
         pole=max(poles);
         
         X=circshift(X(1:end-1),-mpole);
         Y=circshift(Y(1:end-1),-mpole);
         X=[X;X(1)];
         Y=[Y;Y(1)];

