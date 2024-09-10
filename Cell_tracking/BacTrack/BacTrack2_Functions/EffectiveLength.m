function [l,varargout]=EffectiveLength(t,v,varargin);

%Integrates time course of strain rate to find length.

dt=diff(t);
meanV=(v(1:end-1)+v(2:end))/2;

l=zeros(size(v));
l0=3;
l(1)=l0;

for i=2:length(l)
    l(i)=l(i-1)+l(i-1)*meanV(i-1)*dt(i-1);
end

if length(varargin)==1
    varargout=cell(1,2);
    lmax=zeros(size(v));
    lmin=zeros(size(v));
    vstd=varargin{1};
    
    lmax(1)=l(1);
    lmin(1)=l(1);
    for i=2:length(l)
        lmax(i)=l(i-1)+l(i-1)*(meanV(i-1)+vstd(i-1))*dt(i-1);
        lmin(i)=l(i-1)+l(i-1)*(meanV(i-1)-vstd(i-1))*dt(i-1);
    end   
    varargout{1}=lmax;
    varargout{2}=lmin;
end
