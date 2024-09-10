function y=movingaverage(x,s);
%

[mx,lx]=size(x);
y=zeros(mx,lx);

if s>lx
    s=round(lx/2);
end

for i=1:lx
        s2=floor(s/2);
    if i-s2<=0
        y(:,i)=nanmean(x(:,1:i+s2)')';
    elseif i+s2>lx
        y(:,i)=nanmean(x(:,i-s2:end)')';
    else
        y(:,i)=nanmean(x(:,i-s2:i+s2)')';
    end
end