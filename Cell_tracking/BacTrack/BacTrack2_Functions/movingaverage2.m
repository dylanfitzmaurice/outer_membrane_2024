function y=movingaverage2(x,s);

lx=length(x);
y=zeros(size(x));

if s>lx
    s=round(lx/2);
end

for i=1:lx
        s2=floor(s/2);
    if i-s2<=0
        y(i)=mean(x(1:i+s2));
    elseif i+s2>lx
        y(i)=mean(x(i-s2:end));
    else
        y(i)=mean(x(i-s2:i+s2));
    end
end