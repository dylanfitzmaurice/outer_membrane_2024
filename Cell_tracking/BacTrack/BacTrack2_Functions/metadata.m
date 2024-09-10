function time=metadata(filename);

%This function extracts the time in seconds at which images were taken from
%micromanager metadata.
%
%Input:
%filename - a string containing the full pathname of the metadata file.
%
%Output:
%time - a n x T row vector containing the time points, in seconds, at which
%       frames were captured, where n is the number of channels.

fid=fopen(filename);
count=0;
while ~feof(fid);
    line_str=fgetl(fid);
    if strfind(line_str,'ElapsedTime-ms')==4
        count=count+1;
        tpoints(count)=sscanf(line_str,' "ElapsedTime-ms":%f');
    end
end
fclose(fid);

chan=[];
fid=fopen(filename);
while isempty(chan)==1;
    line_str=fgetl(fid);
    if strfind(line_str,'Channels')==4;
        chan=sscanf(line_str,' "Channels":%f');
        break
    end
end
fclose(fid);

tpoints=sort(tpoints);
tpoints=tpoints-tpoints(1);
tpoints=tpoints/1000;

T=floor(length(tpoints)/chan);
time=zeros(chan,T);

for n=1:chan
    eval(['time(n,:)=tpoints(n:chan:T*chan);']);
end