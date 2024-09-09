%Eraseimagepart.m
%Rico Rojas, updated 1/21/19
%This routine lets the user select portions of an image stack to erase (make background).   

%INSTRUCTIONS FOR USE:
%Save the image stack in a directory without any other .tif files in it.  When 
%you run the program, the final image in the image stack will open.  
%Select the regions you want to delete with the cursor and then press Enter.  
%The program writes over the original image stack, so if you want a backup stack, 
%save it in a separate location.

%INPUT:
%dirname:Directory in which the image stack is saved.

%Calls upon:
%norm16bit.m
%%
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='FC1_090324_P3';
dirname={['/Users/dylanfitzmaurice/Documents/MATLAB/Matlab Ready/' basename '/' basename '_1_a']};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curdir=cd;

cd(dirname{1,1});
directory=dir('*.tif');
T=length(directory);
path(dirname{1,1},path)

imagename=directory(1).name; % change to T to show last frame
im1=imread(imagename);
[imM,imN]=size(im1);

ppix=0.0;
im2=norm16bit(im1,ppix);

figure,imshow(im1,stretchlim(im1)*65000)
set(gcf,'Pointer','crosshair')
hold on
axis manual

count=0;
k=0;
while k~=1
    count=count+1;
    k=waitforbuttonpress;
    point1=get(gca,'CurrentPoint');   
    finalRect=rbbox;                   
    %pause
    point2=get(gca,'CurrentPoint');    
    point1=point1(1,1:2);              
    point2=point2(1,1:2);
    point1(point1<1)=1;
    point2(point2<1)=1;
    if point1(2)>imM
        point1(2)=imM;
    end
    if point1(1)>imN
        point1(1)=imN;
    end
    if point2(2)>imM
        point2(2)=imM;
    end
    if point2(1)>imN
        point2(1)=imN;
    end
    p1=min(point1,point2);%Calculate locations
    p2=max(point1,point2);
    offset = abs(point1-point2);%And dimensions
    x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
    y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
    plot(x,y)
    
    rp1(count,:)=round(p1);
    rp2(count,:)=round(p2);
end

for i=1:1

cd(dirname{i,1});
directory=dir('*.tif');
T=length(directory);
path(dirname{i,1},path)

for t=1:T
    t
    for n=1:count
        %Load image
        imagename=directory(t).name;
        
        im=imread(imagename);
        [imM,imN]=size(im);
        
        [imcounts,bins]=hist(double(nonzeros(im1)));
        [~,mpos]=max(imcounts);
        val=bins(mpos);
        if i == 1
        im(rp1(n,2):rp2(n,2),rp1(n,1):rp2(n,1))=val*ones(rp2(n,2)-rp1(n,2)+1,rp2(n,1)-rp1(n,1)+1);
        else
        im(rp1(n,2):rp2(n,2),rp1(n,1):rp2(n,1))=val*zeros(rp2(n,2)-rp1(n,2)+1,rp2(n,1)-rp1(n,1)+1);%was zeros
        end
        delete(imagename);
        imwrite(im,imagename);
    end
end
end

cd('/Users/dylanfitzmaurice/Documents/MATLAB/data');
close all

beep