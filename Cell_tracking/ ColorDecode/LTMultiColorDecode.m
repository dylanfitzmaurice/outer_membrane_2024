%LTMultiColorDecode.m
%3/28/19, Rico Rojas, lines 1-87
%Dylan Fitzmaurice, 08/21/20, lines 88-end
%To be used after BacTrack of 3 different channels.  
%Decodes multiple bacterial strains based on fluorescent labels.
%Save color decode images as an image sequence in a separate directory with
%saved image order: phase, GFP, RFP, CY5, CY7, DAPI.
%For decouple data with CY5, DAPI, GFP images. 

clear, close all

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basename='FC1_090324_P1'; 
legend_str={'Green','Orange','Red','Infrared','Blue','GO','OB','IRB','GIR','GB','OIR''GIRB'}; 
thresholds=[7800 1500 19000 130 19000];%Set thresholds to determine cells with greater intenisty than threshold for each flourescent channel 
ncolors=6; %Number of channels used (phase,GFP,RFP,CY5,CY7,DAPI)

%%
%Directory in which movie is saved
movie_directory=['/Users/dylanfitzmaurice/Documents/MATLAB/Matlab Ready/' basename '/' basename '_1_a' ];%CY5
color_directory=['/Users/dylanfitzmaurice/Documents/MATLAB/Matlab Ready/' basename '/' basename '_col' ];%Directory in which color 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

curdir=cd;
cd(color_directory); 
fluo_directory=dir('*.tif');
cd(movie_directory); 
phase_directory=dir('*.tif');

    load([basename '_CY5' '_BT'],'B','boun','lcents','cellnum','tstamp','l','pixels','T','v','wv','av','tmid', 'ed10','labels','lcell','ncells','cellid','pxls')
    

%Set boundaries (B)
for i=1:size(B,1)
    for j=1:size(B,2)
        if isempty(B{i,j})==0
            nB{i,1}=B{i,j};
            break
        end
    end
end
    
    cd(color_directory); 
    imagename=fluo_directory(1).name;
    TestIm=imread(imagename);
    
    cd(movie_directory); 
    imagename=phase_directory(1).name; % 
    RefIm=imread(imagename);
    
    [output, NewImFT]=dftregistration(fft2(RefIm),fft2(TestIm),10);
    NewIm=abs(ifft2(NewImFT));
    shft(1)=0;
    shft(2)=0;
    
    vav_col=zeros(ncolors,T-1);
    
    cd(color_directory)
    
for i=2:ncolors
    
    imagename=fluo_directory(i).name;
    im=imread(imagename);
    
    im=imtranslate(im,shft);
    intensities_temp=zeros(ncells,1); %Creates a colum of zeros whose rows equal ncells
        for j=1:ncells
            intensities_temp(j)=mean(im(pixels{j,1}));%pixels{j,T})%1 for first frame, T for last
        end
    icell{i}=intensities_temp;
    
    positive_cells{i}=find(intensities_temp>thresholds(i-1)); % Determines which cells are fluorescent    
    v_col{i}=v(positive_cells{i},:);    
    vav_col(i,:)=nanmean(v_col{i});  

    figure,hist(nonzeros(intensities_temp)),pause
    % Shows histogram so you can manual set threshold to select most intense cells
    
    figure,
    imshow(im,[])
    hold on
    for k=1:length(positive_cells{1,i})
       plot(nB{positive_cells{1,i}(k),1}(:,1),nB{positive_cells{1,i}(k),1}(:,2),'-r') %first frame     
    end
    pause
end

%% Grouping dual-labeled cells 

% Grouping the cells based on color combo
GO = intersect(positive_cells{2},positive_cells{3});
OB = intersect(positive_cells{3},positive_cells{6});
IRB = intersect(positive_cells{5},positive_cells{6});
GIR = intersect(positive_cells{2},positive_cells{5});
GB = intersect(positive_cells{2},positive_cells{6});
OIR = intersect(positive_cells{3},positive_cells{5});

% % Finding GO cells
positive_cells{7} = GO;
%deleting green cells from original vector, for GO 
for i=1:size(positive_cells{2},1)% green pos cell size
    for j=1:size(positive_cells{7},1)% GO pos cell size before selection
        if positive_cells{2}(i,1)==positive_cells{7}(j,1) %if green and GO have duplicate cells 
            positive_cells{2}(i,1)=0;%turn GO to nans in green pos cell
        end
    end
end

% % Finding OB cells 
positive_cells{8} = OB;
%deleting orange cells from original vector, for OB 
for i=1:size(positive_cells{3},1)% orange pos cell size
    for j=1:size(positive_cells{8},1)% OB cell size before selection
        if positive_cells{3}(i,1)==positive_cells{8}(j,1)
            positive_cells{3}(i,1)=0;%turn OB to nans in orange pos cell
        end
    end
end

% % Finding IRB cells 
positive_cells{9} = IRB;
%deleting IR cells from original vector, for IRB 
for i=1:size(positive_cells{5},1)% IR pos cell size
    for j=1:size(positive_cells{9},1)% IRB cell size before selection
        if positive_cells{5}(i,1)==positive_cells{9}(j,1) %if infrared and IRB have duplicate cells 
            positive_cells{5}(i,1)=0;%turn IRB to nans in infrared pos cell
        end
    end
end

% % Finding GIR cells 
positive_cells{10} = GIR;
%deleting green cells from original vector, for GIR 
for i=1:size(positive_cells{2},1)% green pos cell size
    for j=1:size(positive_cells{10},1)% GIR cell size before selection
        if positive_cells{2}(i,1)==positive_cells{10}(j,1) %if green and GIR have duplicate cells 
            positive_cells{2}(i,1)=0;%turn GIR to nans in green pos cell
        end
    end
end

% % Finding GB cells 
positive_cells{11} = GB;
%deleting green cells from original vector, for GB 
for i=1:size(positive_cells{2},1)% green pos cell size
    for j=1:size(positive_cells{11},1)% GB cell size before selection
        if positive_cells{2}(i,1)==positive_cells{11}(j,1) %if green and GIR have duplicate cells 
            positive_cells{2}(i,1)=0;%turn GIR to nans in green pos cell
        end
    end
end

% % Finding OIR cells 
positive_cells{12} = OIR;
%deleting orange cells from original vector, for OIR 
for i=1:size(positive_cells{3},1)% orange pos cell size
    for j=1:size(positive_cells{12},1)% OIR cell size before selection
        if positive_cells{3}(i,1)==positive_cells{12}(j,1) %if orange and OIR have duplicate cells 
            positive_cells{3}(i,1)=0;%turn OIR to nans in orange pos cell
        end
    end
end


GIRB = intersect(positive_cells{10},positive_cells{6});

% % Finding GIRB cells 
positive_cells{13} = GIRB;
%deleting GIR cells from original vector, for GIRB 
for i=1:size(positive_cells{10},1)% green pos cell size
    for j=1:size(positive_cells{13},1)% GIR cell size before selection
        if positive_cells{10}(i,1)==positive_cells{13}(j,1) %if green and GIR have duplicate cells 
            positive_cells{10}(i,1)=0;%turn GIR to nans in green pos cell
        end
    end
end

% % deleting orange cells from original vector, for GO 
for i=1:size(positive_cells{3},1)% orange pos cell size
    for j=1:size(positive_cells{7},1)% GO pos cell size before selection
        if positive_cells{3}(i,1)==positive_cells{7}(j,1)%if orange and GO have duplicate cells 
            positive_cells{3}(i,1)=nan;%turn GO to nans in orange pos cell
        end
    end
end

%remove nans
positive_cells{2}=positive_cells{2}(~isnan(positive_cells{2}));%Remove nans from green pos cell
positive_cells{3}=positive_cells{3}(~isnan(positive_cells{3}));%Remove nans from orange pos cell

% % deleting blue cells from original vector, for OB 
for i=1:size(positive_cells{6},1)% blue pos cell size
    for j=1:size(positive_cells{8},1)% OB cell size before selection
        if positive_cells{6}(i,1)==positive_cells{8}(j,1)
            positive_cells{6}(i,1)=nan;%turn OB to nans in blue pos cell
        end
    end
end

%remove nans
positive_cells{3}=positive_cells{3}(~isnan(positive_cells{3}));%Remove nans from orange pos cell 
positive_cells{6}=positive_cells{6}(~isnan(positive_cells{6}));%Remove nans from blue pos cell 

% % deleting blue cells from original vector, for IRB 
for i=1:size(positive_cells{6},1)% blue pos cell size
    for j=1:size(positive_cells{9},1)% IRB cell size before selection
        if positive_cells{6}(i,1)==positive_cells{9}(j,1) %if blue and IRB have duplicate cells 
            positive_cells{6}(i,1)=nan;%turn IRB to nans in blue pos cell
        end
    end
end

%remove nans
positive_cells{5}=positive_cells{5}(~isnan(positive_cells{5}));%Remove nans from infrared pos cell 
positive_cells{6}=positive_cells{6}(~isnan(positive_cells{6}));%Remove nans from blue pos cell 

% % deleting infrared cells from original vector, for GIR 
for i=1:size(positive_cells{5},1)% infrared pos cell size
    for j=1:size(positive_cells{10},1)% GIR cell size before selection
        if positive_cells{5}(i,1)==positive_cells{10}(j,1) %if infrared and GIR have duplicate cells 
            positive_cells{5}(i,1)=nan;%turn GIR to nans in infrared pos cell
        end
    end
end

%remove nans
positive_cells{2}=positive_cells{2}(~isnan(positive_cells{2}));%Remove nans from green pos cell 
positive_cells{5}=positive_cells{5}(~isnan(positive_cells{5}));%Remove nans from infrared pos cell


% % deleting blue cells from original vector, for GIRB 
for i=1:size(positive_cells{6},1)% GIR pos cell size
    for j=1:size(positive_cells{13},1)% GIRB cell size before selection
        if positive_cells{6}(i,1)==positive_cells{13}(j,1) %if GIR and blue have duplicate cells 
            positive_cells{6}(i,1)=nan;%turn blues to nans in GIR pos cell
        end
    end
end

%remove nans
positive_cells{10}=positive_cells{10}(~isnan(positive_cells{10}));%Remove nans from GIR pos cell 
positive_cells{6}=positive_cells{6}(~isnan(positive_cells{6}));%Remove nans from blue pos cell


% % deleting infrared cells from original vector, for GB 
for i=1:size(positive_cells{6},1)% blue pos cell size
    for j=1:size(positive_cells{11},1)% IRB cell size before selection
        if positive_cells{6}(i,1)==positive_cells{11}(j,1) %if blue and GB have duplicate cells 
            positive_cells{6}(i,1)=nan;%turn GB to nans in infrared pos cell
        end
    end
end

%remove nans
positive_cells{2}=positive_cells{2}(~isnan(positive_cells{2}));%Remove nans from green pos cell 
positive_cells{6}=positive_cells{6}(~isnan(positive_cells{6}));%Remove nans from blue pos cell

% % deleting infrared cells from original vector, for OIR 
for i=1:size(positive_cells{5},1)% infrared pos cell size
    for j=1:size(positive_cells{12},1)% OIR cell size before selection
        if positive_cells{5}(i,1)==positive_cells{12}(j,1) %if infrared and OIR have duplicate cells 
            positive_cells{5}(i,1)=nan;%turn OIR to nans in infrared pos cell
        end
    end
end

%remove nans
positive_cells{3}=positive_cells{3}(~isnan(positive_cells{3}));%Remove nans from orange pos cell 
positive_cells{5}=positive_cells{5}(~isnan(positive_cells{5}));%Remove nans from infrared pos cell

%% Remove values of zero

for i=1:13
positive_cells{i}(positive_cells{i}==0)=NaN; % turn zeros to nans
positive_cells{i}=positive_cells{i}(~isnan(positive_cells{i})); % removes final nans
end

%% Delete any duplicates within single colored cells (if any) 

vdposcell1=vertcat(positive_cells{1,2},positive_cells{1,3},positive_cells{1,4},...
    positive_cells{1,5},positive_cells{1,6});

uni_1=unique(vdposcell1)==sort(vdposcell1); %if all 1's no duplicates
if size(uni_1,1)==sum(uni_1) 
    'no duplicates for single colors'
else
    'DUPLICATES!!! for single colors!'
end

%% Delete any duplicates within double colored cells (if any) 

vdposcell2=sort(vertcat(positive_cells{1,7},positive_cells{1,8},positive_cells{1,9},...
    positive_cells{1,10},positive_cells{1,11},positive_cells{1,12},positive_cells{1,13}));


%Find duplicate cells
U = unique(vdposcell2);
nU=U(1<histc(vdposcell2,unique(vdposcell2)));

for k = 7:13 %for double colored cells
    for i=1:size(positive_cells{1,k},1)
        for j = 1:size(nU,1)
            if positive_cells{1,k}(i,1)==nU(j,1)
                positive_cells{1,k}(i,1)=nan;
            end
        end
    end
end

for i=7:13
positive_cells{i}=positive_cells{i}(~isnan(positive_cells{i})); % removes final nans
end

vdposcell3=sort(vertcat(positive_cells{1,7},positive_cells{1,8},positive_cells{1,9},...
    positive_cells{1,10},positive_cells{1,11},positive_cells{1,12},positive_cells{1,13}));

uni_3=unique(vdposcell3)==sort(vdposcell3); %if all 1's no duplicates
if size(uni_3,1)==sum(uni_3) 
    'no duplicates for double colors'
else
    'DUPLICATES!!! for double colors!'
end

%% plot boundaries of stained cells 
% Example: plot cells from dyed with MitoView 720 infrared channel
figure,
imshow(im,[])
col=5;
hold on 
for k=1:length(positive_cells{1,col})
plot(B{positive_cells{1,col}(k),1}(:,1),B{positive_cells{1,col}(k),1}(:,2),'r')
end
hold off
%%

ncolors=13;%total number of combinations
for i = 2:ncolors
    v_col{i}=v(positive_cells{i},:);    
    vav_col(i,:)=nanmean(v_col{i});    
end

v_col{1}=v(positive_cells{1},:);
vav_col(1,:)=nanmean(v_col{1});
Leff_col=zeros(size(vav_col));

for i=1:ncolors
    Leff_col(i,:)=EffectiveLength(tmid,vav_col(i,:));
end

figure
hold on
plot(tmid,Leff_col(2:end,:)) %,linecolor(i))
hold off
xlabel('Time (s)')
ylabel('Effective Population Averaged Length (\mum)')
legend(legend_str)
fig2pretty

cd(curdir)

%% print number of cells decoded
for i=1:12
n(i)=size(positive_cells{1,i},1);
end
n=sum(n)

%% Save Data
save([basename '_BTColDec1']) %CY5

