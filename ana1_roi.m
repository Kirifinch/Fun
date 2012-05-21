% 21.2.2012

% to analyze the quadruple birds

% 5 files had been minsnamed by Mike and were corrected after this program had
% been run on them

% isolated=0 means isolate

clear
clf reset
close all

addpath('/Users/kirill/Documents/LAB_stuff/PET/SPM/spm8')

for filenum=1:12
    switch filenum
        case 1;
            file='/Users/kirill/Documents/LAB_stuff/Matlab_preprocessing_progs_from_Henning/analysis1/data/100311_2450_sumlast9_crop_2x';
            beakup=0;
            isolated=0;
            stimulated=0; % 1 means that Kirill played songs for 20 min each 30 s from tape before the scan
        case 2;
            file='/Users/kirill/Documents/LAB_stuff/Matlab_preprocessing_progs_from_Henning/analysis1/data/100311_2634_sumlast9_crop_2x';
            beakup=0;
            isolated=1;
            stimulated=0;
        case 3;
            file='/Users/kirill/Documents/LAB_stuff/Matlab_preprocessing_progs_from_Henning/analysis1/data/100311_2640_sumlast9_crop_2x';
            beakup=1;
            isolated=1;
            stimulated=1;
        case 4;
            file='/Users/kirill/Documents/LAB_stuff/Matlab_preprocessing_progs_from_Henning/analysis1/data/100311_2646_sumlast9_crop_2x';
            beakup=1;
            isolated=1;
            stimulated=1;
        case 5;
            file='/Users/kirill/Documents/LAB_stuff/Matlab_preprocessing_progs_from_Henning/analysis1/data/100511_2450_sumlast9_crop_2x';
            beakup=1;
            isolated=0;
            stimulated=1;
        case 6;
            file='/Users/kirill/Documents/LAB_stuff/Matlab_preprocessing_progs_from_Henning/analysis1/data/100511_2634_sumlast9_crop_2x';
            beakup=1;
            isolated=1;
            stimulated=1;
        case 7;
            file='/Users/kirill/Documents/LAB_stuff/Matlab_preprocessing_progs_from_Henning/analysis1/data/100511_2640_sumlast9_crop_2x';
            beakup=0;
            isolated=1;
            stimulated=0;
        case 8;
            file='/Users/kirill/Documents/LAB_stuff/Matlab_preprocessing_progs_from_Henning/analysis1/data/100511_2646_sumlast9_crop_2x';
            beakup=0;
            isolated=1;
            stimulated=0;
        case 9;
            file='/Users/kirill/Documents/LAB_stuff/Matlab_preprocessing_progs_from_Henning/analysis1/data/121411_1375_sumlast9_crop_2x';
            beakup=0;
            isolated=0;
            stimulated=1;
        case 10;
            file='/Users/kirill/Documents/LAB_stuff/Matlab_preprocessing_progs_from_Henning/analysis1/data/121411_2703_sumlast9_crop_2x';
            beakup=0;
            isolated=0;
            stimulated=1;
        case 11;
            file='C/Users/kirill/Documents/LAB_stuff/Matlab_preprocessing_progs_from_Henning/analysis1/data/121411_2759_sumlast9_crop_2x';
            beakup=1;
            isolated=0;
            stimulated=0;
        case 12;
            file='/Users/kirill/Documents/LAB_stuff/Matlab_preprocessing_progs_from_Henning/analysis1/data/121411_2762_sumlast9_crop_2x';
            beakup=1;
            isolated=0;
            stimulated=0;
        otherwise; stop
    end
    
    % do not change anything after this: %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    numberformat='int16'; endian='b';
    
    hdr=read_analyze([file '.hdr']);
    xdim=hdr.dime.dim(2);
    ydim=hdr.dime.dim(3);
    zdim=hdr.dime.dim(4);
    
    if xdim~=40 || ydim~=40 || zdim~=40; stop; end
    
    fid = fopen([file '.img'], 'rb');
    x = fread(fid, [xdim * ydim* zdim],numberformat,0,endian); % if Mike saves as nifti
    fclose(fid);
    x=reshape(x,xdim,ydim,zdim);
    
    
    if beakup; x=x(:,:,end:-1:1);  x=x(end:-1:1,:,:); end
    
    figure(1)
    subplot(2,2,1)
    imagesc(x(:,:,round(end/2)));
    
    x_cut=x; x_cut(:,:,1:8)=0; x_cut(:,:,32:40)=0; % prevent neck activation
    maxind=find(x_cut==max(max(max(x_cut)))); % compared with mricron, this is 40-x,y,z
    [xind,yind,zind]=ind2sub(size(x),maxind);
    
    subplot(2,2,2)
    imagesc(squeeze(x(xind,:,:)));
    axis square
    subplot(2,2,3)
    imagesc(squeeze(x(:,yind,:)));
    axis square
    subplot(2,2,4)
    imagesc(squeeze(x(:,:,zind)));
    axis square
    
    isolatedvec(filenum)=isolated;
    stimulatedvec(filenum)=stimulated;
    
    voi1=x(xind-1:xind+1,yind-1:yind+1,zind-1:zind+1); % based on looking at brain and best
    %voi1=x(xind-2:xind+2,yind-2:yind+2,zind-2:zind+2); % also good
    %voi1=x(xind,yind,zind); % high variability
    
    
    %  voi2=x(xind-1:xind+1,yind-1:yind+1,zind+3:zind+5); % based on looking at brain
    voi2=x(xind-1:xind+1,yind-1:yind+1,zind+4:zind+6); % best
    % voi2=x(xind-1:xind+1,yind-2:yind,zind+4:zind+6); % based on looking at brain with lower cb
    
    voi1_mean(filenum)=mean(mean(mean(voi1)));
    voi1_max(filenum)=max(max(max(voi1)));
    voi2_mean(filenum)=mean(mean(mean(voi2)));
    voi2_max(filenum)=max(max(max(voi2)));
    
    disp(' ')
    disp(num2str(filenum))
    disp(num2str(file))
    
    tmp=[xind,yind,zind]; disp(['Coordinates of max: ', num2str(tmp)])
    disp(['max1 = ' num2str(voi1_max(filenum)) ' mean1 = ' num2str(voi1_mean(filenum))])
    disp(['max2 = ' num2str(voi2_max(filenum)) ' mean2 = ' num2str(voi2_mean(filenum))])
    ratio=voi1_mean./voi2_mean;
    disp(['Ratio = ' num2str(ratio(filenum))])
    
    % pause
    
end

% isolated vs. colony %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
subplot(2,2,1)
title('Isolated vs. colony')
plot(ratio,'*')
subplot(2,2,2)
plot(ratio(isolatedvec==0),'k*')
hold on
plot(ratio(isolatedvec==1),'r*')
subplot(2,1,2)
plot(ratio,'*')
hold on
tmp=1:12;
plot(tmp(find(isolatedvec==1)),ratio(isolatedvec==1),'r*')

ratio_nonisolated=ratio(isolatedvec==0)
ratio_isolated=ratio(isolatedvec==1)

% stimulated vs. non-stimulated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
subplot(2,2,1)
title('Stimulated vs. non-stimulated')
plot(ratio,'*')
subplot(2,2,2)
plot(ratio(stimulatedvec==0),'k*')
hold on
plot(ratio(stimulatedvec==1),'r*')
subplot(2,1,2)
plot(ratio,'*')
hold on
tmp=1:12;
plot(tmp(find(stimulatedvec==1)),ratio(stimulatedvec==1),'r*')

ratio_nonstimulated=ratio(stimulatedvec==0)
ratio_stimulated=ratio(stimulatedvec==1)

% statistical testing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% only pairwise scanned birds:
x_stim0=ratio([7,8,1,2]); % 
x_stim1=ratio([3,4,5,6]); % 9,10
[h,p,ci,stats] = ttest(x_stim0,x_stim1);
p

% all male birds, non-paired t-test:
x_stim0=ratio([7,8,1,2,11,12]);  
x_stim1=ratio([3,4,5,6,9,10]); 
[h,p,ci,stats] = ttest2(x_stim0,x_stim1);
p

% all male birds, non-paired t-test:
x_iso0=ratio(isolatedvec==0);
x_iso1=ratio(isolatedvec==1);
[h,p,ci,stats] = ttest2(x_iso0,x_iso1);
p



