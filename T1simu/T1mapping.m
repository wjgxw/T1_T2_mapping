% % T1 mapping calculation with the experimental data
% % created by Angus, wjtcw@hotmail.com
% % 2019.10.11
% the data was simulated with MRiLab IR3D sequence
% the virtual imaging model can be found in /model
clc
clear
close all
% 
ti_num = 5;
tis = zeros(1,ti_num);
NY = 128;
indata = zeros(NY,NY,ti_num);
fid_all = dir('Series*');
for loopi = 1:ti_num
    filename = fid_all(loopi).name;
    load(filename)
    tis(loopi) = VCtl.TI;
    WJG_I = VImg.Mag;
    indata(:,:,loopi) = WJG_I;
    indata(:,:,loopi) = indata(:,:,loopi);
end
temp_data = indata(:,:,:);
a_max = max(temp_data(:))*2;
a_min = 0;
maxT1 = 3;
MASK1 = max(indata(:,:,:),[],3)>(a_max/20);
subplot(floor(sqrt(ti_num))+1,floor(sqrt(ti_num))+1,loopi+1);  imshow(MASK1)
%%%%%%%%%%%%
%% sort data with tis
[tis,idx] = sort(tis);
indata_final = zeros(size(indata));
indata(:,:,:) = indata(:,:,idx);
for loopi = 1:ti_num
    subplot(floor(sqrt(ti_num))+1,floor(sqrt(ti_num))+1,loopi);    
    imagesc(abs( indata(:,:,loopi)));colormap jet;axis off
    title(num2str(tis(loopi)))
end
maxT1 = 3;
TI=tis;
Map = T1Map(indata, tis, maxT1);
figure;imshow(Map.*MASK1,[0,maxT1]);colormap jet;colorbar




    
    