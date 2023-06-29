% % T1 mapping calculation with the experimental data
% % created by Angus, wjtcw@hotmail.com
% % 2022.6.27
clc
clear
close all
addpath(genpath('lib'))
% sems
fid_name = 'data\sems_05.fid';
[RE,IM,NY,NB,NL,HDR] = load_fid(fid_name);
param = getpar(fid_name,'ti','ns');
ti = param.ti;
ti_num = length(ti);
indata = zeros(NY,NY,ti_num);
ns = param.ns;
for slice = 1:ns
    for loopi = 1:ti_num
        WJG_K = RE(:,ns*(loopi-1)+slice:ns*ti_num:end)+1i*IM(:,ns*(loopi-1)+slice:ns*ti_num:end);
        WJG_I = fftshift(ifftn(fftshift(WJG_K)));
        indata(:,:,loopi) = abs(WJG_I);  
        subplot(floor(sqrt(ti_num))+1,floor(sqrt(ti_num))+1,loopi);
        MASK1 = indata(:,:,1)>0.2;
        imshow(abs(WJG_I.*MASK1),[0,5]);colormap jet;
    %     set(gca,'position',[0 0 1 1]);    
    end
    minimum = 0;
    maxT1 = 3;
    MASK = indata(:,:,1)>0.1;
    [Map] = T1Map2(indata, ti, maxT1,MASK);
 
    figure;imshow(Map.*MASK,[]);colormap jet;colorbar
end
% % set(gca,'position',[0 0 1 1]);