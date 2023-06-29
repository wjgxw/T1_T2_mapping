%% gems
clc
clear
% close all
addpath('lib')
fid_name = 'data/gems_01_4T2star.fid';
[RE,IM,NY,NB,NL,HDR] = load_fid(fid_name);
param = getpar(fid_name,'te');
te_num = NB;             %the number of TEs
indata = zeros(NY,NY,te_num);
tes = param.te;
figure;
for loopi = 1:te_num
    WJG_K = RE(:,(loopi-1)*NY+1:loopi*NY)+1i*IM(:,(loopi-1)*NY+1:loopi*NY);
    WJG_I = fftshift(ifftn(fftshift(WJG_K)));
    indata(:,:,loopi) = abs(WJG_I);  
%     tes(loopi) = VCtl.TE;
    subplot(floor(sqrt(te_num))+1,floor(sqrt(te_num))+1,loopi);
    MASK1 = indata(:,:,1)>0.1;
    imagesc(abs(WJG_I.*MASK1));colormap jet
%     imagesc(angle(WJG_I.*MASK1),[-pi,pi]);colormap jet
%     set(gca,'position',[0 0 1 1]);
    
end
minimum = 0;
maxT2 = 0.3;
MASK1 = indata(:,:,1)>0.1;
[Map thresh] = T2Map(indata, tes, minimum, maxT2, 0 );
figure;imagesc(Map(:,:,2).*MASK1);colormap jet;colorbar
set(gca,'position',[0 0 1 1]);