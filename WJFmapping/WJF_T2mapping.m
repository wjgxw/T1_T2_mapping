% % T2 mapping calculation with the experimental data
% % created by Angus, wjtcw@hotmail.com
% % 2018.12.18
clc
clear
% close all
% sems
addpath('lib')
fid_name = '/home/shhcai_p1i/vnmrsys/studies/s_2020121202/sems_01.fid';
[RE,IM,NY,NB,NL,HDR] = load_fid(fid_name);
param = getpar(fid_name,'te','ns','np','nv');
nread = param.np/2;
nphase = param.nv;
tes = param.te;     %can be found in file "parcpar"
te_num = length(tes);             %the number of TEs
indata = zeros(nread,nphase,te_num);

for slicei = 1:param.ns
    figure(1)
    for loopi = 1:te_num
        WJG_K = RE(:,slicei+(loopi-1)*param.ns:param.ns*te_num:end)+1i*IM(:,slicei+(loopi-1)*param.ns:param.ns*te_num:end);
        WJG_I = fftshift(ifftn(fftshift(WJG_K)));
        indata(:,:,loopi) = abs(WJG_I);  
        subplot(floor(sqrt(te_num))+1,floor(sqrt(te_num))+1,loopi);

        MASK1 = indata(:,:,1)>0.2;
        imshow(abs(WJG_I.*MASK1),[]);colormap jet;
    %     set(gca,'position',[0 0 1 1]);    
    end
    minimum = 0;
    maxT2 = 0.25;
    MASK = indata(:,:,1)>0.1;
    [Map thresh] = T2Map(indata, tes, minimum, maxT2, 0 );
    figure(3);imagesc(Map(:,:,2).*MASK,[minimum,maxT2]);colormap jet;colorbar
    T2map = Map(:,:,2).*MASK;
%     M0 = M0Map(indata, tes, minimum, T2map);
%     figure;imagesc(M0.*MASK,[]);colormap jet;colorbar
end
% % set(gca,'position',[0 0 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
