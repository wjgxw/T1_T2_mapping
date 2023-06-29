%% B0 mapping
%% gems
clc
clear
% close all
addpath(genpath('../lib/'))
fid_name = 'K:\oledT1\s_2021010901\B0_20210109_mgems_01.fid';
[RE,IM,NY,NB,NL,HDR] = load_fid(fid_name);
param = getpar(fid_name,'te','lro','nv','ns','ne','te','te2');
ne = param.ne;
ns = param.ns;
tes = param.te:param.te2:ne*param.te2;
te_num = length(tes);
indata = zeros(NY,NY,te_num);

fov=param.lro/100;
res=param.nv;
vox_size = fov/res;
figure(3);
indata = zeros(res,res,te_num);
for loop_ns = 7:ns
    for loop_ne = 1:ne
        WJG_K = RE(:,(loop_ns-1)*ne+loop_ne:ns*ne:end)+1i*IM(:,(loop_ns-1)*ne+loop_ne:ns*ne:end);
        WJG_I = fftshift(fft2(ifftshift(WJG_K)));
        Ph1 = angle(WJG_I);
        Mag1 = abs(WJG_I);     
        unwrappedPhaseImage = sunwrap(  Mag1.*exp(1i*Ph1), 0.01 );
%     unwrappedPhaseImage = unwrapLaplacian(Ph1, [res,res], [vox_size,vox_size]);
    indata(:,:,loop_ne) = unwrappedPhaseImage;
% subplot(floor(sqrt(te_num))+1,floor(sqrt(te_num))+1,loop_ne);
%     imshow(abs( unwrappedPhaseImage),[]);colormap jet 
        
    end
    B0 = (indata(:,:,3)-indata(:,:,2))/(tes(3)-tes(2)+eps)/(2*pi);
%     subplot(floor(sqrt(te_num))+1,floor(sqrt(te_num))+1,loop_ne+1);
    imshow(B0.*(Mag1>0.1),[-100,100]);colormap jet; 
end
imshow(indata(:,:,2)-indata(:,:,1),[]);colormap jet
