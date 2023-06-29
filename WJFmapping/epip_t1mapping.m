clc
clear
close all
addpath(genpath('../../lib'))
% addpath(genpath('D:\Users\How\Desktop\angus\T2star\data\lib\'))
% addpath(genpath('D:\Users\How\Desktop\angus\XMU_weekrep\weekrep-2020.7.24\deghost\niubi'))
fid_name = '/home/shhcai_p1i/vnmrsys/studies/s_2021011103/T1_20210111_epip_03.fid';
[RE,IM,NY,NB,NL,HDR] = load_fid(fid_name);
K_space_data=RE+1i*IM;
param = getpar(fid_name,'nv','np','ns','etl','lro','pro','lpe','ppe','te','ti','nread','nphase','epiref_type','navigator','nnav','nf','WJF_sign','flip1','nseg','image','tr');
num_phase = param.nphase;nread = param.nread/2;num_slice = param.ns;tis = param.ti;
num_te = length(param.te);num_ti = length(param.ti);
num_seg = param.nseg;num_etl = param.etl;nf = param.nf;
if param.navigator == 'y' && param.nnav~=0
    nnav = param.nnav+1;
else
    nnav =0;
end
num_phase = num_phase+nnav;
reconI = zeros(nread,num_phase-nnav,num_te);
reconK = zeros(nread,num_phase-nnav,num_te);
figure();
referenceT1 = zeros(nread,num_phase-nnav,num_slice);
for loopi  =3:num_slice
    for loopj = 1:num_ti
        WJG_K = zeros(nread,num_phase-nnav);
        phi_start = (loopj-1)*num_slice*num_phase;
        ph_sel = (loopj-1)*num_slice*num_seg+loopi:num_slice:(loopj-1)*num_slice*num_seg+num_slice*num_seg;
        tempK = K_space_data(:,ph_sel);
       for segi = 1:param.nseg
            segmentK= tempK(:,segi);
            segmentK = reshape(segmentK,[nread,param.etl+nnav]);
            segmentK(:,2:2:end) = flipud(segmentK(:,2:2:end));         
            segmentK = segmentK(:,nnav+1:end);
            segmentI = fftshift(fft2(ifftshift(segmentK)));  
            WJG_K(:,segi:param.nseg:end) = segmentK;
            WJG_I = fftshift(fft2(ifftshift(WJG_K)));  
       end
         %% deghost
         if (loopj==1)
            [kappaEstSimplex,phiEstSimplex]=WJG_deghostMS(WJG_I,num_seg);
         end
        finalK = zeros(nread,num_phase-nnav);
        for loop_seg = 1:num_seg
            segmentK = WJG_K(:,loop_seg:num_seg:end);
            kspCorrRefLessSimplex = applyFirstOrderPhaseCorr( segmentK, kappaEstSimplex, phiEstSimplex); 
            finalK(:,loop_seg:num_seg:end) = kspCorrRefLessSimplex;
        end
        WJG_I = fftshift(fft2(ifftshift(finalK)));
        subplot(2,2,loopj)
        imagesc(abs(WJG_I),[0,3e4]);colormap jet
        reconI(:,:,loopj) = abs(WJG_I);
    end
     %% mapping
    minimum = 0;
    maxT1 = 2;
    MASK = reconI(:,:,1)>(max(max(abs(reconI(:,:,1))))/10);
    [Map] = T1Map(reconI, tis, maxT1);
    referenceT1(:,:,loopi) = Map.*MASK;
    subplot(2,5,loopi)
    imagesc(Map.*MASK,[minimum,maxT1]);colormap jet;%colorbar
  end



