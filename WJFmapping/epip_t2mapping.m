% clc
% clear
% close all
addpath(genpath('../../lib/'))
fid_name = '/home/shhcai_p1i/vnmrsys/studies/s_2021011103/T2_20210111_epip_02.fid';
[RE,IM,NY,NB,NL,HDR] = load_fid(fid_name);
K_space_data=RE+1i*IM;
param = getpar(fid_name,'nv','np','ns','etl','lro','pro','lpe','ppe','te','ti','nread','nphase','epiref_type','navigator','nnav','nf','WJF_sign','flip1','nseg','image','tr');
num_phase = param.nphase;nread = param.nread/2;num_slice = param.ns;tes = param.te;
num_te = length(param.te);num_ti = length(param.ti);
num_seg = param.nseg;etl = param.etl;nf = param.nf;
if param.navigator == 'y' && param.nnav~=0
    nnav = param.nnav+1;
else
    nnav =0;
end
num_phase = num_phase+nnav;
reconI = zeros(nread,num_phase-nnav,num_te);
reconK = zeros(nread,num_phase-nnav,num_te);
figure()
for loopi  = 1:num_slice
    for loopj = 1:num_te
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
%         subplot(2,2,loopj)
%         imagesc(abs(WJG_I));colormap jet
        reconI(:,:,loopj) = abs(WJG_I);
    end
    %% mapping
    minimum = 0;
    maxT2 = 0.15;
    MASK = reconI(:,:,1)>(max(max(abs(reconI(:,:,1))))/10);
    [Map thresh] = T2Map(reconI, tes, minimum, maxT2, 0 );
    subplot(2,5,loopi)
    imagesc(Map(:,:,2).*MASK,[minimum,maxT2]);colormap jet;axis off
end

