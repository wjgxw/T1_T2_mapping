% epip
% 2019.8.28

clc
clear
close all
addpath('/home/shhcai_p1i/Desktop/seq/UTE_Package_2July2012/mriepdev/recon_matlab/lib')
%% parameters
fid_name = '/home/shhcai_p1i/vnmrsys/studies/s_2020061001/epip_08.fid';

%% start read the data
[RE,IM,NY,NB,NL,HDR] = load_fid(fid_name);
param = getpar(fid_name,'nv','np','ns','lro','pro','lpe','ppe','te','nread','nphase','epiref_type','navigator','nnav');

if param.epiref_type=='manual'
    nphase = param.nphase;
    nread = param.nread/2;
    nslice = param.ns;
    nfigure = nslice;
    nnav = param.nnav+1;
    nphase = nphase+nnav;
    for slicei = 1:nfigure
        WJG_K = RE((slicei-1)*nread*nphase+1:slicei*nread*nphase)+1i*IM((slicei-1)*nread*nphase+1:slicei*nread*nphase);
        WJG_K = reshape(WJG_K,[nread,nphase]);
        WJG_K = WJG_K(:,nnav+1:end);
        WJG_K(:,2:2:end) = flipud(WJG_K(:,2:2:end));        
        WJG_I = fftshift(ifftn(fftshift(WJG_K)));   
%% nyquist ghost


%% image shift
        shift_read = floor(param.pro/(param.lro/nread));
        shift_phase = floor(param.ppe/(param.lpe/nphase));
%         WJG_I = circshift(WJG_I,[-shift_read,-shift_phase]);
        subplot(floor(sqrt(nfigure)+1),floor(sqrt(nfigure)+1),slicei)
        imagesc((abs(WJG_I)));colormap gray
    end    
elseif param.epiref_type=='single'
    nphase = param.nphase;
    nread = param.nread/2;
    nslice = param.ns;
    nfigure = nslice*2;
    nnav = param.nnav+1;
    nphase = nphase+nnav;
    for slicei = 1:nfigure
        WJG_K = RE((slicei-1)*nread*nphase+1:slicei*nread*nphase)+1i*IM((slicei-1)*nread*nphase+1:slicei*nread*nphase);
        WJG_K = reshape(WJG_K,[nread,nphase]);
        WJG_K = WJG_K(:,nnav+1:end);
        WJG_K(:,2:2:end) = flipud(WJG_K(:,2:2:end));        
        WJG_I = fftshift(ifftn(fftshift(WJG_K)));   
%% nyquist ghost


%% image shift
        shift_read = floor(param.pro/(param.lro/nread));
        shift_phase = floor(param.ppe/(param.lpe/nphase));
%         WJG_I = circshift(WJG_I,[-shift_read,-shift_phase]);
        subplot(floor(sqrt(nfigure)+1),floor(sqrt(nfigure)+1),slicei)
        imagesc((abs(WJG_I)));colormap gray
    end
elseif param.epiref_type=='triple'
    disp('Try another epiret_type')
elseif param.epiref_type=='fulltriple'  
    disp('Try another epiret_type')
end
