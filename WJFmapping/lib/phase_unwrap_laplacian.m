function [unwrappedPhase]  = phase_unwrap_laplacian(dataPhase, Params, RefVox, zp)
% [unwrappedPhase]  = phase_unwrap_laplacian(dataPhase, BrainMask, Params, RefVox, zp)
%================================
% HEADER
%================================
%
% Name:             phase_unwrap_laplacian.m
% Author:           Xu Li, PhD
% Email:            xuli@mri.jhu.edu
%
%================================
% PURPOSE
% unwrap phase according to the laplacian method
% i.e. del2(theta) = cos(theta)*del2(sin(theta)) + sin(theta)*del2(cos(theta))
%
%  Updated: 2012-08-17: Xu Li
%  Updated: 2014-01-03: JvB Added waitbars

textWaitbar = 'Performing phase unwrapping';
% multiWaitbar(textWaitbar, 0, 'Color', 'b' );

warning off all

%% zero padding
Nx = (Params.sizeVol(2));
Ny = (Params.sizeVol(1));
Nz = (Params.sizeVol(3));

if nargin < 3
    zp = 2;
    RefVox = [floor(Ny/2), floor(Nx/2), floor(Nz/2)];    
elseif nargin < 4
    zp = 2;
end

FOVx = Params.fov(2);
FOVy = Params.fov(1);
FOVz = Params.fov(3);

Nx_zp = Nx*zp;
Ny_zp = Ny*zp;
Nz_zp = Nz*zp;

% Nx_zp = 2^(round(log2(Nx*zp)));     % round to exponential of 2
% Ny_zp = 2^(round(log2(Ny*zp)));
% Nz_zp = 2^(round(log2(Nz*zp)));
% 
% if Nx_zp > Nx*zp                    % use zero padding ratio not larger than zp
%     Nx_zp = Nx*zp;
% end
% 
% if Ny_zp > Ny*zp
%     Ny_zp = Ny*zp;
% end
% 
% if Nz_zp > Nz*zp
%     Nz_zp = Nz*zp;
% end

zpx = Nx_zp/Nx;                     % final zero padding ratio
zpy = Ny_zp/Ny;
zpz = Nz_zp/Nz;

FOVx_zp = FOVx*zpx;
FOVy_zp = FOVy*zpy;
FOVz_zp = FOVz*zpz;

x1 = floor((Nx_zp - Nx)/2) + 1;     % index where to put the original data
x2 = x1 + Nx - 1;
y1 = floor((Ny_zp - Ny)/2) + 1;
y2 = y1 + Ny - 1;
z1 = floor((Nz_zp - Nz)/2) + 1;
z2 = z1 + Nz - 1;

Phasej = dataPhase;
if isa(dataPhase(1), 'single')
    dataPhase = zeros(Ny_zp, Nx_zp, Nz_zp, 'single');
else
    dataPhase = zeros(Ny_zp, Nx_zp, Nz_zp);    
end
% put the data in the center of the zero padded one
dataPhase(y1:y2, x1:x2, z1:z2) = Phasej;

%% calculate the k^2 in k space

dkx = 1/FOVx_zp;
dky = 1/FOVy_zp;
dkz = 1/FOVz_zp;

kx = linspace(-Nx_zp/2+1, Nx_zp/2, Nx_zp).*dkx;
ky = linspace(-Ny_zp/2+1, Ny_zp/2, Ny_zp).*dky;
kz = linspace(-Nz_zp/2+1, Nz_zp/2, Nz_zp).*dkz;

% Update
multiWaitbar(textWaitbar, 0.2, 'Color', 'b' );

if isa(dataPhase(1), 'single')
    % Form the KSq_Grid using a for loop
    KSq_Grid = zeros(Ny_zp, Nx_zp, Nz_zp, 'single');
    for kk = 1:Nz_zp
        for jj = 1:Nx_zp
            for ii = 1:Ny_zp
                KSq_Grid(ii, jj, kk) = (kx(jj)^2 + ky(ii)^2 + kz(kk)^2);
            end
        end
    end
else
    [KX_Grid, KY_Grid, KZ_Grid] = meshgrid(kx, ky, kz);  % mesh in k space
    KSq_Grid = (KX_Grid.^2 + KY_Grid.^2 + KZ_Grid.^2);

    clear KX_Grid KY_Grid KZ_Grid
end

% Update
multiWaitbar(textWaitbar, 0.45, 'Color', 'b' );

%% calculate 
ThetaSin = sin(dataPhase);
ThetaCos = cos(dataPhase);
clear dataPhase

laplacianSin = fftn(ThetaSin);
laplacianSin = fftshift(laplacianSin).*KSq_Grid;
laplacianSin = ifftshift(laplacianSin);
laplacianSin = ifftn(laplacianSin);
laplacianSin = real(laplacianSin);

% Update
multiWaitbar(textWaitbar, 0.60, 'Color', 'b' );

laplacianCos = fftn(ThetaCos);
laplacianCos = fftshift(laplacianCos).*KSq_Grid;
laplacianCos = ifftshift(laplacianCos);
laplacianCos = ifftn(laplacianCos);
laplacianCos = real(laplacianCos);

% Update
multiWaitbar(textWaitbar, 0.75, 'Color', 'b' );

A = laplacianSin.*ThetaCos - laplacianCos.*ThetaSin;
clear ThetaSin ThetaCos laplacianSin laplacianCos

% Update
multiWaitbar(textWaitbar, 0.85, 'Color', 'b' );

unwrappedPhase = fftshift(fftn(A))./KSq_Grid;
unwrappedPhase(isinf(unwrappedPhase)) = 0;
clear A KSq_Grid

% Update
multiWaitbar(textWaitbar, 0.95, 'Color', 'b' );

unwrappedPhase = real(ifftn(ifftshift(unwrappedPhase)));
unwrappedPhase(isnan(unwrappedPhase)) = 0;

%% shift to a known constant (find the largest phase region and select the middle point for reference)
unwrappedPhase = unwrappedPhase(y1:y2, x1:x2, z1:z2);

c = - unwrappedPhase(RefVox(1), RefVox(2), RefVox(3));

% Update
multiWaitbar(textWaitbar, 0.99, 'Color', 'b' );

unwrappedPhase = unwrappedPhase + c;        % phi' as in Schofield & Zhu, Optics Letters, 2003