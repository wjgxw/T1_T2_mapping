function [Map thresh] = T2Map(indata, tes, minimum, maxT2, display )
% Calculate 1/exponential decay rate for T2 (or T1 by ULFMRI)
%
% indata is assumed to be dimensioned as X * Y * te - i.e., same slice, multiple te
%
% tes is a vector of the te values in msec.
%
% minimum sets the minimum pixel intensity on the input data that will be fitted
%
% maxT2 impose limits on the fitted T2 values, as an outlier suppression
%
% display is set to 0 or 1 to either display the results, or not.
%
% I also created a MATLAB version - though I haven't looked recently at the
% code. This script assumes that you have gotten the data into a MATLAB matrix. 
% As I glance at it now, it is a bit rough, as it seems to have minimal useful comments (see below). 
% It looks to me like it was set up to look at 2D images. You would have to either modify the 
% indexing, or call this function in a loop (by slice location) to do this to a volume.
 
if nargin < 5
    display = false;
end
if nargin < 4
    maxT2 = 2000;
end
if nargin < 3
    minimum = 0;
end
 
if nargin < 2
    disp(sprintf('function [Map thresh] = %s( inData, tes, minimum, maximum, display )',mfilename) );
    disp('Image 1: rho   Image 2: T2   Image 3: R^2');
    return;
end
 
global stop;
rows = size(indata,1);
cols = size(indata,2);
N    = size(indata,3);
ly   = length(tes);
tes = reshape(tes, ly, 1);  % must be a columns vector
 
Map = zeros(rows,cols,2);   %SLL: there are two maps because...
Map(:,:,:) = 0;
tes = [ones(ly,1), tes];  % first columns is ones to allow fit of intercept (rho)
    
indata = abs(indata);
 
% Threshold images before fitting
thresh = ones(rows,cols);
thresh = sum( (abs(indata(:,:,:)) > minimum), 3) > (N-1);
 
inData = abs(indata);   % in case user has input complex data
 
stop = false;
SkippedPixels = 0;
 
h = waitbar(0, 'Fitting T2...','CreateCancelBtn','global stop; stop = true;' );
set(h,'Name','T2fit Progress');
tic;
 
minR2 = 1/maxT2;
AvgCutoff = minimum * sqrt(ly); % Operationally define a cutoff based on vector length
for r=1:rows
    if stop
        break;
    end
    for c=1:cols
        if stop
            break;
        end
        ydata = inData(r,c,:);
        ydata = reshape(ydata,N,1);
        if thresh(r,c) && ( mean(ydata) > AvgCutoff )
            lfit = tes \ log(ydata);            %lfit=-1/T2
            % res = glmfit(tes,log(ydata),[],'constant','on');
            if( lfit(2) < 0 && -lfit(2) > minR2 )
                Map(r,c,1) = lfit(1);
                Map(r,c,2) = -1/lfit(2);
                % Map(r,c,3) = stats(1);
            elseif(-lfit(2) < minR2)
                Map(r,c,1) = lfit(1);
                Map(r,c,2) = maxT2;
            else
                thresh(r,c) = 0;
                SkippedPixels = SkippedPixels + 1;
            end
        else
            thresh(r,c) = 0;
            SkippedPixels = SkippedPixels + 1;
        end
    end
    rate = toc/r;
    remaining = (rows-r) * rate;
    update_string = sprintf('%0.1f s/row for %d rows. Remaining ~ %d min %d s. %0.1f%% pixels used', ...
        rate, r, floor(remaining/60), round(mod(remaining,60)), 100*(1-SkippedPixels/(r*cols)) );
    waitbar(r/rows, h, update_string, 'CreateCancelBtn', 'stop = true;');
end
 
 
if exist('h')
    delete(h);
    clear h;
end
 
if stop
    disp('Calculations canceled by user.');
    return;
end
 
Map(:,:,1) = exp(Map(:,:,1));
 
if display
    clf;
    subplot(3,1,1);
    image(thresh,'CDataMapping','scaled');
    colorbar;
    axis image;
    subplot(3,1,2);
    image(abs(Map(:,:,1)),'CDataMapping','scaled');colormap jet;
    axis image;
    colorbar
    subplot(3,1,3);
    image(abs(Map(:,:,2)),'CDataMapping','scaled');colormap jet;
    axis image;
    colorbar
end