%for M0 fitting
%angus 2019.5.8
function [Map] = M0Map(indata, tes, minimum, T2map)
 % tes is a vector of the te values in msec.
 % indata is assumed to be dimensioned as X * Y * te - i.e., same slice, multiple te
 % T2map
rows = size(indata,1);
cols = size(indata,2);
N    = size(indata,3);
ly   = length(tes);
tes = reshape(tes, ly, 1);  % must be a columns vector
 
Map = zeros(rows,cols);       
indata = abs(indata);
 % Threshold images before fitting
thresh = sum( (abs(indata(:,:,:)) > minimum), 3) > (N-1);
inData = abs(indata);   % in case user has input complex data
for r=1:rows
    for c=1:cols
        if thresh(r,c) 
            ydata = inData(r,c,:);
            ydata = reshape(ydata,N,1);
            T2data = T2map(r,c);
            decay =  exp(-tes./(T2data+eps));
            lfit = ydata\decay;      
            Map(r,c,1) = 1./lfit(1);   
        end
    end
end
 
 






