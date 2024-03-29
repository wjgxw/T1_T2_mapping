% A Robust Methodology for In Vivo T1 Mapping��2010�� MRM
function [Map] = T1Map2(indata, ti, maxT1,MASK)
[row,col,~] = size(indata);
extra.T1Vec = linspace(1,maxT1*1000,100);  % Initial grid points for the T1 search
extra.tVec = ti*1000; % Inversion times (TIs) considered
nlsS = getNLSStruct(extra);
Map = zeros(row,col);
parfor loopi = 1:row
    for loopj = 1:col
        if MASK(loopi,loopj)==1
            S = indata(loopi,loopj,:);
            [T1Est, ~, ~, ~] = rdNlsPr(abs(S), nlsS);
            Map(loopi,loopj) = T1Est;
        end
    end
end