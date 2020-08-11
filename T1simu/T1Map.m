%for T1 fitting
%angus 2019.11.16
%angus modified 2020.8.11
%a*(1-2*exp(-x/T1)) --> a*(1-b*exp(-x/T1))
function [Map] = T1Map(indata, tis,maxT1)
 % tis is a vector of the ti values in msec.
 % indata is assumed to be dimensioned as X * Y * ti - i.e., same slice, multiple ti
rows = size(indata,1);
cols = size(indata,2);
ly   = length(tis);
tis = reshape(tis, ly, 1);  % must be a columns vector
temp_data = indata(:,:,:);
a_max = max(temp_data(:))*2;
a_min = 0;
MASK1 = max(indata(:,:,:),[],3)>(a_max/20);
Map = zeros(rows,cols);  

eq=fittype('a*(1-b*exp(-x/T1))','coefficients',{'a','b','T1'});

parfor rowi= 1:rows
    for coli= 1:cols
        if (MASK1(rowi,coli))
           signal= squeeze(indata(rowi,coli,:));%;[9.855,3.307,2.532,6.296,9.375,11.76]';%
           [mindata,idx] = min(signal);
           temp_signa = signal;
           %% type1
            if mindata>=0 && idx>1
                    temp_signa(1:idx-1) = -signal(1:idx-1);
            end
            fo = fitoptions('method','NonlinearLeastSquares','StartPoint',[max(signal(:)),2,tis(idx)],...
                                     'Lower',[a_min,0,0.01],...
                                     'Upper',[a_max,2,maxT1],'MaxFunEvals',1e10,'MaxIter',1e4);
           [T1f1_1,diff1]=fit(tis,temp_signa,eq,fo);
           %% type2
           temp_signa(1:idx) = -signal(1:idx);
           fo = fitoptions('method','NonlinearLeastSquares','StartPoint',[0,2,tis(idx)],...
                                     'Lower',[a_min,0,0.01],...
                                     'Upper',[a_max,2,maxT1],'MaxFunEvals',1e10,'MaxIter',1e4);
           [T1f1_2,diff2]=fit(tis,temp_signa,eq,fo);
           if (diff1.rmse<diff2.rmse)
               Map(rowi,coli)=T1f1_1.T1;
           else
               Map(rowi,coli)=T1f1_2.T1;
           end 
        end
    end
end






