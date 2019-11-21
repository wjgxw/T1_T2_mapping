%T2 mapping calculation
clc
clear all
close all
row=128;
col = row;
te_num = 5;
indata = zeros(row,col,te_num);
tes = zeros(1,te_num);
for loopi = 1:te_num
    filename = ['Series',num2str(loopi)];
    load(filename)
    II = VImg.Mag;
    indata(:,:,loopi) = II;  
    tes(loopi) = VCtl.TE;
    subplot(2,3,loopi);
    imshow(II,[]);title(num2str(VCtl.TE))
end
minimum = 0;
maxT2 = 2.2;
[Map thresh] = T2Map(indata, tes, minimum, maxT2, 1 );
imshow(Map(:,:,2),[])