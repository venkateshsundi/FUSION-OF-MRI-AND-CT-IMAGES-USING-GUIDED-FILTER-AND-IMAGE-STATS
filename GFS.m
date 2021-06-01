
close all;
clear all;
clc;
  
 % Parameters           
r=25; 
eps=2.1;
cov_wsize=5;

% source image 1
%  I1=double(imread('med256A.jpg'));
  I1=double(imread('source18A.tif'));

 if size(I1,3)==3
    I1=rgb2gray(I1);
 end

% source image 2

%  I2=double(imread('med256B.jpg'));
 I2=double(imread('source18B.tif'));

if size(I2,3)==3
 I2=rgb2gray(I2);
end
 
 I(:,:,1)=I1;
 I(:,:,2)=I2;
 
tic
% Base and detail layers seperation
A1= guidedfilter(double(I1), double(I2), r, eps);
B1=uint8(A1);
D1=double(I1)-A1;
C1=uint8(D1);
A2=guidedfilter(double(I2), double(I1), r, eps);
B2=uint8(A2);
D2=double(I2)-A2;
C2=uint8(D2);

 D(:,:,1)=D1;
 D(:,:,2)=D2;
 % Fusion rule
 xfused=GFS_fusion_rule(I,D,cov_wsize);
  toc
 FF=uint8(xfused);
% Display of images
figure, imshow(I1,[]);
figure, imshow(I2,[]);
figure, imshow(FF,[]);
[p q] = size(FF);
% Mean (F) or average pixel intensity (l)
mean_value = sum(sum(FF))/(p*q)
% Standard deviation (SD or r)
variance_value  = sqrt(sum(sum((double(FF)-mean_value).^2))/(p*q))
% Average gradient (AG)
FF = double(FF);
for i=1:p-1
    for j=1:q-1
        a1(i,j) = (FF(i,j)-FF(i+1,j)).^2 ;
        a2 (i,j) = (FF(i,j)-FF(i,j+1)).^2 ;
        a3 = sqrt(a1+a2)/((p)*(q));
    end
end
average_grad = sum(sum(a3))
% Mutual information (MI) or fusion factor
MIXF=mi(I1,FF); 
MIYF=mi(I2,FF);
Mutual_information = MIXF+MIYF
% Spatial frequency (SF)
for i=2:p
    for j=2:q
        a4 (i,j) = (FF(i,j)-FF(i,j-1)).^2 ;
        a5 (i,j) = (FF(i,j)-FF(i-1,j)).^2 ;
    end
end
rf = sqrt(sum(sum(a4))/((p)*(q)));
cf = sqrt(sum(sum(a5))/((p)*(q)));
Spatial_frequency = sqrt(rf^2+cf^2)    
%Fusion information score

[AX,AY] = gradient(I1);
GA = sqrt(AX.^2+AY.^2);
AA = atan(AY./AX);
GA(isnan(GA))=0;AA(isnan(AA))=0;

[BX,BY] = gradient(I2);
GB = sqrt(BX.^2+BY.^2);
AB = atan(BY./BX);
GB(isnan(GB))=0;AB(isnan(AB))=0;

[FX,FY] = gradient(FF);
GF = sqrt(FX.^2+FY.^2);
AF = atan(FY./FX);
GF(isnan(GF))=0;AF(isnan(AF))=0;
for i=1:464
    for j=1:464
        if(GA(i,j) > GF(i,j) )
            GAF (i,j)= GF(i,j) / GA(i,j);
        else
            GAF (i,j)= GA(i,j) / GF(i,j);
        end
    end
end
GAF(isnan(GAF))=0;
for i=1:464
    for j=1:464
        if(GB(i,j) > GF(i,j) )
            GBF (i,j)= GF(i,j) / GB(i,j);
        else
            GBF (i,j)= GB(i,j) / GF(i,j);
        end
    end
end
GBF(isnan(GBF))=0;
for i=1:464
    for j=1:464
        AAF (i,j)= abs(abs(AA(i,j) - AF(i,j)) - (pi/2))/(pi/2);
    end
end
for i=1:464
    for j=1:464
        ABF (i,j)= abs(abs(AB(i,j) - AF(i,j)) - (pi/2))/(pi/2);
    end
end
for i=1:464
    for j = 1:464
        QAFG =  7.9994/(1+exp(15*(GAF(i,j)-0.5)));
        QAFA =  8.9879/(1+exp(22*(AAF(i,j)-0.8)));
    end
end
for i=1:464
    for j = 1:464
        QBFG =  6.9994/(1+exp(15*(GBF(i,j)-0.5)));
        QBFA =  7.9879/(1+exp(22*(ABF(i,j)-0.8)));
    end
end
QAF = QAFG .* QAFA;
QBF = QBFG .* QBFA;
WA = GA.^2.5;
WB = GB.^2.5;
W = sum(sum(WA + WB));
QABF = sum(sum(QAF.*WA+QBF.*WB))/sum(sum(WA + WB));
FUSION_SCORE = QABF