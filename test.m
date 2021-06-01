

clc
close all
clear all

FF = [1 1 2 1;1 1 1 2;1 1 1 1;1 1 1 1;1 1 1 1]
ag = 0;
FF = double(FF);
[p,q] = size(FF);
for i=1:p-1
    for j=1:q-1
        a1(i,j) = (FF(i,j)-FF(i+1,j)).^2 ;
        a2 (i,j) = (FF(i,j)-FF(i,j+1)).^2 ;
        a3 = sqrt(a1+a2)/(p*q);
    end
end
ag = sum(sum(a3))
