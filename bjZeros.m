function [bjzArr] = bjZeros(mArr,nArr)
%BJZEROS Returns zeros of bessel function; requres xl file with
%precalculated zeros

mathData=xlsread('bjz.xlsx'); 
mMath=(0:size(mathData,1)-1); 
nMath=(1:size(mathData,2)); 
[n2M,m2M]=meshgrid(nMath,mMath); 

bjzArr=interp2(n2M,m2M,mathData,nArr,mArr,'nearest'); 
end

