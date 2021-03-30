function [fullTbl,titles,fullFeatures,fileList] = rphi2tableFun(m1,j1,fname)
%RPHI2TABLEFUN calculates the Bessel transform based on r/phi data

load(fname,'rMesh','phiMesh','fdata','phi1Mesh','r1Mesh','rMax');

dr=rMesh(1,2)-rMesh(1,1); 
dphi=phiMesh(2,1)-phiMesh(1,1); 

% starting bessel transforms
[m2,j2]=meshgrid(m1,j1); 
jz2=bjZeros(m2,j2); 
j12=besselj(m2+1,jz2).^2; 
amj=0*m2; 

fullFeatures=zeros(length(fdata),numel(amj)); 
fullLabels=cell(length(fdata),1);
for ifn=1:length(fdata)
    fullLabels{ifn}=fdata{ifn,2}; 
    srcMesh=fdata{ifn,3}; 
    
integrant=srcMesh.*rMesh;

for im=1:length(m1)
    mC=m1(im); 
    cosC=cos(mC*phi1Mesh); 
    alpC=jz2(:,im); 
    bjC=besselj(mC,alpC*r1Mesh/rMax).';
    amj(:,im)=cosC*integrant*bjC; 
end 
amj=4*dphi*dr*amj/rMax./j12; 
fullFeatures(ifn,:)=amj(:); 

end
tmp=array2table(fullFeatures); 
fullTbl=[tmp,table(fullLabels)];
titles=cell(size(fullFeatures,2),1); 
for ifn=1:length(titles)
    titles{ifn}=['fullFeatures',num2str(ifn)];
end 

fileList=cell(size(fdata,1),2); 
for iln=1:size(fileList,1)
    fileList{iln,1}=fdata{iln,1}; 
end 

end

