% converting image library into polar coordinates (r/phi) representation 
%
% part of the diffractive image classification with subwavelength resolution project
% (c) 2021 V. Podolskiy - University of Massachusetts Lowell

% see A. Ghosh, et.al., ACS Photonics (2021) for more details 

clear
rng(1111)

noiseAdd=0; blurRad=2; % added noise parameters

suf='.L';
fldrIni=['./theory.images'];

%setting up the center of the CCD image
rMax=110; 
nR=100; nPhi=90; 
x0=68; y0=45;

% create a list of all the image files in the folder
listing=dir(fldrIni); 
fdata=cell((length(listing)-2)*10,3); 
ifn=1; 

for il=1:length(listing) 
    dirCur=listing(il).name; 
    if ~(strcmp(dirCur,'.')||strcmp(dirCur,'..'))
        listCur=dir([fldrIni,'\',dirCur]); 
        for ilc=1:length(listCur)
            if ~(listCur(ilc).isdir)
                fdata{ifn,1}=[fldrIni,'/',dirCur,'/',listCur(ilc).name]; 
                fdata{ifn,2}=dirCur; 
                ifn=ifn+1; 
            end 
        end 
    end
end 


% read and process images
% setup r/phi mesh
[x2,y2]=meshgrid((1:137)-x0,(1:182)-y0); 
r2=sqrt(x2.^2+y2.^2); 
phi2=atan2(y2,x2); 

r1Mesh=linspace(0,rMax,nR); 
phi1Mesh=linspace(0,pi,nPhi); 
[rMesh,phiMesh]=meshgrid(r1Mesh, phi1Mesh); 
xMesh=rMesh.*cos(phiMesh); yMesh=rMesh.*sin(phiMesh); 

for ifn=1:length(fdata)
    %read the images; flip them to mimic experimental conditions
    img=double(flipud(imread(fdata{ifn,1}))); 
    src1=img;
    sz1=size(src1); 
    mx1=max(src1(:)); 
    
    %add noise to the images
    if noiseAdd>0
        noiseMat=rand(sz1(1),sz1(2)); 
        noiseMat(noiseMat<1-noiseAdd/blurRad.^2)=0; 
        noiseMat(noiseMat>0)=1;
        noiseMat=imgaussfilt(noiseMat,blurRad); 
        noiseMat=noiseMat/max(noiseMat(:)); 
    else
        noiseMat=0*src1; 
    end 
    src1=src1+mx1*noiseMat; 

    src1(src1>mx1)=mx1; 
    src1(r2>rMax)=0; 
    src1(phi2<0)=0; 

    % remesh the data in r/phi mesh
    srcMesh=interp2(x2,y2,src1,xMesh,yMesh, 'makima',0);
    fdata{ifn,3}=srcMesh;
    
end 
    
save('./rphiImageData.mat', 'rMax','x2','y2','r2','phi2',...
    'r1Mesh','phi1Mesh','rMesh','phiMesh','xMesh','yMesh',...
    'fdata')

