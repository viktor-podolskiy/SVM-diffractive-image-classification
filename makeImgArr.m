% creating the image library for theoretical studies
%
% part of the diffractive image classification with subwavelength resolution project
% (c) 2021 V. Podolskiy - University of Massachusetts Lowell
 
% see A. Ghosh, et.al., ACS Photonics (2021) for more details 



rng(1111)

mMax=30; % parameter defining the spatial resolution of Fourier transforms,  
lam0=1; % operating wavelength for theoretical studies
nImgs=25; % number of images per object type in the library

% geometrical parameters
rMid=0.06; %mean radius
Lam=0.299; LamP=0.33; % grating periods
dr=0.01; % variation of radius
dxy=0.0025; % variation of position


% parameters used to emulate CCD used in experiments
lamO=0.532; % experimental wavelength
p0x=513; p0y=169; % normal incidence: 128
p1x=513; p1y=74;  % position of diffraction max
omg0=2*pi/lam0; 

dx=0.00125; % spatial mesh used for image generation
x1=(0:dx:12*Lam); 
y1=(0:dx:12*LamP); 

sname='./theory.images'; 

fnameLst={'I','C1','C9','L', ...
    'D','C4','S4','S9','S1',...
    'R','S9R','P'}; 

mkdir(sname); 

for sNum=(1:length(fnameLst))
    fname=fnameLst{sNum}; 
    if ~exist([sname,'/',fname],'dir')
        mkdir([sname,'/',fname])
    end 

    for iiter=1:nImgs
        
        [x2,y2]=meshgrid(x1,y1); 

        z2=0*x2; 

        dhM=3*LamP/4; 

        % blank, opaque sample canvas
        ix=(1:11); 
        iy=(1:11); 
        [ix2,iy2]=meshgrid(ix,iy); 
        iz2=1+0*ix2; 

        % creating hole patterns according to experimental geometries
        switch sNum
            case 2
                xc=6; yc=6;
                iz2((ix2==xc) & (iy2==yc))=0; 
            case 3
                xc=6; yc=6;
                iz2(abs(ix2-xc)<=1 & abs(iy2-yc)<=1)=0;
            case 4
                ix=3; 
                for iy=6:10
                    iz2((ix2==ix) & (iy2==iy))=0; 
                end 
            case 5
                for ix=(6:-1:2) 
                    for iy=12-(2*(6-ix)+[1 2])
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end 
            case 6 
                for ix=(1:2)+5
                    for iy=[6 7]
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end 
            case 7 
                for ix=[8 9]
                    for iy=[8 9]
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end 
            case 8
                for ix=[7 8 9]
                    for iy=[7 8 9]
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end 
            case 9 
                ix=8; iy=8;
                iz2((ix2==ix) & (iy2==iy))=0; 
            case 10 
                for ix=(4:6)
                    for iy=[3 4 ]
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end 
            case 11
                for ix=[7 8 9]
                    for iy=[7 8 9]
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end 
                for ix=(4:6)
                    for iy=[3 4 ]
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end 
            case 12
                for ix=[1 2]
                    for iy=[1 2 10 11]
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end 
                for ix=[3 4]
                    for iy=(3:5)
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end 
                for ix=4
                    for iy=[1 2 3 5 8 9]
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end 
                for ix=5
                    for iy=[1 2 3 5 6 8 9]
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end 
                for ix=6
                    for iy=[2 3 5 6]
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end 
                for ix=[9 10]
                    for iy=[1 2 6 7]
                        iz2((ix2==ix) & (iy2==iy))=0; 
                    end 
                end
                ix=11; iy=11; 
                iz2((ix2==ix) & (iy2==iy))=0; 
        end 


        % distribute the openings with randomized positions and sizes
        for ii=1:numel(iz2)
            if iz2(ii)>0
                xCur=ix2(ii)*Lam+2*dxy*(rand-0.5); 
                yCur=iy2(ii)*LamP+2*dxy*(rand-0.5);
                rCur=rMid+2*dr*(rand-0.5);

                z2(((x2-xCur).^2+(y2-yCur).^2)<=rCur.^2)=1; 
            end 
        end 

        % plot the configuration 
        figure(1)
        clf
        subplot(1,2,1)
        surf(x2,y2,z2,'EdgeColor','none')
        colormap gray
        xlim([0 max(x2(:))])
        ylim([0 max(y2(:))])
        xlabel('x')
        ylabel('y')
        box on 
        view(2)
%         keyboard



        % perform Fourier transform to emulate Fourier-microscopy
        % also perform inverse Fourier transfor to assess the adequacy of
        % mMax parameter used in the calculations
        dk=2*pi/Lam/(p0y-p1y); 
        kArr=((-mMax-.5)*omg0:dk:(mMax+.5)*omg0);
        [k2x,k2y]=meshgrid(kArr, kArr); 

        dx=2*pi/(max(k2x(:))-min(k2x(:))); 
        dy=2*pi/(max(k2y(:))-min(k2y(:))); 
        xArr=(1:length(k2x))*dx; 
        yArr=(1:length(k2y))*dy; 

        [xF2,yF2]=meshgrid(xArr,yArr); 

        FRx=exp(-1i*(k2x.').*xF2)*dx/2/pi; FRy=exp(-1i*(k2y.').*yF2)*dy/2/pi; 
        FFx=exp(1i*(k2x).*xF2.')*dk; FFy=exp(1i*(k2y).*yF2.')*dk; 
        zF=interp2(x2,y2,z2,xF2,yF2,'nearest',0); 
        Iij=FFy*(zF)*FFx; 
        Zij=FRy*(Iij)*FRx; 

        % spectral shift and NA-cutoff 
        kx0=omg0*sind(50); ky0=0; suf='L'; %Left incidence direction

        nMax=1.5; noiseL=0.1;  
        dkc=sqrt((k2x-kx0).^2+(k2y-ky0).^2); 
        Iij(dkc>nMax*omg0)=noiseL; 
        Iij(abs(Iij)<noiseL)=noiseL; 

        [pFx,pFy]=meshgrid(1:1024,1:255); 

        p2x=-k2x/dk*lam0/lamO+p0y; 
        p2y=-k2y/dk*lam0/lamO+p0x; 
        IFij=interp2(p2y.',p2x.',Iij.',pFx,pFy,'spline',noiseL); 

        IFij=abs(IFij); 
        IFij=IFij/max(max(IFij)); 

        % plot inverse FT image for verification purposes
        figure(1)
        subplot(1,2,2)
        surf(xF2,yF2,abs(Zij).^2,'EdgeColor','none')
        colormap gray
        view(2)
        xlim([0 max(x2(:))]) 
        ylim([0 max(y2(:))])
        xlabel('x')
        ylabel('y')
        drawnow 

        %plot initial estimate of Fourier transform
        figure(2)
        clf
        imagesc(IFij)
        view(2)
        colormap jet
        drawnow

        %rescale Fourier spectrum, remesh it to mimic CCD used in experimental setup, and implement lens aberration 
        IFij=1024*IFij; 

        [Px,Py]=meshgrid((1:1024),(1:255)); 
        kp2x=dk*(Px-p0x); kp2y=dk*(Py-p0y)+ky0; %L 

        [k2ix,k2iy]=meshgrid(...
            unique([(-15:dk:15),15]),...
            unique([(-30:dk:10),15]) ); 

        DatIJ=interp2(kp2x,kp2y,IFij,k2ix,k2iy,'cubic'); 
        cutLn=25;
        DatIJ(isnan(DatIJ))=cutLn; 

        k2iyp=k2iy.*(1-0.005/21*k2ix.^2); %mimic spherical aberration -L try

        figure(450)


        surf(k2ix,k2iyp,log(abs(DatIJ)), 'LineStyle','none')
        view(2)
        xlim([-15 15])
        ylim([-30 10])
        colormap gray
        axis off
        clrLim=[3 5.5]; 
        caxis(clrLim)
        fnameS=[sname,'/',fname,'/',suf,num2str(randi(1e6)),'.tif']; 

        % image processing
        datFlt=log(abs(DatIJ));
        datFlt(datFlt<clrLim(1))=clrLim(1); 
        datFlt(datFlt>clrLim(2))=clrLim(2); 

        dat16=uint16(65535*(datFlt-clrLim(1))/(clrLim(2)-clrLim(1))); 

%         keyboard 
        imwrite(dat16,fnameS); % save the image

    end 
end 
