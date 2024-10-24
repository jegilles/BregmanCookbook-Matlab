function u=ITV_NB_Deconvolution_3D(f,A,mu,lambda,Niter,typeKernel)

%==========================================================================
% function u=ITV_NB_Deconvolution_3D(f,A,mu,lambda,Niter,typeKernel)
%
% Anisotropic 3D TV Non-Blind Deconvolution
% Version:
% -v1.2 - 06/23/2013
% -v1.0 - 02/29/2012
%
% This function performs the minimization of the Isotropic TV - L2 
% Rudin-Osher-Fatemi model with a convolution kernel for 3D data by 
% Split Bregman Iterations
%
% f = 3D datacube
% A = 3D convolution kernel
% mu = regularization parameter
% lambda = fidelity factor for the split variables
% Niter = maximum number of iterations
% typeKernel: 0 <=> A is considered in the spatial domain
%             1 <=> A is considered in the Fourier domain and must be of
%             the same size of the image with frequency (0,0,0) at location
%             (1,1,1)
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
%
%==========================================================================

[M,N,P]=size(f);

f=double(f);
dx=zeros(M,N,P);
dy=zeros(M,N,P);
dz=zeros(M,N,P);
bx=zeros(M,N,P);
by=zeros(M,N,P);
bz=zeros(M,N,P);
u=f;
Z=zeros(M,N,P);

%Fourier mask + initialization of Fourier constant quantities
if typeKernel==0
    Mask=zeros(M,N,P);
    [H,L,O]=size(A);
    Mask([end+1-floor(H/2):end,1:ceil(H/2)],[end+1-floor(L/2):end,1:ceil(L/2)],[end+1-floor(O/2):end,1:ceil(O/2)]) = A;
    FMask=fftn(Mask);
else
    FMask = A;
end

%Fourier Laplacian mask initialization
D = zeros(M,N,P);
Lap=zeros(3,3,3);
Lap(:,:,1)=[0,0,0;0,1,0;0,0,0];
Lap(:,:,2)=[0,1,0;1,-4,1;0,1,0];
Lap(:,:,3)=[0,0,0;0,1,0;0,0,0];

D([end,1,2],[end,1,2],[end,1,2]) = Lap;
clear Lap;
FD=fftn(D);

%Fourier constant initialization
FW=((mu/lambda)*abs(FMask).^2-real(FD)).^-1;
FF=(mu/lambda)*conj(FMask).*fftn(f);

K=0;
err=sum(sum(sum(f.^2)));
tol=1e-3*err;
while ((err>tol) && (K<Niter)),
    K=K+1;
    
    tx=dx-bx;
    ty=dy-by;
    tz=dz-bz;
    
    up=u;
    %Update u
    u=real(ifftn(FW.*(FF-fftn(tx-tx(:,[1,1:N-1],:)+ty-ty([1,1:M-1],:,:)+tz-tz(:,:,[1,1:P-1])))));
    ux=u-u(:,[1,1:N-1],:);
    uy=u-u([1,1:M-1],:,:);
    uz=u-u(:,:,[1,1:P-1]);
        
    
    %Update dx
    tmpx=ux+bx;
    tmpy=uy+by;
    tmpz=uz+bz;
    
    s=sqrt(tmpx.^2+tmpy.^2+tmpz.^2);  
    thresh=max(Z,s-1/lambda)./max(1e-12,s);
    
    dx=thresh.*tmpx;
    dy=thresh.*tmpy;
    dz=thresh.*tmpz;
    
    %Update bx and by and bz
    bx=tmpx-dx;
    by=tmpy-dy;
    bz=tmpz-dz;
    
    err=sum(sum(sum((u-up).^2)));
end