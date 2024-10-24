function u=Curvelet_NB_Deconvolution_3D(f,A,mu,lambda,Niter,typeKernel)

%==============================================================================
%
% function u=Curvelet_NB_Deconvolution_3D(f,A,mu,lambda,Niter,typeKernel)
%
% This function performs the 3D Nonblind Curvelet deconvolution by the analysis 
% approach.
% Version:
% -v1.2 - 06/23/2013
% -v1.0 - 03/02/2012
%
% Parameters:
% f: input 3D datacube
% A: input 3D blur kernel
% mu,lambda: regularization parameters
% Niter: maximum number of iterations
% typeKernel: 0 <=> A is considered in the spatial domain
%             1 <=> A is considered in the Fourier domain and must be of
%             the same size of the image with frequency (0,0,0) at location
%             (1,1,1)
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
%
%==============================================================================


%structures initialization
[M,N,P]=size(f);
U=fdct3d_forward(zeros(M,N,P));
d=U;
b=U;

%Fourier mask + initialization of Fourier constant quantities
if typeKernel==0
    Mask=zeros(M,N,P);
    [H,L,O]=size(A);
    Mask([end+1-floor(H/2):end,1:ceil(H/2)],[end+1-floor(L/2):end,1:ceil(L/2)],[end+1-floor(O/2):end,1:ceil(O/2)]) = A;
    FMask=fftn(Mask);
else
    FMask = A;
end


%Fourier constant initialization
FW=(mu*abs(FMask).^2+lambda).^-1;
FF=mu*conj(FMask).*fftn(f);
clear FD;
clear FMask;

u=f;
K=0;
%Bregman iterations
err=sum(sum(sum(f.^2)));
tol=1e-3*err;
while ((err>tol) && (K<Niter)),
    K=K+1;
    
    up=u;
    %update u
    tx=fdct3d_inverse(SubCurveletArray(d,b));
    u=real(ifftn(FW.*(FF+lambda*fftn(tx))));
    
    %update d
    U=AddCurveletArray(fdct3d_forward(u),b);
    d=ShrinkComplexCurvelet(U,1/lambda);
    
    %update b
    b=SubCurveletArray(U,d);
    
    err=sum(sum(sum((u-up).^2)));
end
