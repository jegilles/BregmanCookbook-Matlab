function u=L1_L2_2D(f,A,mu,lambda,Niter,typeKernel)

%=======================================================
% function u=L1_L2_2D(f,A,mu,lambda,Niter,typeKernel)
%
% L1-L2 image deconvolution
% Version:
% -v1.0 - 03/14/2016
%
% This function compute the solution (@ = convolution)
% u = arg min ||u||_1+0.5*mu||A@u-f||_2^2
%
% by using the Split Bregman Iteration. See below "classic
% parameter values"
%
% f: blurry image
% A: PSF (can be either in spatial or Fourier domain)
% mu: regularization coefficient
% lambda: "spliting" regularization coefficient
% Niter: maximum number of iteration
% typeKernel: 0=spatial PSF, 1=Fourier PSF
%
% Author: Jerome Gilles
% Institution: SDSU - Dept of Mathematics and Statistics
% email: jgilles@mail.sdsu.edu
% 
% Note: typically mu=1e7, lambda=10,
%       Niter=10 work well
%========================================================
[M,N]=size(f);
f=double(f);
f=(f-min(f(:)))/(max(f(:))-min(f(:)));

%Structures and constants initialization
f=double(f);
d=zeros(M,N);
b=zeros(M,N);
u=f;
Z=zeros(M,N);

K=0;

%Fourier mask + initialization of Fourier constant quantities
if typeKernel==0
    Mask=zeros(M,N);
    [H,L]=size(A);
    Mask([end+1-floor(H/2):end,1:ceil(H/2)],[end+1-floor(L/2):end,1:ceil(L/2)]) = A;
    FMask=fft2(Mask);
else
    FMask = A;
end

%Fourier constant initialization
FW=mu*abs(FMask).^2+lambda;
FW=FW.^-1;
FF=mu*conj(FMask).*fft2(f);

%Bregman iterations
err=norm(f(:),2);
tol=1e-3*err;
while ((err>tol) && (K<Niter)),
    K=K+1;
    t=d-b;

    up=u;
    %Update u
    u=real(ifft2(FW.*(FF+fft2(t))));

    %Update d    
    tmp=u+b;
    d=sign(tmp).*max(Z,abs(tmp)-1/lambda);

    %Update b    
    b=tmp-d;
    
    err=sum(sum((up-u).^2));
end