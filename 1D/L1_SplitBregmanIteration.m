function u=L1_SplitBregmanIteration(f,A,mu,lambda,Niter)

%=================================================
%
% L1 Split Bregman Iteration
% Version:
% -v1.2 - 06/23/2013
% -v1.0 - 03/31/2011
%
% u=L1_SplitBregmanIteration(f,A,mu,lambda,Niter)
%
% This function compute the solution
% u = arg min ||u||_1+0.5*mu||Au-f||_2^2
%
% by using the Split Bregman Iteration
%
% In this version we consider only
% the case where A is a square matrix
%
% f: measured data
% A: some linear operator in its matrix form
% mu: regularization coefficient
% lambda: "spliting" regularization coefficient
% Niter: maximum number of iteration
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
% 
% Note: typically mu=10, lambda=1,
%       Niter=10 work well
%=================================================
N=size(f,1);

d=zeros(N,1);
b=zeros(N,1);
u=zeros(N,1);

Z=zeros(N,1);
Ft=mu*A'*f;
IV=inv(mu*(A'*A)+lambda*eye(N));

err=norm(f,2);
tol=1e-3*err;
K=0;
while ((err>tol) && (K<Niter)),
    K=K+1;
    up=u;
    u=IV*(Ft+lambda*(d-b));
    tmp=u+b;
    d=sign(tmp).*max(Z,abs(tmp)-1/lambda);
    b=tmp-d;
   
    err=norm(u-up,2);
end
