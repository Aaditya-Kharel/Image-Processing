function [sharpness]=solvePDE(x) %x is a vector

k=0.001;
%delta=x(2); 
%lambda=x(3); 
c=x(1);
delta=0;
lambda=1;

%Read image
A=double(imread('inset.jpg'));
A2=A(101:300, 101:300);
[p,q]=size(A2);
B2=blurimage(A2,size(A2),10);


%Convert matrix columns into a vector (u is solution at t=0)
Avector = reshape(A2',[],1);
Bvector = reshape(B2',[], 1);

dx=1;
N=p;

%Forward Difference matrix (one-dimensional)
e=ones(N,1);
Dplus=spdiags([-e,e], [0,1], N, N);
Dplus(N,1)=1;
%Backward Difference matrix (one-dimensional)
Dminus=-(Dplus)';
%Centered Difference matrix (one-dimensional)
Dcenter=(Dplus+Dminus)/2;

sparseI=speye(N,N);

% Create n^2xn^2 Forward differenciation wrt x,
Dxplus=kron(Dplus, sparseI);
% Create n^2xn^2 Backward differentiation wrt x,
Dxminus=kron(Dminus, sparseI);

% Forward differentiation wrt y,
Dyplus=kron(sparseI, Dplus);
% Backward differentiation wrt y
Dyminus=kron(sparseI, Dminus);

% Centered differentiation wrt x,
Dxcenter=kron(Dcenter, sparseI);
% Centered differentiation wrt y
Dycenter=kron(sparseI, Dcenter);

%Time Stepping
u0=Bvector; %sampletext.jpg

T=100;

%parameters
dt=0.01;
%k=0.001;   %smaller values of k
%delta=0;
%lambda=1;
p0=1.05;

%plot
%figure(1);
finalA=reshape(u0, N, N);
%plot(finalA(:,50));

%surf
figure(3);
surf(finalA);
view(0,90);
shading flat;
u=u0;

m_ind=0:p-1;
n_ind=0:q-1;
[m2,n2]=meshgrid(m_ind, n_ind);
lambd=-4*sin(m2*pi/p).^2-4*sin(n2*pi/p).^2;
lambda_reshape=reshape(lambd, numel(lambd), 1);
for i=1:100
    ux=Dxcenter*u;
    uy=Dycenter*u;
    uxy=ux.^2+uy.^2;
    uxy_updated=ux.^2+uy.^2;
    for j=1:length(uxy_updated)
    if (uxy_updated(j)==0)
        uxy_updated(j)=0.0001; %if the value is zero, small non-zero
    end
    end
    gu=(1./(1+k^2*uxy)+delta*(uxy_updated).^(p0-2));
    g0=mean(gu);
    G=spdiags(gu, 0, N^2, N^2);
    A=Dxplus*G*Dxminus+Dyplus*G*Dyminus;
    I=speye(N^2, N^2);
    var=I-dt*A;
    rhs=(u+lambda*(u0-blurimage(u,[N,N],c)));
    %rhs=(u+lambda*blurimage(u0-blurimage(u,[N,N]),[N,N]));
    %eigs(A+100*speye(size(A)),10,'SM')
    [u_new,~,~,~,RESVEC]=minres(var,rhs, 1e-6, 100, @(b)mfun(b,g0,dt,lambda_reshape));
    %RESVEC
    %norm(lambda*(u0-blurimage(u,[N,N]))) look into this
    %norm(rhs)
    %norm(u_new)
    %pause;
    u=u_new;
    finalA=reshape(u, N, N);
    %figure(2);
    %plot(finalA(:,50));
    %pause;
    %figure(4)
    %surf(finalA);
    %view(0,90);
    %shading flat; %get rid of grid lines
end

figure(4)
surf(finalA);
view(0,90);
shading flat; %get rid of grid lines
pause(0.1)

%reshape u into matrix for fft2
%take fft2
%reshape back to vector for 1-norm
%- (1-norm) is what gets returned
u_mat=reshape(u,N,N);
fft_umat=fft2(u_mat);
u_vect=reshape(fft_umat,numel(fft_umat),1);
sharpness=-norm(u_vect,1);
[c,k]
[-sharpness]






