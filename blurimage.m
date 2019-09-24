function B=blurimage(A,sz,c)
if (nargin==1)
    vector=false;
end
if (size(A,2)==1)
    A=reshape(A,sz);
    vector=true;
else
    vector=false;
end
B=zeros(size(A));
N=size(A,1);
w=[ 0:N/2 -(N/2-1):-1 ];
[w2x,w2y]=meshgrid(w);
for k=1:size(A,3)
    F=A(:,:,k);
    TF=fft2(F);
    K=zeros(N);
    for i=1:N
        for j=1:N
            wij=sqrt(w2x(i,j)^2+w2y(i,j)^2);
            % change denominator: larger = less blurring
            K(i,j)=exp(-wij^2*(128/N)^2/abs(c));
            %fourier transform of the Gaussian
        end
    end
    TG=TF.*K;
    G=ifft2(TG);
    B(:,:,k)=G;
end
if vector
    B=reshape(B,numel(B),1);
end

