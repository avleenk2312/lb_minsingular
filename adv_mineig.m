%Graph for example

%size of the matrix
n=10; 
%vector for alpha
ee=(.01:.01:1)'; 
%seed for random number generator
rng(9); 
%random vector for B
r=rand(n-2,1); 
%random vector v
v=randn(n,1); 
%normalized v
v=v/norm(v); 
%v_1
v1=v(1); 
%v_1^2
v12=v1^2;
%collection of four lower bounds
m=zeros(size(ee,1),4); 
for j=1:size(ee,1)
    %alpha
    e=ee(j); 
    %matrix B
    D=diag([0; e; e+r]); 
    %matrix \mathcal{B}
    A=D+v*v'; 
    %sqrt(\alpha)
    es=sqrt(e); 
    %min eigenvalue of \mathcal{B}
    m(j,1)=min(eig(A)); 
    %new lower bound
    m(j,2)=.5*(1+e-sqrt((1+e).^2-4*e*v12)); 
    %Theorem 3.1
    m(j,3)=min(1,e).*(1-sqrt(1-v12)); 
    %Theorem 3.5
    m(j,4)=.5*(e+1-.5*(1+es)*sqrt(1+2*es*(1-2*v12)+e)-.5*abs(1-es)*sqrt(1-2*es*(1-2*v12)+e)); 
end
%m, m1, m2, m3
plot(ee,m(:,1),'-',ee,m(:,2),'-o',ee,m(:,3),'-x',ee,m(:,4),'-.'); xlabel('\alpha');
legend('min eigval','new lower bound','Thm 3.1','Thm 3.5','Location','SouthEast');
