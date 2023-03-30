%Graph for example

%size of the matrix
n=10; 
%vector for alpha
alpha_vector=(.01:.01:1)'; 
%salpha_vectord for random number generator
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
m=zeros(size(alpha_vector,1),4); 
for j=1:size(alpha_vector,1)
    %alpha
    alpha=alpha_vector(j); 
    %matrix B
    D=diag([0; alpha; alpha+r]); 
    %matrix \mathcal{B}
    A=D+v*v'; 
    %sqrt(\alpha)
    sa=sqrt(alpha); 
    %min eigenvalue of \mathcal{B}
    m(j,1)=min(eig(A)); 
    %new lower bound
    m(j,2)=.5*(1+alpha-sqrt((1+alpha).^2-4*alpha*v12)); 
    %Theorem 3.1
    m(j,3)=min(1,alpha).*(1-sqrt(1-v12)); 
    %Theorem 3.5
    m(j,4)=.5*(alpha+1-.5*(1+sa)*sqrt(1+2*sa*(1-2*v12)+alpha)-.5*abs(1-sa)*sqrt(1-2*sa*(1-2*v12)+alpha)); 
end
%m, m1, m2, m3
plot(alpha_vector,m(:,1),'-',alpha_vector,m(:,2),'-o',alpha_vector,m(:,3),'-x',alpha_vector,m(:,4),'-.'); xlabel('\alpha');
legend('min eigval','new lower bound','Thm 3.1','Thm 3.5','Location','SouthEast');
