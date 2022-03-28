% Estimate min singular value of a non-singular matrix A and \tilde{A}
% considered in the "Example 5" of the article
% "New lower bounds for the minimum singular value of a matrix"


% The following are calculated:
%   1. Loglog plot of best result of Theorem 3.9 on A and \tilde{A}
%   2. Loglog plot of best result of Corollary 3.10 on A and \tilde{A}

e=zeros(30,1);
m1=zeros(30,1);
m2=zeros(30,1);
m3=zeros(30,1);
m4=zeros(30,1);
for t=1:30

% the following are blocks for matrix "A" in Example 5 in block form is G=[A B;C D]

A=[t 10 0;3 2 -2;2 0 6];
A11=A(1:1,1:1);
A12=A(1:1,2:3);
A21=A(2:3,1:1);
A22=A(2:3,2:3);
%blocks for tilde(A)
At11=A(1:2,1:1);
At12=A(1:2,2:3);
At21=A(3:3,1:1);
At22=A(3:3,2:3);

% exact min singular value of G
e(t)=min(svd(A));

%best estimates of Theorem 3.9 and Corollary 3.10 on A
[m1(t),m2(t)]=lb(A11,A12,A21,A22);
%best estimates of Theorem 3.9 and Corollary 3.10 on \tilde{A}
[m3(t),m4(t)]=lb(At11,At12,At21,At22);
end
t=1:30;
figure(1)
loglog(t,e,'ro-',t,m3,'g^-',t,m1,'bs-')
title('Loglog Plot of Theorem 3.9 on A and \tilde{A}')
legend({'\sigma_{\min}(A)','Theorem 3.9 on \tilde{A}','Theorem 3.9 on A'},'Location','northeast')
figure(2)
loglog(t,e,'ro-',t,m4,'g^-',t,m2,'bs-')
title('Loglog Plot of Corollary 3.10 on A and \tilde{A}')
legend({'\sigma_{\min}(A)','Corollary 3.10 on \tilde{A}','Corollary 3.10 on A'},'Location','northeast')

function [m,n]=lb(A,B,C,D)
cv=cc([A B]',[C D]');     % sqrt(1-cosΘ)
c1=cc(A,B);               
c2=cc(C,D);
a=minsvd(A); b=minsvd(B); c=minsvd(C); d=minsvd(D);                         %minimum positive singular values of the blocks

% Estimates given by Theorem 3.9
est1=cv*min(min(svd([A,B])),min(svd([C,D])));
est2=Psi([A,B]',[C,D]');
est3=cv*min(max(c1*min(a,b),Psi(A,B)), max(c2*min(c,d), Psi(C,D)));
est4=cv*min(c1,c2)*min([a,b,c,d]);

% For Corollary 3.10
cv_t=cc([A ;C],[B; D]);
c1_t=cc(A',C'); 
c2_t=cc(B',D');

% Estimates given by Corollary 3.9
est1_t=cv_t*min(minsvd([A;C]),minsvd([B;D]));
est2_t=Psi([A;C],[B;D]);
est3_t=cv_t*min(max(c1_t*min(a,c),Psi(A',C')), max(c2_t*min(b,d), Psi(B',D')));
est4_t=cv_t*min(c1_t,c2_t)*min([a,b,c,d]);
m=max([est1,est2,est3,est4]);
n=max([est1_t,est2_t,est3_t,est4_t]);
end
 
function c=cc(A,B)  % compute the function sqrt(c(A,B))
n=size(A,1);
k=rank(A)+rank(B)-n;
[U1,~,~]=svd(A,'econ');      % orthogonal basis for A
[U2,~,~]=svd(B,'econ');      % orthogonal basis for B
e=1e-8; 
r1=rank(A); 
r2=rank(B);
s=svd(U1(:,1:r1)'*U2(:,1:r2));     % singular values of M=U1'*U2
r=sum(s<1-e & s>e);                % number of angles in (0,π/2)
c=sqrt(1+double(r1+r2==2*n));      % c=2 if both matrices are non-singular
if r>0 
   c=sqrt(1-s(k+1)); 
end
end
 

function m=minsvd(A)          % find the minimum positive singular value of A
m=Inf; 
if A==zeros(size(A))          %minimimum singular value or eigenvalue of O matrix is defined to be infinity
   return; 
end
s=svd(A); 
m=min(s(s>1e-8));
end


function d=Psi(A,B)           % calculates the value of Psi(A,B)=sqrt(psi(AA^T,BB^T))
r1=rank(A);        
r2=rank(B);
n=size(A,1);             
k=r1+r2-n;
[U1,~,~]=svd(A,'econ'); 
[U2,~,~]=svd(B,'econ'); 
U1=U1(:,1:r1);
U2=U2(:,1:r2);
s=svd(U1'*U2);
eps=1e-8; 
r=sum(s<1-eps & s > eps); % size(s)
am=minsvd(A);
bm=minsvd(B);
if r>0 
    s2=(1-(s(k+1))^2);                                              % 1-cos(Θ_k+1)^2 = sin(Θ_k+1)^2 
    d=sqrt(0.5*(am^2+bm^2-(0.5*(am+bm)*sqrt((am+bm)^2-4*am*bm*s2))...
        -(0.5*abs(am-bm)*sqrt((am-bm)^2+4*am*bm*s2))));
    return
end
if and(r1==n,r2<n)       % A is non-singular
        d=am;
elseif and(r2==n,r1<n)   % B is non-singular
        d=bm;
elseif r1+r2==2*n        % both A and B are non-singular  
        d=sqrt(am^2+bm^2);
else                     % both A and B are singular and r=0, i.e., all angles are either 0 or π/2
        d=min(am,bm);
end 
end


