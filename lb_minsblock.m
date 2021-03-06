% Estimate min singular value of a non-singular 2 x 2 block matrix G=[A B; C D].
% The following are calculated:
% 1. The exact minimum singular value of G.
% 2. The four "row-formulation" estimates (est1,est2,est3,est4)
% 3. The four "column-formulation" estimates (est1_t,est2_t,est3_t,est4_t)
% 4. The maximum of all the eight estimates is presented

A=randn(2,2); B=randn(2,1); C=randn(1,2); D=randn; %random blocks
G=[A B; C D];             % global matrix
exact=min(svd(G));        % exact min singular value of G

%For Theorem 3.9
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

% Estimates given by Corollary 3.10
est1_t=cv_t*min(minsvd([A;C]),minsvd([B;D]));
est2_t=Psi([A;C],[B;D]);
est3_t=cv_t*min(max(c1_t*min(a,c),Psi(A',C')), max(c2_t*min(b,d), Psi(B',D')));
est4_t=cv_t*min(c1_t,c2_t)*min([a,b,c,d]);

fprintf("The best estimate is: %f for the exact minimum singular value: %f\n", max([est1,est2,est3,est4,est1_t,est2_t,est3_t,est4_t]),exact)

 
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
