
%Page 34 for LeVeque
%This is the code for Example 9 on Page 23
%It plots the exact minimum singular value versus the "best approximation" for minimum singular values
%Best means, the maximum lower bounds out of 3.25-(a),-(b),-(c), 3.26, 3.29-(a),-(b),-(c), 3.30


clear vars

nv=2.^(3:9);
nn=size(nv,2);
ex=zeros(nn,1);
be=zeros(nn,1);

for iter=1:nn

n=nv(iter);
h=1/n;

A=-2*n^2*eye(n);
e = ones(n,1);
B= n^2*full(spdiags([e e],0:1,n,n));
C=B;
D=A;
G=[A B; C D]; % global matrix
exact=min(svd(G));        % exact min singular value of G
ex(iter)=exact;

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
e=[est1,est2,est3,est4,est1_t,est2_t,est3_t,est4_t];
m=max(e);
be(iter)=m;
formula=find(e==m);
str=["(3.25a): first row formula";"(3.25b): second row formula";"(3.25c): third row formula";"(3.26): fourth row formula";...
      "(3.29a): first column formula";"(3.29b): second column formula";"(3.29c): third column formula";"(3.30) fourth column formula"];
fprintf("The best estimate is: %f, for the exact minimum singular value: %f, given by the following formula(s):\n",m,exact)
fprintf("\t\t\t\t %s \n",str(formula))
end

figure(1)
loglog(nv,ex,'ro-',nv,be,'bs-',nv,nv,'gx-')

 
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



