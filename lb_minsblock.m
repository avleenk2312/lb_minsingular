% Estimate min singular value of a non-singular 2 x 2 block matrix [A B; C D].
% The exact, as well as the three estimates (est1, est2, est3) of the min
% singular value, are calculated.

A=randn(2,2); B=randn(2,1); C=randn(1,2); D=randn; %random blocks
G=[A B; C D];       % global matrix
exact=min(svd(G));      % exact min singular value of G
cv=cc([A B]',[C D]');
c1=cc(A,B); 
c2=cc(C,D);
a=minsvd(A); b=minsvd(B); c=minsvd(C); d=minsvd(D); %minimum positive singular values of the blocks
ae=min(eig(A*A')); be=min(eig(B*B')); ce=min(eig(C*C')); de=min(eig(D*D')); %minimum eigenvalues for weyl's 
est1=cv*min(min(svd([A,B])),min(svd([C,D])));
est2=cv*min(max(c1*min(a,b),sqrt(ae+be)), max(c2*min(c,d), sqrt(ce+de));
est3=cv*min(c1,c2)*min([a,b,c,d]);
[exact, est1, est2, est3]

 
function c=cc(A,B)       % compute the function sqrt(c(A,B))
k=rank(A)+rank(B)-size(A,1);
[U1,S1,~]=svd(A,'econ'); 
[U2,S2,~]=svd(B,'econ');
e=1e-13; 
r1=sum(diag(S1)>e); 
r2=sum(diag(S2)>e);
s=svd(U1(:,1:r1)'*U2(:,1:r2));
r=sum(s<1-e & s>e);
c=1; 
if r>0 
   c=sqrt(1-s(k+1)); 
end
end
 
function m=minsvd(A)          % find the minimum positive singular value of A
m=Inf; 
if A==zeros(size(A)) 
   return; 
end
s=svd(A); 
m=min(s(s>1e-13));
end
