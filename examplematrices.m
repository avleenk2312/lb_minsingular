% Estimates for the minimum singular value of non-singular matrices A-K and U considered in Examples 3, 6 and 7 of the
% article: "New lower bounds for the minimum singular value of a matrix"
% Please "uncomment" the matrix you would like to check
% By default the matrix "U" from "Example 3" of the article is uncommented

% The following are calculated:
% 1. The exact minimum singular value of G.
% 2. The four "row-formulation" estimates (est1,est2,est3,est4)
% 3. The four "column-formulation" estimates (est1_t,est2_t,est3_t,est4_t)
% 4. The maximum of all the eight estimates is presente)


index=[1 1;1 2;2 1;2 2];
%%%%%%%%%% Please uncomment one of the following and comment the one you do NOT want to use %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% U-matrix from Example 3 %%%%%%%%%%%%%%%%%%%%%%
G=[10 0 0;4 2 0;1 1 6] %matrix U in article
%%%%%%%%%%%%%%%%% H-matrices from Example 6 %%%%%%%%%%%%%%%%%%%%
%G=[8 -2 -1;-5 7 -3;-3 -4 5] %matrix A in article
%G=[7 -3 -2;-2 5 -1; -3 -4 9] %matrix B in article
%G=[-5 2 -4;3 -6 -2;-1 -4 -8] %matric C in article
%%%%%%%%%%%%%%%%% Matrices from Example 7 %%%%%%%%%%%%%%%%%%%%%%
%G=[10 1 1;1 20 1;1 1 30] %matrix D in article
%G=[10 1 1;1 20 1;10 1 30] %matrix E in article
%G=[10 1 1;1 20 1;20 1 30] %matrix F in article
%G=[10 1 1;10 20 1;20 1 30] %matrix G in article
%G=[3 2 0;1 9 5;0 5 7] %matrix H in article
%G=[2 -1 0;2 1 0;-4 -4 5] %matrix I in article
%G=[5 0 0;-4 9 4;-1 7 9] %matrix J in article
%G=[4 0 0;-1 5 0;0 5 4] %matrix K in article %Lower triangular matrix



for i=1:4 %as there are 4 possible partitions
%partitioning so that "matrix" in block form is G=[A B;C D]
in1=index(i,1);
in2=index(i,2);
A=G(1:in1,1:in2);
B=G(1:in1,in2+1:end);
C=G(in1+1:end,1:in2);
D=G(in1+1:end,in2+1:end);

% exact min singular value of G
exact=min(svd(G));        

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
formula=find(e==m);
str=["(3.27a): first row formula";"(3.27b): second row formula";"(3.27c): third row formula";"(3.28): fourth row formula";...
      "(3.31a): first column formula";"(3.31b): second column formula";"third column formula";"fourth column formula"];
fprintf("Current partition number is %d, with size of the leading block as: %d x %d\n",i,in1,in2);
fprintf("The best estimate is: %f, for the exact minimum singular value: %f, given by the following formula(s):\n",m,exact)
fprintf("\t\t\t\t %s \n",str(formula))
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


