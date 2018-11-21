function [errors,reconstructions]=compute_all_omp(signal,dict,param)
D=dict'; % put atoms in columns;
[N,M]=size(D);
f=signal;
tol=1e-10;

%--------------------------------
% Normalize the dictionary atoms 
%--------------------------------
norm_Of_Atoms = sqrt(sum(D.^2,1));
zeronorm=(norm_Of_Atoms==0);
norm_Of_Atoms(zeronorm)=1;
normalized_D = D ./ repmat(norm_Of_Atoms,N,1);

%--------------------------------
% Compute sparse coefficients
%--------------------------------

coeffs=zeros(M,M);
idx=zeros(M,1);
res=f;

for i=1:M % change 1 and M here to allow user to specify minimum and maximum number of coefficients
    proj = normalized_D' * res;
    [~,pos] = max(abs(proj));
    pos = pos(1);
    idx(i) = pos;
    a = normalized_D(:,idx(1:i))\f;
    res = f - normalized_D(:,idx(1:i)) * a;
    coeffs(idx(1:i),i)=a;
    if sum(res.^2) < tol
        coeffs(idx(1:i),(i+1):M)=repmat(a,1,M-i);
        break;
    end
end

%--------------------------------
% Renormalize the coefficients
%--------------------------------

coeffs = coeffs ./ repmat(norm_Of_Atoms',1,M);

%--------------------------------
% Compute reconstruction errors
%--------------------------------

reconstructions=D*coeffs;
diffs=repmat(f,1,M)-reconstructions;
errors=sqrt(sum(diffs.*diffs,1));
errors=errors/norm(f,2);