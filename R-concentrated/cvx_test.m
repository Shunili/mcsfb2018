


n = 100;
A = randn(2*n,n);
b = randn(2*n,1);
tic
cvx_begin
   variable x(n)
   minimize( norm( A*x-b ) )
cvx_end

cvx_time=toc

tic
x_ls = A\b
mat_time=toc
norm(A*x_ls-b)
[x,x_ls]