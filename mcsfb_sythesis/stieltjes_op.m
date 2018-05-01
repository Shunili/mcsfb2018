function [z, z_dif, num_iter] = stieltjes_op(G, b, xi,t, k, tol)
%function [z, z_dif, res, num_iter] = stieltjes_approx(A, b, f, t, k, fAb, tol)
%
% STIELTJES_APPROX: the main algorithm for approximating f(A)*b
%
% A: the matrix. Can also be a function that returns A times a vector
% b: the vector
% f: the function handle, e.g., @sqrt
% t: the knots
% k: maximum number of iterations
% fAb: f(A)*b, optional
% tol: tolerance, used to check for convergence
%
% z: the computed vector z_k that approximates f(A)*b
% z_dif: the differences between z_k and z_{k-1} for all iterations
% res: the differences between z_k and f(A)*b for all iteraions, if
%      fAb is provided
% num_iter: number of actual iterations
%

n= size(xi,1);
z_dif = zeros(k,1);
sigma = zeros(n,k+1);
eta   = zeros(n,k+1);
mu    = zeros(n,k+2);
mu1   = zeros(n,k+2);
h = reshape(diff(t)/2, n,1);
g = reshape((t(1:end-1)+t(2:end))/2, n,1);


beta = sqrt(n*pi);
mu(:,1) = 1/beta;
gamma = pi * sum(xi(:,1).*mu(:,1));
v = b/beta;
z = gamma * v;

for j = 1:k
  sigma(:,1) = h/2.*mu(:,2) + g.*mu(:,1);
  for p = 1:j
    sigma(:,p+1) = h/2.*(mu(:,p)+mu(:,p+2)) + g.*mu(:,p+1);
  end
  
  alpha = sum( sigma(:,1).*mu(:,1) + h/4.*mu(:,1).*mu(:,2) );
  for p = 1:j
    alpha = alpha + (sigma(:,p+1)'*mu(:,p+1))/2;
  end
  alpha = pi * alpha;

  if j ~= 1
    for p = 0:j
      eta(:,p+1) = sigma(:,p+1)-alpha*mu(:,p+1)-beta*mu0(:,p+1);
    end
  else
    for p = 0:j
      eta(:,p+1) = sigma(:,p+1)-alpha*mu(:,p+1);
    end
  end
  
  beta1 = sum(eta(:,1).^2) + sum((eta(:,2)+h/2.*mu(:,1)).^2)/2;
  for p = 2:j
    beta1 = beta1 + (eta(:,p+1)'*eta(:,p+1))/2;
  end
  beta1 = sqrt(pi*beta1);
  
  mu1(:,2) = (eta(:,2) + h/2.*mu(:,1)) / beta1;
  for p = [0 2:j]
    mu1(:,p+1) = eta(:,p+1) / beta1;
  end
  

    if j ~= 1
      v1 = (G.L*v - alpha*v - beta*v0) / beta1;
    else
      v1 = (G.L*v - alpha*v) / beta1;
    end

  gamma = sum( xi(:,1).*mu1(:,1) );
  for p = 1:min(j,3)
    gamma = gamma + sum(xi(:,p+1).*mu1(:,p+1))/2;
  end
  gamma = pi * gamma;
  
  z1 = z + gamma*v1;
  
  % check residual
  num_iter = j;
  z_dif(j) = norm(z-z1)/norm(z1);
  if z_dif(j) < tol
    z = z1;
    return;
  end

  % update
  beta = beta1;
  v0 = v;
  v = v1;
  mu0 = mu;
  mu = mu1;
  z = z1;
end
