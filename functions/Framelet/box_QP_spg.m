function x = box_QP_spg(A,AT,c,rho,l,u,x,eps)
% solve the box QP problem via SPG
%  min x'*Q*x/2 - c'*x
%  s.t. l <= x <= u,
% where Q = A'*A + rho I.

iter_in = 1;
% parameters for the SPG method %
gamma = 1e-4;
M = 20;
alphamin = 1e-8;
alphamax = 1;

% Initialization
f_1 = -inf*ones(M,1);
[n,m] = size(x);
[f, Qx] = func(A,AT,c,rho,x);
g = Qx - c;
lambda_1 = zeros(n,m);
I = find(g>0);
lambda_1(I) = g(I);
lambda_2 = lambda_1 - g;
err = - l*sum(lambda_1(:)) + g(:)'*x(:) + u*sum(lambda_2(:));
alphas = 1;
f_1(1) = f;
k = 1;

% Main algorithm
while (err/max(abs(f),1)) >= eps
  f_max = max(f_1);
  d = min(max(x-alphas*g, l),u) - x;
  delta = d(:)'*g(:);
  [x_new,f_new,g_new] = linesearch(A,AT,c,rho,f_max,gamma,delta,x,d);  
  f_1(mod(k,M)+1) = f_new;
  dx = x_new - x;
  dg = g_new - g;
  xdotg = dx(:)'*dg(:);
  xsqr =  dx(:)'*dx(:);
  alphas = max(alphamin,min(alphamax,xsqr/xdotg));
  x = x_new;
  g = g_new;
  f = f_new;  
  lambda_1 = zeros(n,m);
  I = find(g>0);
  lambda_1(I) = g(I);
  lambda_2 = lambda_1 - g; 
  err = - l*sum(lambda_1(:)) + g(:)'*x(:) + u*sum(lambda_2(:));
  iter_in = iter_in + 1;
  k = k + 1;
end



