function x = Lower_QP_spg(A,AT,c,rho,l,x,eps)
% solve the lower bounded QP problem via SPG
%  min x'*Q*x/2 - c'*x
%  s.t. l <= x,
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
err = abs(g(:)'*(x(:)-l));
feas = -min([g(:);0])/max(norm(g(:)),1);
alphas = 1;
f_1(1) = f;
k = 1;

% Main algorithm
while (err/max(abs(f),1) >= eps || feas >= 1e-4)
  f_max = max(f_1);
  d = max(x-alphas*g, l) - x;
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
  err = abs(g(:)'*(x(:)-l));
  feas = -min([g(:);0])/max(norm(g(:)),1);
  iter_in = iter_in + 1;
  k = k + 1;  
end




