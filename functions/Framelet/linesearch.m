% nonmonotone line search used in the SPG method
function [x_new,f_new,g_new] = linesearch(A,AT,c,rho,f_max,gamma,delta,x,d)

global nf;

alphas = 1;

while (1)  
  x_new = x + alphas*d;
  [f_new, Qx] = func(A,AT,c,rho,x_new);
  nf  = nf + 1;
  if (f_new <= f_max + gamma*alphas*delta) || (alphas <= 1e-8)
    break
  else
    alphas = alphas/2;
  end  
end
g_new = Qx - c;