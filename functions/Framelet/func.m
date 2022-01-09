function [f,Qx] = func(A,AT,c,rho,x)
Qx = AT(A(x)) + rho *x;
f = sum(sum(x.*(0.5*Qx - c)));