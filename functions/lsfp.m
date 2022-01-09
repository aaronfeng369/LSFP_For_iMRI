%=====================================================================
% argmin f(x) + g(B(x))
% where 
% x = L + S, f(x) = ||A(L + S) - b||^2
% B = [I,0,framlet_L,0;0,\nabla_t,0,framlet_S;]^T
% g = (|.|_*,|.|_1,|.|_1,|.|_1,|.|_1,|.|_1)
% The method in the following uses the primal dual fixed point method: 
% http://math.sjtu.edu.cn/faculty/xqzhang/publications/CHZ_IP.pdf
%=====================================================================
function [L,S,f] = lsfp(param)
M=param.E'*param.d;
[nx,ny,nt]=size(M);
M=reshape(M,[nx*ny,nt]);
% M = zeros(size(M));
S = zeros(nx*ny,nt);
L = zeros(nx*ny,nt);
pt_L = zeros(size(L));
pt_S = zeros(size(S));
%initialization
frame = 1; Level = 1;
[D,R] = GenerateFrameletFilter(frame);
p_L = Framelet_D(L,nx,ny,nt,D,Level);
p_S = Framelet_D(S,nx,ny,nt,D,Level);
itr = 0;
gamma = 0.5;
lambda = 1.5e-1;  % to be modified!
c = lambda/gamma;
fprintf('\n ********** L+S reconstruction **********\n')
while(1),
	itr = itr + 1;
    M0  = M;
    gradient = param.E'*(param.E*reshape(L + S,[nx,ny,nt]) - param.d);
    gradient = reshape(gradient,[nx*ny,nt]);
    %===================================================================================
    Y_L = L - gamma*gradient - gamma*pt_L - gamma*Framelet_T(p_L,nt,R,Level);
    Y_S = S - gamma*gradient - gamma*Wtxs(pt_S) - gamma*Framelet_T(p_S,nt,R,Level);
    Par_L = c*Y_L + pt_L;
    Par_S = c*Wxs(Y_S) + pt_S;
    Par_pl = Cell_Add(Framelet_D(Y_L,nx,ny,nt,D,Level),p_L,c,1);  
    Par_ps = Cell_Add(Framelet_D(Y_S,nx,ny,nt,D,Level),p_S,c,1);      
    %===================================================================================

    %===============================================================
    [Ut,St,Vt] = svd(Par_L,0);
    pt_L = Ut*diag(Project_inf(diag(St),St(1)*param.lambda_L))*Vt';
    pt_S = Project_inf(Par_S,param.lambda_S);
    p_L = Project_inf_Cell(Par_pl,param.lambda_spatial_L);
    p_S = Project_inf_Cell(Par_ps,param.lambda_spatial_S);
    %===============================================================
    
    %================================================================================    
    L = L - gamma*gradient - gamma*pt_L - gamma*Framelet_T(p_L,nt,R,Level);
    S = S - gamma*gradient - gamma*Wtxs(pt_S) - gamma*Framelet_T(p_S,nt,R,Level);
    %================================================================================
    
    %=============================================
    f(itr) = Compute_loss(param,L,S,nx,ny,nt,D,Level);
    % print cost function and solution update
    fprintf(' ite: %d , cost: %f3\n', itr,f(itr)); 
    % stopping criteria 
    if (itr > param.nite)
        break;
    end
    %=============================================
end
[Ut,St,Vt]=svd(L,0);
L = Ut(:,1)*St(1)*Vt(:,1)';
L = reshape(L,nx,ny,nt);
S = reshape(S,nx,ny,nt);
end

function res = Wxs(x)
res = x(:,[2:end,end]) - x;
end

function res = Wtxs(x)
res = x(:,[1,1:end - 1]) - x;
res(:,1) = -x(:,1);
res(:,end) = x(:,end - 1);
end

function W = Framelet_D(x,nx,ny,nt,D,Level)
W = cell(1,nt);
for i = 1:nt
    A = reshape(x(:,i),nx,ny);
    W{i} = FraDecMultiLevel(A,D,Level);
end
end

function WT = Framelet_T(W,nt,R,Level)
for i = 1:nt
    temp = FraRecMultiLevel(W{i},R,Level);
    WT(:,i) = temp(:);
end
end


function s = Cell_Add(C1,C2,c1,c2)
s = C1;
[m n] = size(C1{1}{1});
for i = 1:length(s)
   for j = 1:length(s{1})
       for k = 1:m
           for r = 1:n
               s{i}{j}{k,r} = c1*C1{i}{j}{k,r} + c2*C2{i}{j}{k,r};
           end
       end
   end
end
end



function s = Project_inf(x,c)
s = (x/c ./ max(abs(x/c),1))*c;
end

function s = Project_inf_Cell(X,c)
s = X;
[m n] = size(X{1}{1});
for i = 1:length(s)
   for j = 1:length(s{1})
       for k = 1:m
           for r = 1:n
               temp = X{i}{j}{k,r};
               s{i}{j}{k,r} = (temp/c ./ max(abs(temp/c),1))*c;
           end
       end
   end
end
end

function s = Norm_Cell(X)
s = 0;
[m n] = size(X{1}{1});
for i = 1:length(X)
   for j = 1:length(X{1})
       for k = 1:m
           for r = 1:n
               s = s + norm(X{i}{j}{k,r}(:),1);
           end
       end
   end
end
end


function s = Compute_loss(param,L,S,nx,ny,nt,D,Level)
s = 0;
resk = param.E*reshape(L+S,[nx,ny,nt])-param.d;
[Ut,St,Vt] = svd(L,0);
C = diff(S,1,2);
s = 1/2*norm(resk(:),2)^2 + param.lambda_L*sum(diag(St))+param.lambda_S*norm(C(:),1);
% s = s + param.lambda_spatial*(Norm_Cell(Framelet_D(L,nx,ny,nt,D,Level)));
% s = s + param.lambda_spatial*(Norm_Cell(Framelet_D(S,nx,ny,nt,D,Level)));
end