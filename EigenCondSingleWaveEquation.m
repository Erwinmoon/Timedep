clc
clear
%% 
%目标 ： 计算用新谱方法计算以下方程的迭代矩阵的最大特征值条件数 ：
%  du/dt = du/dx 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  其中时间方向利用一阶显示欧拉方法 ， 时间方向用新谱方法


%%
% 参数设置
N = 63;
dt = 1 / N^2;
% dt = 0.0001;

%% 
% 数值初值条件a_0 = chebcoeffs(chebfun('1 - x' , N + 1));

%%
% 微分矩阵
D_0 = convertmat(N + 1, 1, 0) * diffmat(N + 1, 1);

%% 
% 基底转换矩阵
S_0 = convertmat(N + 1, 0, 0);
%
%% 
% 边界条件矩阵
P0 = eye(N + 1 , N + 1);
P0(N + 1 , 1 : N) = ones(1 , N);

%%
% 迭代矩阵
P =  P0 ^ (-1) * S_0^(-1) * dt * D_0;
Q = P';

%% 
% 左特征值
[xp, yp] = eig(P);
eigenvaluep = diag(yp);

[xq  , yq] = eig(Q);
eigenvalueq = diag(yq);

% max eigenvalue condition number
[ip, iq] = maxeigenvalue(eigenvaluep, eigenvalueq);

rvector = xp(: , ip);
lvector = xq(: , iq);

keigenmax = norm(lvector,2) * norm(rvector,2) / abs(lvector' * rvector);

% other eigenvalues condition number
keigen = zeros(5 , length(P));
for i = 1 : 1 : length(P)
    [ip, iq] = maxeigenvalue(eigenvaluep, eigenvalueq);
    if eigenvaluep(ip) == 0
        break;
    end
    rvector = xp(: , ip);
    lvector = xq(: , iq);
    keigen(1 , i) = real(ip);
    keigen(2 , i) = real(iq);
    keigen(3 , i) = eigenvaluep(ip);
    keigen(4 , i) = conj(eigenvalueq(iq));
    keigen(5 , i) = real(norm(xq(: , iq),2) * norm(xp(: , ip),2) / abs(xq(: , iq)' * xp(: , ip)));
    eigenvaluep(ip) = 0;
    eigenvalueq(iq) = 0;
end

%%
% 其他特征值的条件数
function [ip, iq] = maxeigenvalue(eigenvaluep, eigenvalueq)
    rlamda = max(abs(eigenvaluep));
    for ip = 1 : 1 : length(eigenvaluep)
       if rlamda ==  abs(eigenvaluep(ip))
           break;
       end
    end
    llamda = max(abs(eigenvalueq));
    for iq = 1 : 1 : length(eigenvalueq)
       if llamda ==  abs(eigenvalueq(iq))
           break;
       end
    end
    if norm(conj(eigenvalueq(iq)) - (eigenvaluep(ip)), 2) > 1e-12
        eigenvalueq(iq) = 0;
        llamda = max(abs(eigenvalueq));
        for iq = 1 : 1 : length(eigenvalueq)
            if llamda ==  abs(eigenvalueq(iq))
                  break;
            end
        end
    end
end