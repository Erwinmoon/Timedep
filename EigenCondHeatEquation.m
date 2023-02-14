clc
clear
%% 
%目标 ： 计算用新谱方法计算以下方程的迭代矩阵的最大特征值条件数 ：
%  du/dt = d2u/dx2 
%  u(+-1 , t) = 0 
%  u(x , 0) = f(x) 
%  其中时间方向利用一阶显示欧拉方法 ， 时间方向用新谱方法


%% 
% 参数设置
N = 63;
% dt =  03 / N^(4);
dt = 1/N^4;

%% 
% 微分矩阵
D_2 = convertmat(N + 1, 2, 1) * diffmat(N + 1, 2);

%% 
% 基底转换矩阵
S = convertmat(N + 1, 0, 1);

%% 
% 边值条件
B = [bc(N + 1, 'd', 'l'); bc(N+1 , 'd', 'r')];
BC = [1 ; 1];
P0 = eye(N + 1 , N + 1);
P0(N : N + 1 , 1 : N + 1) = B;

%% 
% 求解
P = dt * P0 ^(-1) * S^(-1) * D_2;
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
    keigen(1 , i) = ip;
    keigen(2 , i) = iq;
    keigen(3 , i) = eigenvaluep(ip);
    keigen(4 , i) = eigenvalueq(iq);
    keigen(5 , i) = norm(xq(: , iq),2) * norm(xp(: , ip),2) / abs(xq(: , iq)' * xp(: , ip));
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
    if abs(eigenvalueq(iq) - eigenvaluep(ip)) > 1e-12
        eigenvalueq(iq) = 0;
        llamda = max(abs(eigenvalueq));
        for iq = 1 : 1 : length(eigenvalueq)
            if llamda ==  abs(eigenvalueq(iq))
                  break;
            end
        end
    end
end