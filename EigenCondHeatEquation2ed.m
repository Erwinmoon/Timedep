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
eigenvalue = diag(yp);

rlamda = max(abs(eigenvalue));
for i = 1 : 1 : length(P)
   if rlamda ==  abs(eigenvalue(i))
       break;
   end
end
eigenvalue(i) = 0;

rlamda = max(abs(eigenvalue));
for i = 1 : 1 : length(P)
   if rlamda ==  abs(eigenvalue(i))
       break;
   end
end

rvector = xp(: , i);

% 右特征值
[xq  , yq] = eig(Q);
eigenvalue = diag(yq);
llamda = max(abs(eigenvalue));

for j = 1 : 1 : length(Q)
   if llamda ==  abs(eigenvalue(j))
       break;
   end
end
eigenvalue(j) = 0;

llamda = max(abs(eigenvalue));
for j = 1 : 1 : length(P)
   if llamda ==  abs(eigenvalue(j))
       break;
   end
end
lvector = xq(: , j);

keigen = norm(lvector,2) * norm(rvector,2) / abs(lvector' * rvector)

