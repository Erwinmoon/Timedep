clc;
clear;
set(gcf , 'Color' , 'w')
%% 
%目标 ： 计算用新谱方法计算以下方程的迭代矩阵的谱和伪谱 ：
%  d2u/dt2 = d2u/dx2 
%  u(+-1 , t) = 0 
%  u(x , 0) = f(x) 
%  其中时间方向利用中心差分方法 ， 时间方向用新谱方法


%% 
% 参数设置
N = 500;
dt = 0.00000001;

%% 
% 微分矩阵
D_2 = convertmat(N + 1, 2, 1) * diffmat(N + 1, 2);

%% 
% 基底转换矩阵
S = convertmat(N + 1, 0, 1);

%% 
% 边界条件矩阵
b = [1 ; -1];
B = zeros(2 , N + 1);
for n = 1 : 1 : N + 1
   B( : , n) = cos((n - 1) * acos(b)); 
end
P0 = eye(N + 1 , N + 1);
P0(N : N + 1 , 1 : N + 1) = B;

%%
% 迭代矩阵
P = dt * P0 ^ (-1) * S^(-1) * D_2;

%% 
% 迭代矩阵的谱
lamda = eig(P);
x_P = real(lamda);
y_P = imag(lamda);
plot(x_P,y_P,'.')
hold on

%% 
% 迭代矩阵的伪谱
opts.npts = 100;
opts.ax = [-200 20 -110 110];
opts.levels = [-0.75:0.25 :0.25]; 
B = eigtool(P, opts)