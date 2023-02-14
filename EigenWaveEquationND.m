clc;
clear;
%% 
%目标 : 计算用新谱方法计算以下方程的迭代矩阵的谱和伪谱 ：
%  d2u/dt2 = d2u/dx2 
%  u(1 , t) = 0 
%  u_x(-1 , t) = 0 
%  u(x , 0) = f(x) 
%  其中时间方向利用中心差分方法 ， 时间方向用新谱方法

%% 
% 参数设置
N = 8;
dt = 5.1 / N^(2);

%% 
% 微分矩阵
D_2 = convertmat(N + 1, 2, 1) * diffmat(N + 1, 2);

%% 
% 基底转换矩阵
S = convertmat(N + 1, 0, 1);

%% 
% 边界条件
B = [bc(N , 'n', 'l'); bc(N + 1, 'd', 'r')];
BC = [0 ; 0];
P0 = eye(N + 1 , N + 1);
P0(N : N + 1 , 1 : N + 1) = B;

%%
% 迭代矩阵
P2 = dt^2 * P0 ^ (-1) * S^(-1) * D_2;

%% 
% 迭代矩阵的谱
lamda = eig(P2);
x_P = real(lamda);
y_P = imag(lamda);
plot(x_P,y_P,'.')
hold on

%% 
% 迭代矩阵的伪谱
opts.npts = 100;
opts.ax = [-2.5 1 -1.5 1.5];
opts.levels = [-10:-1]; 
eigtool(P2, opts)