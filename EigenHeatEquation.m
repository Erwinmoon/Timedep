clc;
clear;
set(gcf , 'Color' , 'w')

%% 
%目标 ： 计算用新谱方法计算以下方程的迭代矩阵的谱和伪谱 ：
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
% P = P0 ^(-1) * S^(-1) * D_2;

%% 
% 迭代矩阵的谱
lamda = eig(P);
x_P = real(lamda);
y_P = imag(lamda);
plot(x_P,y_P,'.')
hold on

%% 
% 迭代矩阵的伪谱
opts.npts = 500;
opts.ax = [-0.08 0.02 -0.05 0.05];
opts.levels = [-2.4:0.2:-1]; 
eigtool(P, opts)