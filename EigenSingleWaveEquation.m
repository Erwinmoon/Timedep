clc
clear
clf

figure(1)
set(gcf , 'Color' , 'w')
lw = 1;
FS = 'FontSize'; fs = 15;
IN = 'interpret'; LX = 'latex';
%% 
%目标 ： 计算用新谱方法计算以下方程的迭代矩阵的谱和伪谱 ：
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
% 数值初值条件
a_0 = chebcoeffs(chebfun('1 - x' , N + 1));

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

%% 
% 迭代矩阵的谱
% lamda = eig(P);
% x_P = real(lamda);
% y_P = imag(lamda);
% % plot(x_P,y_P,'.k' , 'LineWidth',2)
% scatter(x_P,y_P,100,'.k')
% % axis([-4000 0 -2000 2000]),axis square
% hold on

%% 
% 迭代矩阵的伪谱
opts.npts = 500;
opts.ax = [-0.63 0.1 -0.32 0.32];
opts.levels = [-10:-1]; 
eigtool(P, opts)

%%
%BDF稳定域
% z = exp(1i*pi*(0:200)/100);
% d = 1-1./z;
% plot([-40 40],[0 0],'k'),
% hold on, 
% plot([0 0],[-40 40],'k')
% r = 0; 
% for i = 1:1, r = r+(d.^i)/i;end % orders 1-6
% for i = 2:6, r = r+(d.^i)/i; plot(r , 'k'), end
% axis([-4 4 -4 4]), axis square, grid on
% % text(1.5,1,'1')
% text(0.5,1.5,'2')
% text(0.0,1.6,'3')
% text(-0.6,1.7,'4')
% text(-1.5,1.5,'5')
% text(-1.7 , 0.5 , '6')
% text(0.3,0.2,'$\Delta t = 4.5 \times 10^{-4}$',FS,fs,IN,LX)