clc
clear
%% 
%目标 ： 计算用新谱方法计算以下方程的迭代矩阵的谱和伪谱 ：
%  du/dt = du/dx 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  其中时间方向利用一阶显示欧拉方法 ， 时间方向用旧谱方法


%%
% 参数设置
N = 300;
dt = 1 / N^2;
% dt = 0.0001;
% dt = 1;

% 微分矩阵
[D_0 , xx] = cheb(N);


%% 
% 边界条件矩阵
BC = zeros(1 , N + 1);
BC(1 ,1) = 0;
D_0(1 , : ) = BC;
% D_0 = D_0(2 : N + 1 , 2 : N + 1);

%% 
% 迭代矩阵的谱
lamda = eig(D_0*dt);
x_P = real(lamda);
y_P = imag(lamda);
plot(x_P,y_P,'.')
hold on

%% 
% 迭代矩阵的伪谱
opts.npts = 500;
opts.ax = [-0.63 0.1 -0.32 0.32];
opts.levels = [-1:-1]; 
eigtool(D_0 * dt, opts)

% %%
% %BD4稳定域
% z = exp(1i*pi*(0:200)/100);
% d = 1-1./z;
% plot([-40 40],[0 0]), hold on, plot([0 0],[-40 40])
% r = 0; 
% for i = 1:3, r = r+(d.^i)/i; plot(r), end % orders 1-6
% axis([-50 10 -30 30]), axis square, grid on
% title('backward differentiation')

