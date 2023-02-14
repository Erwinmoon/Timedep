clc
clear
clf

figure(1)
set(gcf , 'Color' , 'w')
lw = 1;
FS = 'FontSize'; fs = 16;
IN = 'interpret'; LX = 'latex';
%% 
%目标 ： 计算用新谱方法计算以下方程的迭代矩阵的谱和伪谱 ：
%  du/dt = du/dx 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  其中时间方向利用一阶显示欧拉方法 ， 时间方向用新谱方法


%%
% 参数设置
j = 1;
for N = 10 : 10 : 1000
    dt = 1;
    % dt = 0.0001;

    %%
    % 微分矩阵
    D = diffmat(N + 1,3);

    %% 
    % 基底转换矩阵
    S = convertmat(N + 1, 0, 2);
    %
    %% 
    % 边界条件矩阵
    P = eye(N + 1 , N + 1);
    B = [bc(N + 1, 'd', 'l'); bc(N+1 , 'd', 'r');bc(N, 'n', 'l')];
    P(N - 1 : N + 1 , :) = B;
    %%
    % 迭代矩阵
    Q =  P ^ (-1) * S^(-1) * dt * D;

    %% 
    % 迭代矩阵的谱
    lamda = eig(Q);
    
    %%
    % 画图参数
    C(j) = max(abs(lamda))/N^6;
    NN(j) = (N);
    j = j + 1;
end

plot(NN , C , 'k' , 'LineWidth',1)
xlabel('n')
ylabel('$\frac{\mathbf{\rho}}{\mathrm{n}^6}$', IN , LX ,'rotation',0, FS ,fs ,...
    'position' , [-120 0.075])
axis([10,1000,0.055,0.095])
set(gca,'yTick',[0.055:0.01: 0.095])
pbaspect([2 1 1])
