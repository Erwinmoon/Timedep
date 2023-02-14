clc
clear
%% 
%目标 ： 利用新谱方法求解如下单行波方程组 ：
%  du/dt = du/dx 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  其中时间方向利用BD4 ， 
%  C = 3.45

% 改迭代和初值

% 参数
figure(44)
hold on
lw = 1;
set(gcf , 'Color' , 'w')
N = 79; 
M = N;
[~ , x] = cheb(N);
x = x;
dt = 3.5 / N^2;

D = convertmat(N + 1, 1, 0) * diffmat(N + 1, 1);
S = convertmat(N + 1, 0, 0);
P = eye(N + 1);
P(N + 1 , : ) = 1;

v1 = chebcoeffs(chebfun(@(x)exp(-200 * (x)^2) , N + 1));

tmax = 0.16; 
tplot = 10*dt;
plotgap = round(tplot/dt);
dt = tplot/plotgap;
nplots = round(tmax/tplot);
plotdata = [clenshaw(x , v1)'; zeros(nplots,N + 1)]; 
tdata = 0;
% clf, drawnow, h = waitbar(0,'please wait...');
T = 0:10*dt:10*nplots*dt;

P =  P ^ (-1) * S^(-1) * dt * D;

% 解方程
for i = 1:nplots, waitbar(i/nplots)
    for n = 1:plotgap
        v2   =  P * v1 + v1;
        v1 = v2;
    end
    plotdata(i+1,:) = clenshaw(x , v1)'; 
    tdata = [tdata; dt*i*plotgap];
end

% % 画图:
% clf, drawnow, waterfall(x,tdata,plotdata)
% axis([-1 1 0 tmax -2 2]), view(10,70), grid off
% colormap([0 0 0]), ylabel t, zlabel u, close(h)
% title('C= 3.41')

axis([-1,1,0,tmax,-1,1])
set(gca,'YTick',[0:tmax:tmax]) 
set(gca,'xTick',[-1:2:1]) 
set(gca,'zTick',[-1:1:1]) 
xlabel('x' , 'position' , [0 0 -1.2])
ylabel('t' , 'position' , [1.1 tmax/2 -1.2])
% zlabel('u' ,'rotation',0, 'position' , [-1.1 0 0.3])
pbaspect([2 2 0.5])
m = mesh(x,T,plotdata),view(10,40)
m.EdgeColor = [0 0 0]
