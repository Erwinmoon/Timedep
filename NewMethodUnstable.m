clc
clear
%% 
%Ŀ�� �� �������׷���������µ��в������� ��
%  du/dt = du/dx 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  ����ʱ�䷽������BD4 �� 
%  C = 3.4188 �� 3.4189

% �ĵ����ͳ�ֵ

% ����
N = 80; 
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
clf, drawnow, h = waitbar(0,'please wait...');

P =  P ^ (-1) * S^(-1) * dt * D;

% �ⷽ��
for i = 1:nplots, waitbar(i/nplots)
    for n = 1:plotgap
        v2   =  P * v1 + v1;
        v1 = v2;
    end
    plotdata(i+1,:) = clenshaw(x , v1)'; 
    tdata = [tdata; dt*i*plotgap];
end

% ��ͼ:
clf, drawnow, waterfall(x,tdata,plotdata)
axis([-1 1 0 tmax -2 2]), view(10,70), grid off
colormap([0 0 0]), ylabel t, zlabel u, close(h)
title('C= 3.41')