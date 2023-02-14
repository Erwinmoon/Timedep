clc
clear
%% 
%Ŀ�� �� ���������׷����������·��̵ĵ���������׺�α�� ��
%  du/dt = du/dx 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  ����ʱ�䷽������һ����ʾŷ������ �� ʱ�䷽���þ��׷���


%%
% ��������
N = 63;
% dt = 0.01 / N;
dt = 1/N^4;
% ΢�־���
[D_0 , xx] = cheb(N);


%% 
% �߽���������
BC = zeros(2 , N + 1);
BC(1 ,1) = 0;
B(2 , :) = D_0(N + 1 , : );
D_0 = D_0^2;
D_0(1 , : ) = BC(1 , : );
D_0(N + 1 , : ) = BC(2 , :);

%% 
% �����������
% lamda = eig(D_0(2 : N , 2 : N) * dt);
lamda = eig(D_0 * dt);
% m = min(lamda)
x_P = real(lamda);
y_P = imag(lamda);
plot(x_P,y_P,'.')
hold on

%% 
% ���������α��
opts.npts = 500;
opts.ax = [-0.08 0.02 -0.05 0.05];
opts.levels = [-10:-1]; 
% eigtool(D_0(2 : N , 2 : N)  * dt, opts)
eigtool(D_0 * dt, opts)
% 
% %%
% %BD4�ȶ���
% z = exp(1i*pi*(0:200)/100);
% d = 1-1./z;
% plot([-40 40],[0 0]), hold on, plot([0 0],[-40 40])
% r = 0; 
% for i = 1:3, r = r+(d.^i)/i; plot(r), end % orders 1-6
% axis([-15 35 -25 25]), axis square, grid on
% title('backward differentiation')
