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
N = 300;
dt = 1 / N^2;
% dt = 0.0001;
% dt = 1;

% ΢�־���
[D_0 , xx] = cheb(N);


%% 
% �߽���������
BC = zeros(1 , N + 1);
BC(1 ,1) = 0;
D_0(1 , : ) = BC;
% D_0 = D_0(2 : N + 1 , 2 : N + 1);

%% 
% �����������
lamda = eig(D_0*dt);
x_P = real(lamda);
y_P = imag(lamda);
plot(x_P,y_P,'.')
hold on

%% 
% ���������α��
opts.npts = 500;
opts.ax = [-0.63 0.1 -0.32 0.32];
opts.levels = [-1:-1]; 
eigtool(D_0 * dt, opts)

% %%
% %BD4�ȶ���
% z = exp(1i*pi*(0:200)/100);
% d = 1-1./z;
% plot([-40 40],[0 0]), hold on, plot([0 0],[-40 40])
% r = 0; 
% for i = 1:3, r = r+(d.^i)/i; plot(r), end % orders 1-6
% axis([-50 10 -30 30]), axis square, grid on
% title('backward differentiation')

