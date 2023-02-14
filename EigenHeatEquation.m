clc;
clear;
set(gcf , 'Color' , 'w')

%% 
%Ŀ�� �� ���������׷����������·��̵ĵ���������׺�α�� ��
%  du/dt = d2u/dx2 
%  u(+-1 , t) = 0 
%  u(x , 0) = f(x) 
%  ����ʱ�䷽������һ����ʾŷ������ �� ʱ�䷽�������׷���


%% 
% ��������
N = 63;
% dt =  03 / N^(4);
dt = 1/N^4;

%% 
% ΢�־���
D_2 = convertmat(N + 1, 2, 1) * diffmat(N + 1, 2);

%% 
% ����ת������
S = convertmat(N + 1, 0, 1);

%% 
% ��ֵ����
B = [bc(N + 1, 'd', 'l'); bc(N+1 , 'd', 'r')];
BC = [1 ; 1];
P0 = eye(N + 1 , N + 1);
P0(N : N + 1 , 1 : N + 1) = B;

%% 
% ���
P = dt * P0 ^(-1) * S^(-1) * D_2;
% P = P0 ^(-1) * S^(-1) * D_2;

%% 
% �����������
lamda = eig(P);
x_P = real(lamda);
y_P = imag(lamda);
plot(x_P,y_P,'.')
hold on

%% 
% ���������α��
opts.npts = 500;
opts.ax = [-0.08 0.02 -0.05 0.05];
opts.levels = [-2.4:0.2:-1]; 
eigtool(P, opts)