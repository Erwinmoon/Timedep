clc;
clear;
%% 
%Ŀ�� : ���������׷����������·��̵ĵ���������׺�α�� ��
%  d2u/dt2 = d2u/dx2 
%  u(1 , t) = 0 
%  u_x(-1 , t) = 0 
%  u(x , 0) = f(x) 
%  ����ʱ�䷽���������Ĳ�ַ��� �� ʱ�䷽�������׷���

%% 
% ��������
N = 8;
dt = 5.1 / N^(2);

%% 
% ΢�־���
D_2 = convertmat(N + 1, 2, 1) * diffmat(N + 1, 2);

%% 
% ����ת������
S = convertmat(N + 1, 0, 1);

%% 
% �߽�����
B = [bc(N , 'n', 'l'); bc(N + 1, 'd', 'r')];
BC = [0 ; 0];
P0 = eye(N + 1 , N + 1);
P0(N : N + 1 , 1 : N + 1) = B;

%%
% ��������
P2 = dt^2 * P0 ^ (-1) * S^(-1) * D_2;

%% 
% �����������
lamda = eig(P2);
x_P = real(lamda);
y_P = imag(lamda);
plot(x_P,y_P,'.')
hold on

%% 
% ���������α��
opts.npts = 100;
opts.ax = [-2.5 1 -1.5 1.5];
opts.levels = [-10:-1]; 
eigtool(P2, opts)