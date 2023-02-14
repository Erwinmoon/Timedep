clc;
clear;
set(gcf , 'Color' , 'w')
%% 
%Ŀ�� �� ���������׷����������·��̵ĵ���������׺�α�� ��
%  d2u/dt2 = d2u/dx2 
%  u(+-1 , t) = 0 
%  u(x , 0) = f(x) 
%  ����ʱ�䷽���������Ĳ�ַ��� �� ʱ�䷽�������׷���


%% 
% ��������
N = 500;
dt = 0.00000001;

%% 
% ΢�־���
D_2 = convertmat(N + 1, 2, 1) * diffmat(N + 1, 2);

%% 
% ����ת������
S = convertmat(N + 1, 0, 1);

%% 
% �߽���������
b = [1 ; -1];
B = zeros(2 , N + 1);
for n = 1 : 1 : N + 1
   B( : , n) = cos((n - 1) * acos(b)); 
end
P0 = eye(N + 1 , N + 1);
P0(N : N + 1 , 1 : N + 1) = B;

%%
% ��������
P = dt * P0 ^ (-1) * S^(-1) * D_2;

%% 
% �����������
lamda = eig(P);
x_P = real(lamda);
y_P = imag(lamda);
plot(x_P,y_P,'.')
hold on

%% 
% ���������α��
opts.npts = 100;
opts.ax = [-200 20 -110 110];
opts.levels = [-0.75:0.25 :0.25]; 
B = eigtool(P, opts)