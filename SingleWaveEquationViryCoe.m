clc
clear
%% 
%Ŀ�� �� �������׷���������µ��в������� ��
%  du/dt = -c(x) * du/dx
%  c(x) = (sin(x - 1))^(2) + 1 / 5
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  ����ʱ�䷽������һ����ʾŷ������ �� ʱ�䷽�������׷���


%%
% ��������
BASE = 200;
y_num = zeros(1.01 * BASE + 1 , 2);
NN = [BASE , 1.01 * BASE];
for kk = 1 : 1 : 2
    
N = NN(1 , kk);
dt = 1 / N^(2);
T = 1;
t_end = floor(T / dt);
x = chebfun('x' , N + 1);

%% 
% ��ֵ��ֵ����
a_0 = chebcoeffs(chebfun('exp(-100 * (x )^2)' , N + 1));

%%
% ΢�־���
D_0 = convertmat(N + 1, 1, 0) * multmat(N + 1 ,1 / 5  + (sin(x))^(2) , 1) *...
diffmat(N + 1, 1);

%% 
% ����ת������
S_0 = convertmat(N + 1, 0, 0);
%
%% 
% �߽���������
P0 = eye(N + 1 , N + 1);
P0(N + 1 , 1 : N) = ones(1 , N);

%%
% ���
P =  P0 ^ (-1) * S_0^(-1) * dt * D_0;
for jj = 1 : 1 : t_end
    a_1 = P * a_0 + a_0;
    a_0 = a_1;
%     
    % 3��R-K����
%     b = a_0 + P * a_0;
%     c = (3 / 4) * a_0 + (1 / 4) * (b + P * b);
%     a_1 = (1 / 3) * a_0 + (2 / 3) * (c + P * c);
%     a_0 = a_1;
    
    % 4��R-K����
%     k1 = a_0;
%     k2 = a_0 + (1 / 2) * P * k1;
%     k3 = a_0 + (1 / 2) * P * k2;
%     a_1 = a_0 + P * k3;
%     a_0 = a_1;
end

%% 
% ͼ��
M = 1000;
xx = linspace(-1 , 1 , M);
A = zeros(M , N + 1);
A( : , 1) = ones(M , 1);
for n = 1 : 1 : N
   A( : , n + 1) = cos(n * acos(xx')); 
end

y_num(1 : length(a_1) , kk) = a_1;

end
%%
% ���
err = norm(y_num( 1 : BASE + 1 , 1) - y_num(1 : BASE + 1 , 2) , 2)