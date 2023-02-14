clc
clear
%% 
%Ŀ�� �� �������׷���������µ��в������� ��
%  du/dt = du/dx 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  �ж����ַ���֮��Ĳ�� RK3


%%
% ��������
N = 300;
dt =  1 / N^(4);
T = 0.1;
t_end = 50000;
set(gcf , 'Color' , 'w')

%% 
% ��ֵ��ֵ����
a0 = chebcoeffs(chebfun('sin(2 * pi * x)' , N + 1));
g0 = a0;

%%
% ΢�־���
D = convertmat(N + 1, 2, 1) * diffmat(N + 1, 2);

%% 
% ����ת������
S = convertmat(N + 1, 0, 1);
Sa2 = S; %A2
%
%% 
% �߽���������
B = [bc(N + 1, 'd', 'l'); bc(N + 1, 'd', 'r')];
BC = [1 ; 1];
P0 = eye(N + 1 , N + 1);
P0(N : N + 1 , 1 : N + 1) = B;
S(N : N + 1 , 1 : N + 1) = B;%A1

%%
% ���
Pa1 =  S^(-1) * dt * D;
Pa2 =  P0^(-1) * Sa2^(-1) * dt * D;
for jj = 1 : 1 : t_end
    
    % 3��R-K����
    b = a0 + Pa1 * a0;
    c = (3 / 4) * a0 + (1 / 4) * (b + Pa1 * b);
    a1 = (1 / 3) * a0 + (2 / 3) * (c + Pa1 * c);
    a0 = a1;
    
    b = g0 + Pa2 * g0;
    c = (3 / 4) * g0 + (1 / 4) * (b + Pa2 * b);
    g1 = (1 / 3) * g0 + (2 / 3) * (c + Pa2 * c);
    g0 = g1;
    
    err(jj) = norm(a1-g1 , inf);
end

JJ = 1 : 1 : t_end;
plot(JJ , log10(err))



