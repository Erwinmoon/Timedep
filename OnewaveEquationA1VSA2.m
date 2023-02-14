clc
clear
%% 
%目标 ： 利用新谱方法求解如下单行波方程组 ：
%  du/dt = du/dx 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  判断两种方法之间的差别  AB4


%%
% 参数设置
N = 300;
dt = 1 / N^(4);
T = 1;
t_end = 50000;
set(gcf , 'Color' , 'w')

%% 
% 数值初值条件
a0 = chebcoeffs(chebfun('exp(-200 * (x)^2)' , N + 1));
coea1 = zeros(N + 1 , 4);
coea1( : , 1) = a0;
coea2 = coea1;

%%
% 微分矩阵
D = convertmat(N + 1, 1, 0) * diffmat(N + 1, 1);

%% 
% 基底转换矩阵
S = convertmat(N + 1, 0, 0);
Sa2 = S; %A2
%
%% 
% 边界条件矩阵
P0 = eye(N + 1 , N + 1);
P0(N + 1 , : ) = 1; 
S(N + 1 , 1 : N+1) = ones(1 , N+1);%A1

%%
% 求解

Pa1 =  S^(-1) * dt * D;
Qa1 = (eye(N + 1) - (12 / 25) * Pa1) ^ (-1);
Pa2 = P0^(-1) * Sa2^(-1) * dt * D;
Qa2 = (eye(N + 1) - (12 / 25) * Pa2) ^ (-1);

for j = 1 : 1 : 3
    
    % 3阶R-K方法
    b = coea1( : , j) + Pa1 * coea1( : , j);
    c = (3 / 4) * coea1( : , j) + (1 / 4) * (b + Pa1 * b);
    coea1( : , j + 1) = (1 / 3) * coea1( : , j) + (2 / 3) * (c + Pa1 * c);
        
end
coea2 = coea1;

for jj = 1 : 1 : t_end
    
    % BDF4
    IVa1 = coea1( : , 4) + (1 / 24) * (55 * Pa1 * coea1( : , 4) - 59 * Pa1 * coea1( : , 3)  ...
        + 37 * Pa1 * coea1( : , 2) - 9 * Pa1 * coea1( : , 1));    
    coea1( : , 1 : 3) = coea1( : , 2 : 4);
    coea1( : , 4) = IVa1;   
    
    IVa2 = coea2( : , 4) + (1 / 24) * (55 * Pa2 * coea2( : , 4) - 59 * Pa2 * coea2( : , 3)  ...
        + 37 * Pa2 * coea2( : , 2) - 9 * Pa2 * coea2( : , 1));    
    coea2( : , 1 : 3) = coea2( : , 2 : 4);
    coea2( : , 4) = IVa2;  
    
    err(jj) = norm(IVa1-IVa2 , inf);
end


JJ = 1 : 1 : t_end;
plot(JJ , log10(err))