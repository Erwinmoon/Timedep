clc
clear
lw = 1;
fz = 20;

%% 
%目标 ： 利用新谱方法求解如下单行波方程组 ：
%  du/dt = du/dx 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  判断两种方法之间的差别  AB4


%%
% 参数设置
N = 300;
dt = 0.1 / N^(2);
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
    
    % AB4
    IVa1 = coea1( : , 4) + (1 / 24) * (55 * Pa1 * coea1( : , 4) - 59 * Pa1 * coea1( : , 3)  ...
        + 37 * Pa1 * coea1( : , 2) - 9 * Pa1 * coea1( : , 1));    
    coea1( : , 1 : 3) = coea1( : , 2 : 4);
    coea1( : , 4) = IVa1;   
    
    IVa2 = coea2( : , 4) + (1 / 24) * (55 * Pa2 * coea2( : , 4) - 59 * Pa2 * coea2( : , 3)  ...
        + 37 * Pa2 * coea2( : , 2) - 9 * Pa2 * coea2( : , 1));    
    coea2( : , 1 : 3) = coea2( : , 2 : 4);
    coea2( : , 4) = IVa2;  
    
    err1(jj) = norm(IVa1-IVa2 , inf);
end


JJ = 1 : 1 : t_end;
semilogy(JJ , err1 , '-k','LineWidth',lw)
hold on


%% 
%目标 ： 利用新谱方法求解如下单行波方程组 ：
%  du/dt = d2u/dx2 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  判断两种方法之间的差别 RK3


%%
% 参数设置
N = 300;
dt =  1 / N^(4);
t_end = 50000;

%% 
% 数值初值条件
a0 = chebcoeffs(chebfun('sin(2 * pi * x)' , N + 1));
g0 = a0;

%%
% 微分矩阵
D = convertmat(N + 1, 2, 1) * diffmat(N + 1, 2);

%% 
% 基底转换矩阵
S = convertmat(N + 1, 0, 1);
Sa2 = S; %A2
%
%% 
% 边界条件矩阵
B = [bc(N + 1, 'd', 'l'); bc(N + 1, 'd', 'r')];
BC = [1 ; 1];
P0 = eye(N + 1 , N + 1);
P0(N : N + 1 , 1 : N + 1) = B;
S(N : N + 1 , 1 : N + 1) = B;%A1

%%
% 求解
Pa1 =  S^(-1) * dt * D;
Pa2 =  P0^(-1) * Sa2^(-1) * dt * D;
for jj = 1 : 1 : t_end
    
    % 3阶R-K方法
    b = a0 + Pa1 * a0;
    c = (3 / 4) * a0 + (1 / 4) * (b + Pa1 * b);
    a1 = (1 / 3) * a0 + (2 / 3) * (c + Pa1 * c);
    a0 = a1;
    
    b = g0 + Pa2 * g0;
    c = (3 / 4) * g0 + (1 / 4) * (b + Pa2 * b);
    g1 = (1 / 3) * g0 + (2 / 3) * (c + Pa2 * c);
    g0 = g1;
    
    err2(jj) = norm(a1-g1 , inf);
end

JJ = 1 : 1 : t_end;
semilogy(JJ , err2, '--k', 'LineWidth',lw)

legend('transport','heat', 'Location','east')
xlabel('number of steps')
ylabel('error')
set(gcf , 'Color' , 'w')
axis([0 50000 1e-21 1e-14])
set(gca,'yTick',[1e-21 1e-20 1e-19 1e-18 1e-17 1e-16 1e-15 1e-14]) 
set(gca,'xTick',[0 10000 20000 30000 40000 50000]) 
pbaspect([16 9 1])

% set(gca,'ytick',[1 2 3 4]);
% set(gca,'yticklabel',{'-21','-20','-19','-18'});




