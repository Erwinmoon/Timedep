clc;
clear;
%% 
%目标 ： 测试滚动误差u_t = u_xx
%% 
% 参数设置
figure(44)
pbaspect([8 6 1])
hold on
lw = 1;
set(gcf , 'Color' , 'w')
N = 30;
dt = 0.00000001;
t_end = 100000;
%%
% 微分矩阵
[D , xx] = cheb(N);
BC = zeros(2 , N + 1);
BC(1 , 1) = 1;
BC(2 , N + 1) = 1;
D = D^2;
D(1 , :) = BC(1 , :);
D(N + 1 , :) = BC(2 , :);

%% ------------------RK3-------------------------------------------
% 记录
err = zeros(t_end , 1);
J = 1 : 1 : t_end;

% 数值初值条件
a_0 = sin(2 * pi * xx);


P = dt * D;
P = sparse(P);

% 求解
for j = 1 : 1 : t_end
    % 3阶R-K方法
    b = a_0 + P * a_0;
    c = (3 / 4) * a_0 + (1 / 4) * (b + P * b);
    a_1 = (1 / 3) * a_0 + (2 / 3) * (c + P * c);
    a_0 = a_1;
    
    y_num = a_1;
    y_exact = exp(- 4 * pi * pi * dt * j) * sin( 2 * pi * xx);
    
    err(j , 1) = norm(y_num - y_exact , inf);
end
    plot(J , err ,':k' ,'LineWidth',lw)

%% ----------------------BDF3--------------------------------------
% 记录
err = zeros(t_end - 3 , 1);
J = 4 : 1 : t_end;
% 数值初值条件
a_0 = sin(2 * pi * xx);
coe = zeros(N + 1 , 4);
coe( : , 1) = a_0;

% 求解

% 迭代矩阵
P = dt * D;
Q = (eye(N + 1) - (12 / 25) * P) ^ (-1);

% 多步法初值设置（共3步）
for j = 1 : 1 : 3
    
    % 3阶R-K方法
    b = coe( : , j) + P * coe( : , j);
    c = (3 / 4) * coe( : , j) + (1 / 4) * (b + P * b);
    coe( : , j + 1) = (1 / 3) * coe( : , j) + (2 / 3) * (c + P * c);
    
end

% BD3求解方程
for j = 1 : 1 : t_end - 3
    IV = Q * ((48 / 25) * coe( : , 4) - (36 / 25) * coe( : , 3) + ...
        (16 / 25) * coe( : , 2) - (3 / 25) * coe( : , 1));
    
    coe( : , 1 : 3) = coe( : , 2 : 4);
    coe( : , 4) = IV;
 
    y_exact = exp(- 4 * pi * pi * dt * (j+3)) * sin( 2 * pi * xx);
    y_num = IV;
    err(j , 1) = norm(y_num - y_exact , inf);
end    
plot(J , err , '--k','LineWidth',lw)
    
%% ---------------------AB4-----------------------------------
err = zeros(t_end - 3 , 1);
J = 4 : 1 : t_end;
% 数值初值条件
a_0 = sin(2 * pi * xx);
coe = zeros(N + 1 , 4);
coe( : , 1) = a_0;

% 迭代矩阵
P =  dt * D;

% 多步法初值设置（共4步）
for j = 1 : 1 : 3
    
    % 3阶R-K方法
    b = coe( : , j) + P * coe( : , j);
    c = (3 / 4) * coe( : , j) + (1 / 4) * (b + P * b);
    coe( : , j + 1) = (1 / 3) * coe( : , j) + (2 / 3) * (c + P * c);
        
end

% A-B4求解方程
for j = 1 : 1 : t_end - 3
    IV = coe( : , 4) + (1 / 24) * (55 * P * coe( : , 4) - 59 * P * coe( : , 3)  ...
        + 37 * P * coe( : , 2) - 9 * P * coe( : , 1));
    
    coe( : , 1 : 3) = coe( : , 2 : 4);
    coe( : , 4) = IV;   
    y_exact = exp(- 4 * pi * pi * dt * (j+3)) * sin( 2 * pi * xx);
    y_num = IV;
    err(j , 1) = norm(y_num - y_exact , inf);
end    
plot(J , err , 'k','LineWidth',lw+1)


legend('RK3','BDF3','AB4', 'Location','east')
xlabel('number of steps')
ylabel('error')
set(gca,'yTick',[0:10^(-11):8*10^(-11)]) 
axes('position',[0.2 0.6 0.36 0.3])
plot(J , err  , 'k','LineWidth',lw+1)