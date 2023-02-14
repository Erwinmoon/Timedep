clc;
clear;
%% 
%Ŀ�� �� ���Թ������u_t = u_x
%% 
% ��������
figure(44)
pbaspect([8 6 1])
hold on
lw = 1;
set(gcf , 'Color' , 'w')
N = 200;
dt = 3 / 500^2;
t_end = 10000;
[D , xx] = cheb(N);
D(1 , :) = zeros(1 , N + 1);
D(1 , 1) = 1;

%% ------------------RK3-------------------------------------------
% ��¼
err = zeros(t_end , 1);
J = 1 : 1 : t_end;

% ��ֵ��ֵ����
a_0 = exp(-200 .* (xx).^2);

% ����ת������
P = dt * D;

% ���
for j = 1 : 1 : t_end
    % 3��R-K����
    b = a_0 + P * a_0;
    c = (3 / 4) * a_0 + (1 / 4) * (b + P * b);
    a_1 = (1 / 3) * a_0 + (2 / 3) * (c + P * c);
    a_0 = a_1;
    a_1(1) = 0;
    
    y_num = a_1;
    y_exact = exp(-200 * (xx + dt * j) .^ 2);
    
    err(j , 1) = norm(y_num - y_exact , inf);
end
    plot(J , err  ,':k', 'LineWidth',lw)
    
%% ----------------------BDF3--------------------------------------
% ��¼
err = zeros(t_end - 2 , 1);
J = 3 : 1 : t_end;
% ��ֵ��ֵ����
a_0 = exp(-200 .* (xx).^2);
coe = zeros(N + 1 , 4);
coe( : , 1) = a_0;

% ���

% ��������
P = dt * D;
Q = (eye(N + 1) - (6 / 11) * P) ^ (-1);

% �ಽ����ֵ���ã���3����
for j = 1 : 1 : 2
    
    % 3��R-K����
    b = coe( : , j) + P * coe( : , j);
    c = (3 / 4) * coe( : , j) + (1 / 4) * (b + P * b);
    coe( : , j + 1) = (1 / 3) * coe( : , j) + (2 / 3) * (c + P * c);
    
end

% BD3��ⷽ��
for j = 1 : 1 : t_end - 2
    IV = Q * ((18 / 11) * coe( : , 3) - (9 / 11) * coe( : , 2) + (2 / 11) *...
        coe( : , 1));
    
    coe( : , 1 : 2) = coe( : , 2 : 3);
    coe( : , 3) = IV;
 
    y_exact = exp(-200 * (xx + dt * (j + 2)).^2);
    y_num = IV;
    err(j , 1) = norm(y_num - y_exact , inf);
end 
plot(J , err , '--k' , 'LineWidth',lw)
    
%% ---------------------AB4-----------------------------------
err = zeros(t_end - 3 , 1);
J = 4 : 1 : t_end;
% ��ֵ��ֵ����
a_0 = exp(-200 .* (xx).^2);
coe = zeros(N + 1 , 4);
coe( : , 1) = a_0;

% ��������
P = dt * D;

% �ಽ����ֵ���ã���4����
for j = 1 : 1 : 3
    
    % 3��R-K����
    b = coe( : , j) + P * coe( : , j);
    c = (3 / 4) * coe( : , j) + (1 / 4) * (b + P * b);
    coe( : , j + 1) = (1 / 3) * coe( : , j) + (2 / 3) * (c + P * c);
        
end

% A-B4��ⷽ��
for j = 1 : 1 : t_end - 3
    IV = coe( : , 4) + (1 / 24) * (55 * P * coe( : , 4) - 59 * P * coe( : , 3)  ...
        + 37 * P * coe( : , 2) - 9 * P * coe( : , 1));
    
    coe( : , 1 : 3) = coe( : , 2 : 4);
    coe( : , 4) = IV;   
    y_exact = exp(-200 * (xx + dt * (j + 3)).^2);
    y_num = IV;
    err(j , 1) = norm(y_num - y_exact , inf);
end    
plot(J , err  , 'k','LineWidth',lw+1)

legend('RK3','BDF3' ,'AB4', 'Location','east')
xlabel('number of steps')
ylabel('error')
axes('position',[0.2 0.6 0.36 0.3])
plot(J , err  , 'k','LineWidth',lw+1)