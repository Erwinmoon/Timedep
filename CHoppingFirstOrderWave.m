clc
clear
%% 
%目标 ： 利用新谱方法求解如下单行波方程组 ：
%  du/dt = -c(x)*du/dx , c(x) = 1/5 + (sin(x))^2 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  其中时间方向利用一阶显示欧拉方法 ， 时间方向用新谱方法


%%
% 参数设置
N = 400;
M = 2000;
xx = linspace(-1 , 1 , 200)';
tdata = 0;
dt = 1 / 500^(2);
T = 1;
t_end = floor(T / dt);
% tol = eps;
tol = 1e-10;
plateau = zeros(t_end , 1);
plateau1 = zeros(t_end , 1);

%% 
% 数值初值条件
a_0 = chebcoeffs(chebfun('exp(-400 * (x-0.75)^2)' , N));
x = chebfun('x' , M + 1);
plotdata(1 , : ) = clenshaw(xx , a_0)';
%%
% 微分矩阵
D_0 = 3 * convertmat(M + 1, 1, 0) * multmat(M + 1 ,1 / 5  + (sin((x-1)^2))^(2) , 1) *...
diffmat(M + 1, 1);


%% 
% 基底转换矩阵
S_0 = convertmat(M + 1, 0, 0);
%
%% 
% 边界条件矩阵
% P0 = eye(M + 1 , M + 1);
% P0(M + 1 , 1 : M) = ones(1 , M);

%%
% 求解
P =  S_0^(-1) * dt * D_0;
a_1 = zeros(N, 1);
jj = 1;
while jj < t_end +1
    a_1 = P(1:N , 1:N) * a_0(1:N) + a_0(1:N); 
    a_1(N) = -sum(a_1(1:N-1));
    
    if mod(jj , 10000) == 0
        plotdata(jj/100+1,:) = clenshaw(xx , a_1)'; 
        tdata(jj/100+1) = dt*jj;
    end
    
    cutoff = CHopping(a_1, tol);
    N = cutoff(2);
    plateau(jj) = N;
    plateau1(jj) = cutoff(1);
    if N < length(a_0) + 1
        a_0 = a_1(1:N);
        jj = jj + 1;
    else
        a_0 = [a_0;zeros(N - length(a_0) , 1)];
    end     
    
end

%% 长度
% 
JJ = zeros(2*t_end , 1);
plotplateau = zeros(2 * t_end , 1);
plotplateau1 = zeros(2 * t_end , 1);
for ii = 1 : 1 : t_end
    JJ(2 * ii) = ii * dt;
    JJ(2 * ii + 1) = ii * dt;
    plotplateau(2*ii) = plateau(ii);
    plotplateau(2*ii-1) = plateau(ii);
    plotplateau1(2*ii) = plateau1(ii);
    plotplateau1(2*ii-1) = plateau1(ii);
end
figure(1)
set(gcf , 'Color' , 'w')
plot(JJ(1 : 2*t_end) , plotplateau,'-k' ,'LineWidth',1)
hold on
plot(JJ(1 : 2*t_end) , plotplateau1,':k' ,'LineWidth',1.5)
legend('with plateau','plateau dropped', 'Location','northwest')
axis([0,1,80,240])
xlabel('t' )
ylabel('n' , 'rotation',0 , 'position' , [-0.1 160] )

% axis([0 t_end 205 280])
figure(2)
set(gcf , 'Color' , 'w')
h = waterfall(xx,tdata,plotdata)
h.LineWidth = 1;
axis([-1,1,0,T,-0.5,1])
set(gca,'YTick',[0:0.2:T]) 
set(gca,'xTick',[-1:1:1]) 
set(gca,'zTick',[-0.5 0 1]) 
xlabel('x' , 'position' , [-0.05 0 -1.2])
ylabel('t' , 'position' , [1.2 T/2 + 0.2 -1.4])
view(10,70), grid off
colormap([0 0 0])
%% 
% 图像
% M = 1000;
% xx = linspace(-1 , 1 , M);
% A = zeros(M , N);
% A( : , 1) = ones(M , 1);
% for n = 1 : 1 : N - 1
%    A( : , n + 1) = cos(n * acos(xx')); 
% end
% y_num = A * a_0;
% plot(xx , y_num)

% y_num = clenshaw(xx' , a_1);

function cutoff = CHopping(coeffs, tol) 
% cutoff has two number
% for coeffs = (c_0 , c_2 , ... ,c_j , ... , c_j2 , ... , c_n)
% (c_j , ... , c_j2) is the plateau , then cutoff = [j ; j2])


if ( tol >= 1 ) % input not enough
    cutoff = [1 ; 1];
    return
end

n = length(coeffs);
cutoff = [2 * n ; 2 * n]; %assume there is no plateau
if ( n < 17 )  %length of coeffs is too short to build a plateau
    return
end

b = abs(coeffs);
m = b(end)*ones(n, 1);
for j = n-1:-1:1
    m(j) = max(b(j), m(j+1));
end   
if ( m(1) == 0 )
    cutoff = [1 ; 1];
    return
end
envelope = m/m(1);

for j = 2:n
    j2 = round(1.25*j + 5); 
    if ( j2 > n )
        % there is no plateau: exit
        return
    end      
    e1 = envelope(j);
    e2 = envelope(j2);
    r = 3*(1 - log(e1)/log(tol));
    plateau = (e1 == 0) | (e2/e1 > r);
    if ( plateau )
        cutoff = [j ; j2];
        break
    end
end
end