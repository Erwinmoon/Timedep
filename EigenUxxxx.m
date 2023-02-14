clc
clear
clf

figure(1)
set(gcf , 'Color' , 'w')
lw = 1;
FS = 'FontSize'; fs = 16;
IN = 'interpret'; LX = 'latex';

%% 
%Ŀ�� �� ���������׷����������·��̵ĵ���������׺�α�� ��
%  du/dt = du/dx 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  ����ʱ�䷽������һ����ʾŷ������ �� ʱ�䷽�������׷���


%%
% ��������
j = 1;
for N = 10 : 10 : 1000
    dt = 1;
    % dt = 0.0001;

    %%
    % ΢�־���
    D = diffmat(N + 1, 4);

    %% 
    % ����ת������
    S = convertmat(N + 1, 0, 3);
    %
    %% 
    % �߽���������
    P = eye(N + 1 , N + 1);
    B = [bc(N + 1, 'd', 'l'); bc(N+1 , 'd', 'r');bc(N, 'n', 'l');bc(N, 'n', 'r')];
    P(N - 2 : N + 1 , :) = B;
    %%
    % ��������
    Q =  P ^ (-1) * S^(-1) * dt * D;

    %% 
    % �����������
    lamda = eig(Q);
    
    %%
    % ��ͼ����
    C(j) = max(abs(lamda))/N^8;
    NN(j) = (N);
    j = j + 1;
end

plot(NN , C , 'k' , 'LineWidth',1)
xlabel('n')
ylabel('$\frac{\mathbf{\rho}}{\mathrm{n}^8}$', IN , LX ,'rotation',0, FS ,fs ,...
    'position' , [-120 0.016])
axis([10,1000,0.008,0.024])
set(gca,'yTick',[0.008:0.004: 0.024])
pbaspect([2 1 1])
