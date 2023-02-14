clc
clear
%% 
%Ŀ�� �� ���������׷����������·��̵ĵ���������������ֵ������ ��
%  du/dt = d2u/dx2 
%  u(+-1 , t) = 0 
%  u(x , 0) = f(x) 
%  ����ʱ�䷽������һ����ʾŷ������ �� ʱ�䷽�������׷���


%% 
% ��������
N = 63;
% dt =  03 / N^(4);
dt = 1/N^4;

%% 
% ΢�־���
D_2 = convertmat(N + 1, 2, 1) * diffmat(N + 1, 2);

%% 
% ����ת������
S = convertmat(N + 1, 0, 1);

%% 
% ��ֵ����
B = [bc(N + 1, 'd', 'l'); bc(N+1 , 'd', 'r')];
BC = [1 ; 1];
P0 = eye(N + 1 , N + 1);
P0(N : N + 1 , 1 : N + 1) = B;

%% 
% ���
P = dt * P0 ^(-1) * S^(-1) * D_2;
Q = P';

%% 
% ������ֵ
[xp, yp] = eig(P);
eigenvalue = diag(yp);

rlamda = max(abs(eigenvalue));
for i = 1 : 1 : length(P)
   if rlamda ==  abs(eigenvalue(i))
       break;
   end
end
eigenvalue(i) = 0;

rlamda = max(abs(eigenvalue));
for i = 1 : 1 : length(P)
   if rlamda ==  abs(eigenvalue(i))
       break;
   end
end

rvector = xp(: , i);

% ������ֵ
[xq  , yq] = eig(Q);
eigenvalue = diag(yq);
llamda = max(abs(eigenvalue));

for j = 1 : 1 : length(Q)
   if llamda ==  abs(eigenvalue(j))
       break;
   end
end
eigenvalue(j) = 0;

llamda = max(abs(eigenvalue));
for j = 1 : 1 : length(P)
   if llamda ==  abs(eigenvalue(j))
       break;
   end
end
lvector = xq(: , j);

keigen = norm(lvector,2) * norm(rvector,2) / abs(lvector' * rvector)

