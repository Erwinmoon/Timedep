clc
clear
%% 
%Ŀ�� �� ���������׷����������·��̵ĵ���������������ֵ������ ��
%  du/dt = du/dx 
%  u(1 , t) = 0 
%  u(x , 0) = f(x) 
%  ����ʱ�䷽������һ����ʾŷ������ �� ʱ�䷽�������׷���


%%
% ��������
N = 63;
dt = 1 / N^2;
% dt = 0.0001;

%% 
% ��ֵ��ֵ����a_0 = chebcoeffs(chebfun('1 - x' , N + 1));

%%
% ΢�־���
D_0 = convertmat(N + 1, 1, 0) * diffmat(N + 1, 1);

%% 
% ����ת������
S_0 = convertmat(N + 1, 0, 0);
%
%% 
% �߽���������
P0 = eye(N + 1 , N + 1);
P0(N + 1 , 1 : N) = ones(1 , N);

%%
% ��������
P =  P0 ^ (-1) * S_0^(-1) * dt * D_0;
Q = P';

%% 
% ������ֵ
[xp, yp] = eig(P);
eigenvaluep = diag(yp);

[xq  , yq] = eig(Q);
eigenvalueq = diag(yq);

% max eigenvalue condition number
[ip, iq] = maxeigenvalue(eigenvaluep, eigenvalueq);

rvector = xp(: , ip);
lvector = xq(: , iq);

keigenmax = norm(lvector,2) * norm(rvector,2) / abs(lvector' * rvector);

% other eigenvalues condition number
keigen = zeros(5 , length(P));
for i = 1 : 1 : length(P)
    [ip, iq] = maxeigenvalue(eigenvaluep, eigenvalueq);
    if eigenvaluep(ip) == 0
        break;
    end
    rvector = xp(: , ip);
    lvector = xq(: , iq);
    keigen(1 , i) = real(ip);
    keigen(2 , i) = real(iq);
    keigen(3 , i) = eigenvaluep(ip);
    keigen(4 , i) = conj(eigenvalueq(iq));
    keigen(5 , i) = real(norm(xq(: , iq),2) * norm(xp(: , ip),2) / abs(xq(: , iq)' * xp(: , ip)));
    eigenvaluep(ip) = 0;
    eigenvalueq(iq) = 0;
end

%%
% ��������ֵ��������
function [ip, iq] = maxeigenvalue(eigenvaluep, eigenvalueq)
    rlamda = max(abs(eigenvaluep));
    for ip = 1 : 1 : length(eigenvaluep)
       if rlamda ==  abs(eigenvaluep(ip))
           break;
       end
    end
    llamda = max(abs(eigenvalueq));
    for iq = 1 : 1 : length(eigenvalueq)
       if llamda ==  abs(eigenvalueq(iq))
           break;
       end
    end
    if norm(conj(eigenvalueq(iq)) - (eigenvaluep(ip)), 2) > 1e-12
        eigenvalueq(iq) = 0;
        llamda = max(abs(eigenvalueq));
        for iq = 1 : 1 : length(eigenvalueq)
            if llamda ==  abs(eigenvalueq(iq))
                  break;
            end
        end
    end
end