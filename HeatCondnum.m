clc;
clear;
%% 
%Ŀ�� �� ���Թ������u_t = u_x
%% 
% ��������
N = 30;
dt = 0.00000001;

E = eye(N + 1);
O = zeros(N + 1 , N + 1);
BC = [bc(N + 1, 'd', 'l'); bc(N + 1, 'd', 'r')];
%% ------------------RK3-------------------------------------------
% ΢�־���
Drk = convertmat(N + 1, 1, 0) * diffmat(N + 1, 2);

% ����ת������
Srk = convertmat(N + 1, 0, 1);
Prk = eye(N + 1 , N + 1);
Prk(N :N + 1 , 1 : N + 1) = BC;

Ark = Srk * Prk;
Brk = Drk * dt;

a = norm(Ark^(-1) , inf)
b = cond(Ark , inf)
c = norm(Brk , inf)
d = KreissConstant(Ark^(-1) * Brk)

    
%% ----------------------BDF3--------------------------------------
% ΢�־���
Dbdf = convertmat(N + 1, 1, 0) * diffmat(N + 1, 2);

% ����ת������
Sbdf = convertmat(N + 1, 0, 1);

% �߽���������
Pbdf = eye(N + 1 , N + 1);
Pbdf(N :N + 1 , 1 : N + 1) = BC;

Sbdf = Sbdf * Pbdf;

Abdf = [Sbdf - 6 * dt * Dbdf , E , E;O , E , O; O, O , E];
Bbdf = [18*Sbdf , -9*Sbdf , 2*Sbdf ;E , O , O ; O , E , O];

a = norm(Abdf^(-1) , inf)
b = cond(Abdf , inf)
c = norm(Bbdf , inf)  
d = KreissConstant(Abdf^(-1) * Bbdf)
%% ---------------------AB4-----------------------------------
% ΢�־���
Dab = convertmat(N + 1, 1, 0) * diffmat(N + 1, 2);

% ����ת������
Sab = convertmat(N + 1, 0, 1);

% �߽���������
Pab = eye(N + 1 , N + 1);
Pab(N :N + 1 , 1 : N + 1) = BC;

Sab = Sab * Pab;
Dab = dt * Dab;

Aab = [Sab , O , O , O ; O , E , O ,O; O , O , E , O; O , O , O , E];
Bab = [Sab + (55/24) * Dab , -(59/24) * Dab , (37 / 24) * Dab , -(9 / 24) * Dab ; ... 
    E , O , O , O ; O , E , O , O; O , O , E , O];

a = norm(Aab^(-1) , inf)
b = cond(Aab , inf)
c = norm(Bab , inf) 
d = KreissConstant(Aab^(-1) * Bab)


function ka = KreissConstant(A)
    ka = 1;
    for power = 1 : 0.5 : 16
        eps = exp(-power);
        ka = max(abs((pspr_2way(A , eps) - 1) / eps) , ka);
    end      
end