% This is an example to illustrate the algorithm in "Resilient Reinforcement Learning and Robust Output Regulation Under Denial-of-Service Attacks"
% authored by Weinan Gao, Chao Deng, Yi Jiang, and Zhong-Ping Jiang published in Automatica, vol. 142, article 110366, 2022,


function [] = Gao_Auto2022()
clear; clc; close all;
%% Initialization
global A B C F K E G2 learn_index alpha attack_index attack_xi attack_e Azeta Bzeta Kzeta


D1 = 1;  H1 = 3; 
D2 = 1;  H2 = 3; 
w0 = 314;
A = [0 1; 0 -D1/(2 * H1)]; 
B = [0; w0/(2 * H1)];
C = [1 0 ];
F = -1;
E = 0;
G2 = 1;

Azeta = [0 1; 0 -D2/(2*H2)]; 
Bzeta = [0; w0/(2*H2)];
Kzeta = lqr(Azeta, Bzeta, 40*eye(2), 1);

Q = eye(3); R = eye(1); % Weights

alpha = 0.1;

barA = [ [A zeros(2,1)]; [G2 * C E]];
barB = [B; 0];
barC = [C 0];
[Ks,Ps] = lqr(barA,barB,Q,R);
K = zeros(1,3);

T_star =  (max(eig(Ps))*(4*norm(Ks)^2+4*norm(Ps*[zeros(2,3);[G2*barC]])))...
   /((min(eig(Q)))*min(eig(Ps)));

T_DOS = round(T_star)+1; % DoS Parameters

T=.01;  % length of each learning time interval

kappa=T; eta=1; % DoS Parameters

temp = norm(C)*sqrt(1+8*exp(kappa*T_star*min(eig(Q))/max(eig(Ps))+eta));
gamma_e = temp/ sqrt( min( min(eig(Q))* min(eig(Ps)) / max(eig(Ps)), (4*norm(Ks)^2+4*norm(Ps)*norm(G2*C))) ); % IOS gain


%% Initial conditions
x0 = [0;-1]; 
v0 = .065;
u0 = 0;
ze0 = [1;-1];
z0 = 0;
xi0 = [x0;z0];


%% DoS Settings
T_actual=50;

h=[];
for i=1:20
    h = [h; (i+2)*T_actual; (i+2)*T_actual+randi([1,15])];
end


    h = intersect( h(h<=1000), h(h>=100))* T;
    h = [19*T;20*T; 39*T; 40*T; 59*T; 60*T; 79*T; 80*T; 99*T; 100*T;h]; % attack transition   
    h2 = []; 

    h3 = [13	12	0	1	6	10	13	8	16	14	19	10	6	2	12	15	8	1	5	3	5	8 ...
             10	9	17	10	18	12	19	4	13	5	13	13	1	5	4	13	16	6	15	13 ...
              0	12	7	18	0	9	8	9];
for i=1:50
     h2 = [h2; i * T_DOS; i * T_DOS + h3(i)];
end
     h2 = intersect( h2(h2 <= 400), h2(h2 >= 300))* T;
     h = [h; h2];
    
learn_index = 1;

Dxixi_V = []; XiXi = []; UXi = []; VXi = [];
Xscope1 = []; Vscope1 = []; errorK = []; errorP = [];
P_1 = zeros(4);
Time1 = [];
ZeScope1 = []; 

%% Data Collection
for i=1:100
    
    attack_xi = [x0; z0];
    attack_e = C * x0+F * v0;
    attack_index = 0;
    if (rem(i,20) == 0)
    attack_index = 1;
    end
    
    [t,X] = ode45(@mysys,[(i-1) * T,i * T],[x0; z0; v0; ze0; kron(xi0',xi0')'; kron(xi0,u0); kron(xi0',v0')']); 
    if (attack_index == 0)
        Xi = X(end,1:3);
        Xi_pre = X(1,1:3);
        Dxixi_V = [Dxixi_V; kron(Xi,Xi) - kron(Xi_pre,Xi_pre)];
        XiXi = [XiXi; X(end,6+1:6+9) - X(1,6+1:6+9)];
        UXi = [UXi; X(end,15+1:15+3) - X(1,15+1:15+3)];
        VXi = [VXi; X(end,18+1:18+3) - X(1,18+1:18+3)];
    end
    x0 = X(end,1:2)';
    z0 = X(end,3)';
    v0 = X(end,4)';
    ze0 = X(end,5:6)';
    Xscope1 = [Xscope1; X(:,1:3)];
    Vscope1 = [Vscope1; X(:,4)];
    ZeScope1 = [ZeScope1; X(:,5:6)];
    Time1 = [Time1; t];
end
Dxixi = Dxixi_V(:, [1,2,3,5,6,9]); 
XiXi_V = XiXi(:, [1,2,3,5,6,9]); 
if (rank([Dxixi, UXi, VXi])<12)
    rank([Dxixi, UXi, VXi])
error('The rank conditon is not satisfied!')
end

%% Hybrid Iteration
errorP = []; errorK = [];
q = 0;
P0 = 0.001*eye(3);
P1 = P0;
epsilonk = 0.05;
k_underline = 1;
k_q_bar = (2 * (q+1) / epsilonk / min(eig(Q)) )^2 + 2;

for i=1:500
        
    Y = Dxixi_V * P1(:);
    X = [XiXi_V,2 * UXi,2 * VXi];
    pp = inv(X' * X) * X' * Y;
    K = [pp(7:9)]';   
    H = [pp(1) pp(2)/2 pp(3)/2; 
       pp(2)/2 pp(4) pp(5)/2;
       pp(3)/2 pp(5)/2 pp(6)];   
    P1_1 = P1 + epsilonk * (H+Q-R*K'*K);
    errorP=[errorP; norm(P1-Ps) ];
    errorK=[errorK; norm(K-Ks) ];
    if (norm(P1_1) >= 2 *(q+1)) || (min(eig(P1_1))<0)
    P1_1 = P0;
    q = q+1;
    epsilonk = epsilonk / 2; 
    k_underline = i;
    k_q_bar = (2 * (q+1) / epsilonk / min(eig(Q)) )^2 + 2;
    display('q is updated by');
        q
    end

    if (min( min( real( eig(2 * K' * K - H))), min( real( eig(P1))) ) > 0) || (i > k_underline + k_q_bar + 1)         
        iter_VI = i
        break;
    end 
P1 = P1_1;
end

for i = iter_VI : iter_VI+1e3
    QK = Q + K'*R*K;
    Y = -XiXi * QK(:);
    X = [Dxixi,-XiXi * kron(eye(3),K') - UXi,VXi];
    pp = pinv(X)*Y;
    P1 = [pp(1) pp(2)/2 pp(3)/2; 
       pp(2)/2 pp(4) pp(5)/2;
       pp(3)/2 pp(5)/2 pp(6)];
    K = inv(R) * [pp(7:9)]' / 2;
    errorK = [errorK; norm(K-Ks)];
    errorP = [errorP; norm(P1-Ps)];

if (norm(P1 - P1_1)<0.5) 
    disp('meet the converngent requirement and the norm of the P_k-P_{k-1} is')
    i
    disp('The learned control gain and value are')
    vpa(K,6)
    P1
    disp('The optimal control gain and value are')
    vpa(Ks,6)
    Ps
    break
    
end 
P1_1 = P1;
end

index = [1:1:iter_VI iter_VI+1:length(errorP)]; 

figure;
subplot(2,1,1); plot(index,errorP(index),'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
xlabel('Number of Iterations','Fontsize',10);
h1 = legend('$|P_k-P^*|$'); set(h1,'Interpreter','latex')
subplot(2,1,2); plot(index,errorK(index),'--rs','LineWidth',2, 'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
xlabel('Number of Iterations','Fontsize',10);
h1 = legend('$|K_k-K^*|$'); set(h1,'Interpreter','latex')


learn_index = 0;
for i=101:1000
    attack_xi = [x0; z0];
    attack_e = C * x0 + F * v0;
    attack_index = 0;
    if (rem(i, T_DOS)==1)
    attack_index = 1;
    end
    [t, X] = ode45(@mysys,[(i-1)*T,i*T],[x0;z0;v0;ze0;kron(xi0',xi0')';kron(xi0,u0);kron(xi0',v0')']);
    x0 = X(end,1:2)';
    z0 = X(end,3)';
    v0 = X(end,4)';
    ze0 = X(end,5:6)';
    Xscope1 = [Xscope1; X(:,1:3)];
    Vscope1 = [Vscope1; X(:,4)];
    ZeScope1 = [ZeScope1; X(:,5:6)];
    Time1 = [Time1; t];
end

   figure;

   subplot(2,1,1)
   for i=1:length(h)/2
   X = [h(2*i-1), h(2*i), h(2*i), h(2*i-1)];  
   Y = [-0.5, -0.5, 1.5, 1.5];  
   fill(X,Y,[0.85, 0.85, 0.85]); hold on;
   end
   plot(Time1, Xscope1(:,1) / 100 + F * Vscope1(:,1) / 100,'LineWidth', 2);
   hold on; 
   grid;
   xlabel('Time(sec)');
   ylabel('e(t)');
   
   subplot(2,1,2)
   plot(Time1, ZeScope1(:,1), 'LineWidth', 2);
   hold on;
   grid;
   xlabel('Time(sec)');
   ylabel('\zeta_1(t)');   
end



%% System subfunction
function dX=mysys(t,X)
global A B C F K E G2 learn_index attack_index attack_xi attack_e Azeta Bzeta Kzeta
E1 = 1.2; E2 = 5; X_E = 15;

x = X(1 : 2);
z = X(3);
v = X(4);
ze = X(5 : 6);
xi = [x; z];

noise = .5 * sum(sin([1:1:100]*t)); 

if (learn_index == 0)
noise = 0;
end

if (attack_index == 0)
    e = C * x + F * v;
    u=-K * xi + noise;
end

if(attack_index == 1)
    e = attack_e;
    u = -K * attack_xi + noise; 
end


Phi = -E1 * E2 / X_E * (sin(x(1)-ze(1))-sin(v));
dx = A * x + B * (u + Phi);
dz = E * z + G2 * e;
w = u + Phi;
dv = E * v;
dze = (Azeta - Bzeta * Kzeta) * ze - Bzeta * Phi;
dxixi = kron(xi',xi')';
duxi = kron(xi',w')';
dvxi = kron(xi',v')';

dX = [dx;  %2
    dz; % 1
    dv; %1
    dze; % 2
    dxixi; % 9
    duxi; %3
    dvxi; %3
    ];
end