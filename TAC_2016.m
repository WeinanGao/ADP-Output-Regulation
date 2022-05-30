% This is an example to illustrate the algorithm in "Adaptive dynamic programming and adaptive
% optimal output regulation of linear systems," by Weinan Gao and Zhong-Ping Jiang
% published in IEEE Transactions on Automatic Control, 2016.
function []=TAC_2016()
clear all; clc;  close all;
%% initialization
global A B K E Xstar Ustar X_trial
A=[0 1;0 0];  B=[0; 0.1];  C=[1 0];        E=[0 1;-1 0];       F=[-1 0];
Q=100*eye(2);     R=eye(1);
K=[30 40]; % initial control policy
Xscope1=[]; Vscope1=[]; Xscope2=[]; errorK=[]; errorP=[]; P_1=zeros(2);
Vscope2=[]; Time1=[]; Time2=[]; EE=[]; UE=[]; VE=[]; Dee=[];
x0=[1;2];       v0=[1;2];     u0=0;       %initial condition
T=.1;  % length of each learning time interval
%% Select trials of regulator eq.: X_trial(:,:,i)="X_i" in the paper
X_trial(:,:,1)=[1 0;0 0];  % C*X_trial(:,:,1)+F=0
X_trial(:,:,2)=[0 0;1 0];  % C*X_trial(:,:,2)=0
X_trial(:,:,3)=[0 0;0 1];  % C*X_trial(:,:,3)=0
%% Collect data
for i=1:20
    [t,X]=ode45(@mysys,[(i-1)*T,i*T],[x0;v0;zeros(30,1)]); %34   
    for j=1:size(X_trial,3)
        E_cur=X(end,1:2)-X(end,3:4)*X_trial(:,:,j)';
        E_pre=X(1,1:2)-X(1,3:4)*X_trial(:,:,j)';
        Dee(i,:,j)=kron(E_cur,E_cur)-kron(E_pre,E_pre);
        EE(i,:,j)=X(end,10*(j-1)+5:10*(j-1)+8)-X(1,10*(j-1)+5:10*(j-1)+8);
        UE(i,:,j)=X(end,10*(j-1)+9:10*(j-1)+10)-X(1,10*(j-1)+9:10*(j-1)+10);
        VE(i,:,j)=X(end,10*(j-1)+11:10*(j-1)+14)-X(1,10*(j-1)+11:10*(j-1)+14);
    end
    x0=X(end,1:2)';    v0=X(end,3:4)';
    Xscope1=[Xscope1;X(:,1:2)];
    Vscope1=[Vscope1;X(:,3:4)];
    Time1=[Time1;t];
end

for k=1:size(X_trial,3)
    if (rank([EE(:,:,k), UE(:,:,k), VE(:,:,k)]) < 9)
        error('The rank conditon is not satisfied!')
    end
end

Dee=Dee(:,[1,2,4],:);
[Kstar,Pstar]=lqr(A,B,Q,R); % Optimal values

%% Compute Apporixmate Optimal Cost and Control Gains
for j=1:30
    QK=Q+K'*R*K;
    Y=-EE(:,:,1)*QK(:);
    X=[Dee(:,:,1),-EE(:,:,1)*kron(eye(2),K')-UE(:,:,1),VE(:,:,1)];
    pp=inv(X'*X)*X'*Y;
    j
    P=[pp(1) pp(2)/2; pp(2)/2 pp(3)]
    K=inv(R)*[pp(4:5)]'/2

    errorK=[errorK;norm(K-Kstar)];
    errorP=[errorP;norm(P-Pstar)];

    if (norm(P-P_1)<1e-2) 
        disp('The solution to the Riccati eq')
            Pstar
            Kstar
  %% find Sylvester map: S(X_i)
        SX1=([pp(6:7) pp(8:9)]/2/P)'; %X_trial(:,:,1)*E-A*X_trial(:,:,1)  
        
        Y=-EE(:,:,2)*QK(:);
        X=[Dee(:,:,2),-EE(:,:,2)*kron(eye(2),K')-UE(:,:,2),VE(:,:,2)];
        pp=inv(X'*X)*X'*Y;
        newP=[pp(1) pp(2)/2; pp(2)/2 pp(3)];
        SX2=([pp(6:7) pp(8:9)]/2/newP)';%        X_trial(:,:,2)*E-A*X_trial(:,:,2)
        
        Y=-EE(:,:,3)*QK(:);
        X=[Dee(:,:,3),-EE(:,:,3)*kron(eye(2),K')-UE(:,:,3),VE(:,:,3)];
        pp=inv(X'*X)*X'*Y;
        newP2=[pp(1) pp(2)/2; pp(2)/2 pp(3)];
        SX3=([pp(6:7) pp(8:9)]/2/newP2)';%        X_trial(:,:,3)*E-A*X_trial(:,:,3)
        break
    end 
    P_1=P;
end
%% The solution of the regulator eq
    newB=inv(P)'*K'*R;  
    M=[SX2(1,1) SX3(1,1) -newB(1) 0
       SX2(1,2) SX3(1,2) 0 -newB(1)
       SX2(2,1) SX3(2,1) -newB(2) 0
       SX2(2,2) SX3(2,2) 0 -newB(2)];
    N=-[SX1(1,1) SX1(1,2) SX1(2,1) SX1(2,2)]';
    capX=M\N;    
    disp('The solution to the regulator eq')
    Xstar=X_trial(:,:,1)+capX(1)*X_trial(:,:,2)+capX(2)*X_trial(:,:,3)
    Ustar=capX(3:4,:)'
   
    for i=21:100
    [t,X]=ode45(@mysys2,[(i-1)*T,i*T],[x0;v0]);
    x0=X(end,1:2)';
    v0=X(end,3:4)';
    Xscope2=[Xscope2;X(:,1:2)];
    Vscope2=[Vscope2;X(:,3:4)];
    Time2=[Time2;t];
    end

    %% Plot figures
    Xscope=[Xscope1;Xscope2];
    Vscope=[Vscope1;Vscope2];  
    Time=[Time1;Time2];
    plot(Time,Xscope(:,1),Time,Vscope(:,1),'--r','LineWidth',2);
    hold on; 
    xlabel('Time');
    h = legend('Output','Reference');
    grid;
    annotation(gcf,'textarrow',[0.370535714285714 0.291666666666667],...
    [0.623015873015873 0.579365079365079],'TextEdgeColor','none',...
    'String',{'Controller Updated'});
    
    figure;
    subplot(1,2,1);plot(1:length(errorP),errorP,'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
    xlabel('Number of Iteration','Fontsize',10);    
    h = legend('$\|\bar P_j-\bar P^*\|$'); set(h,'Interpreter','latex')
    subplot(1,2,2);plot(1:length(errorK),errorK,'--rs','LineWidth',2, 'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
    xlabel('Number of Iteration','Fontsize',10);   
    h = legend('$\|\bar K_j-\bar K^*\|$'); set(h,'Interpreter','latex')
end

%% System subfunction
function dX=mysys(t,X)
global A B K E X_trial
x=X(1:2);  v=X(3:4);
e0=x-X_trial(:,:,1)*v; e1=x-X_trial(:,:,2)*v; e2=x-X_trial(:,:,3)*v;
u=-K*x+sin(10*t); % sin(10*t) is an exploration noise
dx=A*x+B*u;
dv=E*v;
de0e0=kron(e0',e0')'; due0=kron(e0',u')'; dve0=kron(e0',v')';
de1e1=kron(e1',e1')'; due1=kron(e1',u')'; dve1=kron(e1',v')';
de2e2=kron(e2',e2')'; due2=kron(e2',u')'; dve2=kron(e2',v')';
dX=[dx; dv; 
    de0e0; due0;  dve0; 
    de1e1; due1;  dve1;
    de2e2; due2;  dve2; 
    ];
end

function dX=mysys2(t,X)
global A B K E Xstar Ustar
x=X(1:2); v=X(3:4);
u=-K*(x-Xstar*v)+Ustar*v;
dx=A*x+B*u;   dv=E*v;
dX=[dx; dv];
end
