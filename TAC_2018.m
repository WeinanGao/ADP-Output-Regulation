% This is an example to illustrate the algorithm in 
% Weinan Gao, Zhong-Ping Jiang, Frank L. Lewis, and Yebin Wang.  
% Leader-to-formation stability of multi-agent systems: An adaptive optimal control approach
% IEEE Transactions on Automatic Control, vol. 63, no. 10, pp. 3581-3587, 2018. doi: 10.1109/TAC.2018.2799526.

function []=TAC_2018()
clear all;
clc;
close all;
% initialization
global bA bB bD bD2 bK0 E Ke H Hrho BKlearn mu

seed=100*rand(5,100)-50;

A1=[0 1;0 0]; B1=[0 1]'; C1=[1 0]; D1=[0 0;0 0.5]; F1=[-2 0];
A2=[0 1;0 0]; B2=[0 1]'; C2=[1 0]; D2=[0 0;0 1];  F2=[-3 0];
A3=[0 1;0 0]; B3=[0 1]'; C3=[1 0]; D3=[0 0;0 1.5]; F3=[-4 0];
A4=[0 1;0 0]; B4=[0 1]'; C4=[1 0]; D4=[0 0;0 2];  F4=[-5 0];

E=[0 1;-1 0]; H=[1 1]';  Ke=[1.5,0.5];
G1=E; G2=[0 1]';

bA1=[A1 zeros(2); G2*C1 G1]; bB1=[B1;zeros(2,1)]; bC1=[C1 0 0]; 
bA2=[A2 zeros(2); G2*C2 G1]; bB2=[B2;zeros(2,1)]; bC2=[C2 0 0]; 
bA3=[A3 zeros(2); G2*C3 G1]; bB3=[B3;zeros(2,1)]; bC3=[C3 0 0]; 
bA4=[A4 zeros(2); G2*C4 G1]; bB4=[B4;zeros(2,1)]; bC4=[C4 0 0];
bK10=[0.2 1 -0.1 -0.1];
%bK10=lqr(bA1,bB1,eye(4),1);
bA=[bA1 -1/2*[zeros(2,4); G2*C2 zeros(2)] zeros(4) zeros(4);
         zeros(4) bA2 zeros(4) zeros(4);
         zeros(4) -[zeros(2,4); G2*C2 zeros(2)] bA3 zeros(4);
         -1/3*[zeros(2,4); G2*C1 zeros(2)] -1/3*[zeros(2,4); G2*C2 zeros(2)] -1/3*[zeros(2,4); G2*C3 zeros(2)] bA4];
bB=blkdiag(bB1,bB2,bB3,bB4);
bD=blkdiag([D1;zeros(2)],[D2;zeros(2)],[D3;zeros(2)],[D4;zeros(2)]);
bD2=blkdiag([zeros(2);G2*F1],[zeros(2);G2*F2],[zeros(2);G2*F3],[zeros(2);G2*F4]);
bK0=blkdiag(bK10,bK10,bK10,bK10);

H1=[1 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0];

H2=[0 0 0 0
    0 1 0 0 
    0 -1 1 0
    0 0 0 0 ];

H3=[0 0 0 0
    0 0 0 0
    0 0 0 0
    -1 -1 0 2];

mu=1;

N=50; % number of samples
T=.2;
xi0=(-8:7)';
zeta0=[1;2;1;2;1;2;1;2];
v0=[1 1]'; vstar0=[5 0]';

E1E1=[]; E1U=[]; E1V=[]; E1E2=[];
E2E2=[]; E2U=[]; E2V=[]; 
E3E3=[]; E3U=[]; E3V=[]; E3E2=[];
E4E4=[]; E4U=[]; E4V=[]; E4E1=[]; E4E2=[]; E4E3=[];
Dee1=[]; Dee2=[]; Dee3=[]; Dee4=[];
Xscope1=[]; Vscope1=[]; Time1=[]; Uscope1=[]; ZetaScope1=[]; VstarScope1=[];
for time_i=1:N
    % Ini=[x0;v0;vstar0;zeros(1,16*16)'; zeros(1,16*4)'; zeros(1,16*2)']; %XiXi UXi VXi
    Ini=[xi0;zeta0;v0;vstar0;zeros(1,36+36+36+36)']; 
    if (mod(time_i,3)==1)
        Hrho=H1;
    end
    if (mod(time_i,3)==2)
        Hrho=H2;
    end
    if (mod(time_i,3)==0)
        Hrho=H3;
    end
    [t,X]=ode45(@mysys,[(time_i-1)*T,time_i*T], Ini);
     E0=X(end,1:4); E0_pre=X(1,1:4);  Dee1=[Dee1;kron(E0,E0)-kron(E0_pre,E0_pre)];
     E0=X(end,5:8); E0_pre=X(1,5:8);  Dee2=[Dee2;kron(E0,E0)-kron(E0_pre,E0_pre)];
     E0=X(end,9:12); E0_pre=X(1,9:12);  Dee3=[Dee3;kron(E0,E0)-kron(E0_pre,E0_pre)];
     E0=X(end,13:16); E0_pre=X(1,13:16);  Dee4=[Dee4;kron(E0,E0)-kron(E0_pre,E0_pre)];
%     
     E1E1=[E1E1;X(end,28+1:28+16)-X(1,28+1:44)];     E1U=[E1U;X(end,44+1:44+4)-X(1,44+1:48)];
     E1V=[E1V;X(end,48+1:48+16)-X(1,48+1:64)];    
%     
     E2E2=[E2E2;X(end,64+1:64+16)-X(1,64+1:80)];     E2U=[E2U;X(end,80+1:80+4)-X(1,80+1:84)];
     E2V=[E2V;X(end,84+1:84+16)-X(1,84+1:100)];    
%     
     E3E3=[E3E3;X(end,100+1:100+16)-X(1,100+1:116)];     E3U=[E3U;X(end,116+1:116+4)-X(1,116+1:120)];
     E3V=[E3V;X(end,120+1:120+16)-X(1,120+1:136)];   
%     
     E4E4=[E4E4;X(end,136+1:136+16)-X(1,136+1:152)];    E4U=[E4U;X(end,152+1:152+4)-X(1,152+1:156)];
     E4V=[E4V;X(end,156+1:156+16)-X(1,156+1:end)];  
    
    xi0=X(end,1:16)';  zeta0=X(end,17:24)';  v0=X(end,25:26)';  vstar0=X(end,27:28)';
    Xscope1=[Xscope1;X(:,1:16)];    ZetaScope1=[ZetaScope1;X(:,17:24)]; Vscope1=[Vscope1;X(:,25:26)]; 
    VstarScope1=[VstarScope1; X(:,27:28)];    Time1=[Time1;t];
    noise=[];
    for i=1:length(t)
 noise =[noise; 0.3*[sum(sin([38.1558   76.5517   79.5200   18.6873]*t(i)));
 sum(sin([17.9703   15.5098  -33.7388  -38.1002   100   45.9744  -15.9614    8.5268  -27.6188   25.1267]*t(i)));
 sum(sin([7.9703   5.5098  -3.7388  -8.1002   10   5.9744  -5.9614    .5268  -7.6188   5.1267]*t(i)));
 sum(sin([137.9703   125.5098  -323.7388  -8.1002   1000   425.9744  -.9614    .5268  -.6188   .1267]*t(i)))]'];
    end
    Uscope1=[Uscope1;[-bK0*X(:,1:16)'+noise']'];
end

if (rank([E1E1 E1U E1V])<30)
    error('The rank condition 1 does not hold!');
end

if (rank([E2E2 E2U E2V])<30)
    error('The rank condition 2  does not hold!');
end

if (rank([E3E3 E3U E3V])<30)
    error('The rank condition 3  does not hold!');
end

if (rank([E4E4 E4U E4V])<30)
    error('The rank condition 4  does not hold!');
end

Dee1=Dee1(:,[1,2,3,4,6,7,8,11,12,16]); 
Dee2=Dee2(:,[1,2,3,4,6,7,8,11,12,16]); 
Dee3=Dee3(:,[1,2,3,4,6,7,8,11,12,16]); 
Dee4=Dee4(:,[1,2,3,4,6,7,8,11,12,16]);

bK20=bK10; bK30=bK10; bK40=bK10; P_1=100*eye(4);
errorP1=[]; errorP2=[]; errorP3=[]; errorP4=[];

%% i=1
[Ks,Ps,Es]=lqr(bA1,bB1,eye(4),1); 
for k=1:10
QK=eye(4)+bK10'*bK10;
 Y=-E1E1*QK(:);
 X=[Dee1,-E1E1*kron(eye(4),bK10')-E1U,-E1V];
 pp=(X'*X)\X'*Y;
 Ph=lyap((bA1-bB1*bK10)',QK)
 P1=[pp(1) pp(2)/2 pp(3)/2 pp(4)/2; 
       pp(2)/2 pp(5) pp(6)/2 pp(7)/2;
       pp(3)/2 pp(6)/2 pp(8) pp(9)/2;
       pp(4)/2 pp(7)/2 pp(9)/2 pp(10)]
    bK10=[pp(11:14)]'/2;
    errorP1=[errorP1; norm(P1-Ps) ];
    if (norm(P1-P_1)<1e-3) 
        
        disp('meet the converngent requirement and the norm of the P_1^k-P_1^{k-1} is')
        k
        break
    end
    P_1=P1;
end
 

%% i=2
[Ks,Ps,Es]=lqr(bA2,bB2,eye(4),1);
P_1=100*eye(4);
for k=1:10
QK=eye(4)+bK20'*bK20;
 Y=-E2E2*QK(:);
 X=[Dee2,-E2E2*kron(eye(4),bK20')-E2U,-E2V];
 pp=(X'*X)\X'*Y;
 Ph=lyap((bA2-bB2*bK20)',QK)
 P1=[pp(1) pp(2)/2 pp(3)/2 pp(4)/2; 
       pp(2)/2 pp(5) pp(6)/2 pp(7)/2;
       pp(3)/2 pp(6)/2 pp(8) pp(9)/2;
       pp(4)/2 pp(7)/2 pp(9)/2 pp(10)]
    errorP2=[errorP2; norm(P1-Ps) ];   
    bK20=[pp(11:14)]'/2
    if (norm(P1-P_1)<1e-3) 
        
        disp('meet the converngent requirement and the norm of the P_2^k-P_2^{k-1} is')
        k
         
        break
    end
    P_1=P1;
end
% 
%% i=3
 [Ks,Ps,Es]=lqr(bA3,bB3,eye(4),1);
P_1=100*eye(4);
for k=1:10
QK=eye(4)+bK30'*bK30;
 Y=-E3E3*QK(:);
 X=[Dee3,-E3E3*kron(eye(4),bK30')-E3U,-E3V];
 pp=(X'*X)\X'*Y;
 Ph=lyap((bA3-bB3*bK30)',QK);
 P1=[pp(1) pp(2)/2 pp(3)/2 pp(4)/2; 
       pp(2)/2 pp(5) pp(6)/2 pp(7)/2;
       pp(3)/2 pp(6)/2 pp(8) pp(9)/2;
       pp(4)/2 pp(7)/2 pp(9)/2 pp(10)];
   errorP3=[errorP3; norm(P1-Ps) ];
    bK30=[pp(11:14)]'/2
    if (norm(P1-P_1)<1e-3) 
        
        disp('meet the converngent requirement and the norm of the P_3^k-P_3^{k-1} is')
        k
       
        break
    end
    P_1=P1;
end


%% i=4
 [Ks,Ps,Es]=lqr(bA4,bB4,eye(4),1);
P_1=100*eye(4);
for k=1:20
QK=eye(4)+bK40'*bK40;
 Y=-E4E4*QK(:);
 X=[Dee4,-E4E4*kron(eye(4),bK40')-E4U,-E4V];
 pp=(X'*X)\X'*Y;
 Ph=lyap((bA4-bB4*bK40)',QK);
 P1=[pp(1) pp(2)/2 pp(3)/2 pp(4)/2; 
       pp(2)/2 pp(5) pp(6)/2 pp(7)/2;
       pp(3)/2 pp(6)/2 pp(8) pp(9)/2;
       pp(4)/2 pp(7)/2 pp(9)/2 pp(10)];
  errorP4=[errorP4; norm(P1-Ps) ];
    bK40=[pp(11:14)]'/2
    if (norm(P1-P_1)<1e-3) 
        
        disp('meet the converngent requirement and the norm of the P_4^k-P_4^{k-1} is')
        k

        break
    end
    P_1=P1;
end
% 
 BKlearn=blkdiag(bK10,bK20,bK30,bK40);

%% Noise off

for time_i=N+1:3*N
    Ini=[xi0;zeta0;v0;vstar0]; 
    if (mod(time_i,3)==1)
        Hrho=H1;
    end
    if (mod(time_i,3)==2)
        Hrho=H2;
    end
    if (mod(time_i,3)==0)
        Hrho=H3;
    end
    [t,X]=ode45(@mysys2,[(time_i-1)*T,time_i*T], Ini);
    xi0=X(end,1:16)';  zeta0=X(end,17:24)';  v0=X(end,25:26)';  vstar0=X(end,27:28)';
    Xscope1=[Xscope1;X(:,1:16)];    ZetaScope1=[ZetaScope1;X(:,17:24)]; Vscope1=[Vscope1;X(:,25:26)]; 
    VstarScope1=[VstarScope1; X(:,27:28)];    Time1=[Time1;t];
    Uscope1=[Uscope1;[-BKlearn*X(:,1:16)']'];
end

%% Figures

   plot(Time1,Uscope1(:,1),'r','LineWidth',2);hold on;
   plot(Time1,Uscope1(:,2),'g','LineWidth',2);hold on; 
   plot(Time1,Uscope1(:,3),'b','LineWidth',2);hold on;
   plot(Time1,Uscope1(:,4),'k','LineWidth',2);hold on;
   legend('u_1','u_2','u_3','u_4');
   xlabel('Time(s)','Fontsize',10);
 
   figure;
   h2=axes('position',[0.1 0.1 0.85 0.85]);
   axis(h2);
   grid;
   plot(Time1,Xscope1(:,1),'r','LineWidth',2);hold on;
   plot(Time1,Xscope1(:,5),'g','LineWidth',2);hold on; 
   plot(Time1,Xscope1(:,9),'b','LineWidth',2);hold on;
   plot(Time1,Xscope1(:,13),'k','LineWidth',2);hold on;
   plot(Time1,Vscope1(:,1),'m','LineWidth',2);hold on;
   plot(Time1,-F1(1)*VstarScope1(:,1),'r:','LineWidth',2);hold on;
   plot(Time1,-F2(1)*VstarScope1(:,1),'g:','LineWidth',2);hold on;
   plot(Time1,-F3(1)*VstarScope1(:,1),'b:','LineWidth',2);hold on;
   plot(Time1,-F4(1)*VstarScope1(:,1),'k:','LineWidth',2);hold on;
   plot(Time1,VstarScope1(:,1),'m:','LineWidth',2);hold on;
   %axis([0 30 -20 35]);
   legend('y_1','y_2','y_3','y_4','y_0','y_1^*','y_2^*','y_3^*','y_4^*','y_0^*');
   xlabel('Time(s)','Fontsize',10);
%    hold on;
%    h1=axes('position',[0.45 0.65 0.3 0.3]);
%    axis(h1);
%    plot(Time1,Xscope1(:,1),'r','LineWidth',2);hold on;
%    plot(Time1,Xscope1(:,5),'g','LineWidth',2);hold on; 
%    plot(Time1,Xscope1(:,9),'b','LineWidth',2);hold on;
%    plot(Time1,Xscope1(:,13),'k','LineWidth',2);hold on;
%    plot(Time1,Vscope1(:,1),'m-.','LineWidth',2);hold on;
%    plot(Time1,VstarScope1(:,1),'c:','LineWidth',2);hold on;
%    axis([24 26 1 6]);

% Create textarrow
annotation('textarrow',[0.480745223477982 0.397454031117397],...
    [0.178195260976718 0.249480249480249],'TextEdgeColor','none',...
    'Interpreter','latex',...
    'String',{'Control Policy Updated'});

%    
 figure;
 subplot(2,2,1);plot(1:length(errorP1),errorP1,'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
 xlabel('Number of Iteration','Fontsize',10);
 h = legend('$\|P_{1}^{(k)}-P_1^*\|$'); set(h,'Interpreter','latex')
%  subplot(2,4,2);plot(1:length(errorKA),errorKA,'--rs','LineWidth',2, 'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
%  xlabel('Number of Iteration','Fontsize',10);
%  h = legend('$\|\bar K_{1k}-\bar K_1^*\|$'); set(h,'Interpreter','latex')
 
 subplot(2,2,2);plot(1:length(errorP2),errorP2,'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
 xlabel('Number of Iteration','Fontsize',10);
 h = legend('$\|P_{2}^{(k)}-P_2^*\|$'); set(h,'Interpreter','latex')
 
 subplot(2,2,3);plot(1:length(errorP3),errorP3,'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
 xlabel('Number of Iteration','Fontsize',10);
 h = legend('$\|P_{3}^{(k)}-P_3^*\|$'); set(h,'Interpreter','latex')
 
 subplot(2,2,4);plot(1:length(errorP4),errorP4,'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
 xlabel('Number of Iteration','Fontsize',10);
 h = legend('$\|P_{4}^{(k)}-P_4^*\|$'); set(h,'Interpreter','latex')
 
end

function dX=mysys(t,X)
    
global bA bB bD bD2 bK0 E Ke H Hrho mu
    
xi1=X(1:4); xi2=X(5:8); xi3=X(9:12); xi4=X(13:16);
zeta1=X(17:18);zeta2=X(19:20);zeta3=X(21:22);zeta4=X(23:24);
xi=[xi1;xi2;xi3;xi4];
zeta=[zeta1;zeta2;zeta3;zeta4];
v=X(25:26); vstar=X(27:28);
psi1=[v;zeta1]; psi2=[v;zeta2]; psi3=[v;zeta3]; psi4=[v;zeta4];

noise =0.3*[sum(sin([38.1558   76.5517   79.5200   18.6873]*t));
            sum(sin([17.9703   15.5098  -33.7388  -38.1002   100   45.9744  -15.9614    8.5268  -27.6188   25.1267]*t));
            sum(sin([7.9703   5.5098  -3.7388  -8.1002   10   5.9744  -5.9614    .5268  -7.6188   5.1267]*t));
            sum(sin([137.9703   125.5098  -323.7388  -8.1002   1000   425.9744  -.9614    .5268  -.6188   .1267]*t))];
%tildew = 0; %
tildew=1*sin(3*t);
u=-bK0*xi + noise;
%dxi=bA*xi+bB*u+[bD bD2]*[v;v;v;v;zeta];
dxi1=bA(1:4,1:4)*xi1+bB(1:4,1)*u(1)+bD(1:4,1:2)*v+bD2(1:4,1:2)*zeta1;
dxi2=bA(5:8,5:8)*xi2+bB(5:8,2)*u(2)+bD(5:8,3:4)*v+bD2(5:8,3:4)*zeta2;
dxi3=bA(9:12,9:12)*xi3+bB(9:12,3)*u(3)+bD(9:12,5:6)*v+bD2(9:12,5:6)*zeta3;
dxi4=bA(13:16,13:16)*xi4+bB(13:16,4)*u(4)+bD(13:16,7:8)*v+bD2(13:16,7:8)*zeta4;

dxi=[dxi1;dxi2;dxi3;dxi4];

dv=E*v+H*(-Ke*(v-vstar)+tildew);
dvstar=E*vstar;
dzeta=(kron(eye(4),E)-mu*kron(Hrho,eye(2)))*zeta+mu*kron(Hrho,eye(2))*kron([1;1;1;1],v);


dxi1xi1=kron(xi1',xi1')'; % 16
dxi1u1=kron(xi1',u(1))'; %4 
dxi1v=kron(xi1',psi1')'; % 16

dxi2xi2=kron(xi2',xi2')'; % 16
dxi2u2=kron(xi2',u(2))'; %4 
dxi2v=kron(xi2',psi2')'; % 16

dxi3xi3=kron(xi3',xi3')'; % 16
dxi3u3=kron(xi3',u(3))'; %4 
dxi3v=kron(xi3',psi3')'; % 16

dxi4xi4=kron(xi4',xi4')'; % 16
dxi4u4=kron(xi4',u(4))'; %4 
dxi4v=kron(xi4',psi4')'; % 16

dX=[dxi; % 16   
    dzeta; % 8
    dv; % 2
    dvstar; % 2  % 16+8+4=28      
        dxi1xi1; % 16
        dxi1u1; % 4
        dxi1v; % 16 % 16+4+16=36
        dxi2xi2; % 16
        dxi2u2; %4 
        dxi2v; % 16 16+4+16=36
        dxi3xi3; % 16
        dxi3u3; %4 
        dxi3v; % 8 %36
        dxi4xi4; % 16
        dxi4u4; %4 
        dxi4v; % 16 % 36
        ];
end



function dX=mysys2(t,X)
    
global bA bB bD bD2 bK0 E Ke H Hrho mu BKlearn
    
xi1=X(1:4); xi2=X(5:8); xi3=X(9:12); xi4=X(13:16);
zeta1=X(17:18);zeta2=X(19:20);zeta3=X(21:22);zeta4=X(23:24);
xi=[xi1;xi2;xi3;xi4];
zeta=[zeta1;zeta2;zeta3;zeta4];
v=X(25:26); vstar=X(27:28);
psi1=[v;zeta1]; psi2=[v;zeta2]; psi3=[v;zeta3]; psi4=[v;zeta4];


tildew=sin(3*t);
u=-BKlearn*xi;
dxi1=bA(1:4,1:4)*xi1+bB(1:4,1)*u(1)+bD(1:4,1:2)*v+bD2(1:4,1:2)*zeta1;
dxi2=bA(5:8,5:8)*xi2+bB(5:8,2)*u(2)+bD(5:8,3:4)*v+bD2(5:8,3:4)*zeta2;
dxi3=bA(9:12,9:12)*xi3+bB(9:12,3)*u(3)+bD(9:12,5:6)*v+bD2(9:12,5:6)*zeta3;
dxi4=bA(13:16,13:16)*xi4+bB(13:16,4)*u(4)+bD(13:16,7:8)*v+bD2(13:16,7:8)*zeta4;

dxi=[dxi1;dxi2;dxi3;dxi4];
dv=E*v+H*(-Ke*(v-vstar)+tildew);
dvstar=E*vstar;
dzeta=(kron(eye(4),E)-mu*kron(Hrho,eye(2)))*zeta+mu*kron(Hrho,eye(2))*kron([1;1;1;1],v);


dX=[dxi; % 16   
    dzeta; % 8
    dv; % 2
    dvstar; % 2  % 16+8+4=28      
        ];
end



