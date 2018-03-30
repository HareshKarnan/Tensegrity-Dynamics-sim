% class k prism reduced order control
function H = class_k_prism_control()
clear all;clc;close all;warning off
% class k tensegrity dynamics with a prism example. 
[N,Cb,Cs] = tenseg_prism_sidestring(7.05,6.18);
% f1 = figure

% tenseg_plot(N,Cb,Cs,f1,[1 2 3 4 5 6])
%%
B = N*Cb';S=N*Cs';b0=sqrt(diag(diag(B'*B)));s0=sqrt(diag(diag(S'*S)));

% fix the bottom 3 nodes to ground
P = [1 0 0;
     0 1 0;
     0 0 1;
     0 0 0;
     0 0 0;
     0 0 0];
D = [N(:,1) N(:,2) N(:,3)];

% external force on the structure = no external forces ! mwahaha. 
W = zeros(3,6);
L = [1 0 0;
     0 1 0;
     0 0 1];
R = [0 0 0;
     0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1];
% optimization program to find the target node
th =0;phi=75;
v = [cosd(phi)*cosd(th);cosd(phi)*sind(th);sind(phi)];
[Yt,sucflag] = getoptimized(th,phi,1);
if(sucflag<0)
    fprintf('Control cannot proceed. Illegal target');
    stop
end
Ntar = [N(:,1:3) Yt];

% Ntar = [ -3.09             0          3.09          4.58         -1.66         -1.68;
%          -1.78          3.57         -1.78          1.04         -2.41          4.83;
%              0             0             0          5.52          7.66          5.54];
Yt = Ntar(:,4:6);
Btar = Ntar*Cb';
b0tar=sqrt(diag(diag(B'*B)));
% tenseg_plot([N(:,1:3) Yt],Cb,Cs,f1,[1 2 3 4 5 6])
% string stiffness
k = 100;
damp = 0.1;
a=2;
b=3;


[U,Ez,V] = svd(P);
U1 = U(:,1:size(P,2));U2 = U(:,size(P,2)+1:end); E0 = Ez(1:size(P,2),:);
Nd = 0*N;
x0 = [N*U2 , Nd*U2];
dt = 0.1; tf = 20;
t = 0:dt:tf;


Inp.a=a;Inp.b=b;Inp.Cs = Cs; Inp.Cb=Cb; Inp.W=W; Inp.s0=s0; Inp.b0=b0;Inp.tf=tf;Inp.P=P;Inp.D=D;Inp.B=B;Inp.S=S;Inp.Yt=Yt;Inp.L=L;Inp.R=R;
Inp.k=k; Inp.damp=damp;
y = ode4plus_barlength(@classkcontrol,t,x0,Inp);
X1 = D*V*E0^-1;

%% figure
E=[];ble=[];angle_error=[];
for i = 1:size(y,1)
    X2 = (y(i,1:size(y,2)/2));
    X2 = reshape(X2,3,size(U2,2));
    N = [X1 X2]*U^-1;
    Ntrace(:,:,i) = N;
    Btrace(:,:,i) = N*Cb';
    Strace(:,:,i) = N*Cs';
    RBtrace(:,:,i) = 0.5*N*abs(Cb');
    RStrace(:,:,i) = 0.5*N*abs(Cs');
    ble = [ble sqrt(diag(B'*B))-sqrt(diag(Btrace(:,:,i)'*Btrace(:,:,i)))];
    mnb=N(:,4:6)-Yt;
    E = [E;mnb(1,:) mnb(2,:) mnb(3,:) ];
    topstr = Strace(:,1:3,i);
    angle_error(i,:) = [dot(topstr(:,1),v) dot(topstr(:,2),v) dot(topstr(:,3),v)];
end
hold off;
%
figure
for k=1:9
plot(t,E(:,k));hold on;
end
title('Errors in nodes');
legend('x1','x2','x3','y1','y2','y3','z1','z2','z3')
% figure
% for k2=1:3
% plot(t,ble(k2,:)); hold on;
% end
% title('Errors in bl');
figure
for k=1:3
    plot(t,angle_error(:,k)); hold on;
end
title('Angle error');
H.N = Ntrace;
H.Nhist = Ntrace;
H.t=t;
tensegst.C_b = Cb;
tensegst.C_s = Cs;
tenseg_animation(H,tensegst,'hello.avi')
H.B = Btrace;
H.S = Strace;
H.RB = RBtrace;
H.RS = RStrace;
H.tf = tf;
H.dt = t;
hold off;
end

