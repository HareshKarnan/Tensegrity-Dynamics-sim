function dy = classkcontrol_reducedorder(t,x0,Inp)

Cs = Inp.Cs;Cb = Inp.Cb;W = Inp.W;s0 = Inp.s0;b0 = Inp.b0;P=Inp.P;D=Inp.D;Yt=Inp.Yt;L=Inp.L;R=Inp.R;a=Inp.a;b=Inp.b;
k = Inp.k;damp=Inp.damp;

persistent gammashist
persistent strlength
persistent Omghist
strlength(:,1) = diag(s0);
gammashist(:,1)=zeros(size(Cs,1),1);

X2 = x0(1:numel(x0)/2);
X2d = x0(numel(x0)/2+1:end);

[U,Ez,V] = svd(P);
rankEz = rank(Ez);
U1 = U(:,1:rankEz);
U2 = U(:,rankEz+1:end);
E0 = Ez(1:rankEz,1:rankEz);

n = size(Cb,2); % number of nodes
nb = size(Cb,1); % number of bars
ns = size(Cs,1); % number of strings
X2 = reshape(X2,3,size(U2,2));
X2d = reshape(X2d,3,size(U2,2));
l_hat = b0;
X1 = D*V*E0^-1;
N = [X1 X2]*U';
Nd = [0*X1 X2d]*U';

m=1; % mass of each bar = 1 kg
%%  bar length correction
[N,Nd]=bar_length_correction_ClassK(N,Nd,Cb,P,D,l_hat);
X2d = Nd*U2;
X2 = N*U2;

%% calculate lagrange multiplier
B=N*Cb';Bd=Nd*Cb';S=N*Cs';Sd=Nd*Cs';
Cr=(1/2)*abs(Cb);
m_hat = m.*eye(nb);
M = (1/12)*Cb'*m_hat*Cb + Cr'*m_hat*Cr;
Minv = 3*Cb'*m_hat^-1*Cb + 4*Cr'*m_hat^-1*Cr;

% gamma_hat = diag(gammashist(:,end));
gamma_hat = find_gamma(strlength(:,end),S,Sd,damp,k,ns);

if t==0
% gamma_hat = zeros(ns,ns);
% Omghist(:,:,1) = zeros(3,size(P,2));
end


%% control to find gamma
Omg = calclagmat(N,Cb,B,Bd,Cs,gamma_hat,Minv,P,W,l_hat,m_hat);

% Omg = Omghist(:,:,end);
EYE = eye(size(Cb,1));
biglambda = []; tau = [];
for i=1:size(Cb,1)
    biglambda = [biglambda;(-1/(2*l_hat(i,i)^2))*B(:,i)'*S*diag(Cs*Cb'*EYE(:,i))];
    tau = [tau;(1/(2*l_hat(i,i)^2))*B(:,i)'*(W+Omg*P')*Cb'*EYE(:,i)+ (1/(12*l_hat(i,i)^2))*m*norm(Bd(:,i))^2];
end
Acon = a*eye(size(L,1));
Bcon = b*eye(size(L,1));

Mi = U2*(U2'*M*U2)^-1*U2';

BTu = L*W*Mi*R+Acon*L*X2d*U2'*R+Bcon*(L*[X1 X2]*U'*R-Yt);
BIGTAU=[];meu=[];
EYE = eye(size(R,2));

for i=1:size(R,2)
BIGTAU = [BIGTAU;L*(S*diag(Cs*Mi*R*EYE(:,i))+B*diag(Cb*Mi*R*EYE(:,i))*biglambda)];
meu = [meu;BTu*EYE(:,i) - L*B*diag(Cb*Mi*R*EYE(:,i))*tau];
end

options = optimoptions('lsqlin','Algorithm','interior-point','Display','off');
gammas = lsqlin(BIGTAU,meu,[],[],[],[],zeros(ns,1),[],[],options);
% gammas = pinv(BIGTAU)*meu;
% [gammas,fval,exitflag] = linprog(ones(ns,1),[],[],BIGTAU,meu,zeros(ns,1),[]);
% if exitflag<=0
% gammas = gammashist(:,end);
% end
gammashist = [gammashist gammas];

% find restlength of the strings
sinit=[];
for i=1:ns
    sinit(i,1) = norm(S(:,i)) - (gammas(i)*norm(S(:,i))-norm(Sd(:,i))*damp)*k^-1;
end
% plot(t,norm(S(:,1))-sinit(1),'-o','Color','b');hold on;
% plot(t,norm(S(:,4))-sinit(4),'-o','Color','r');hold on;
% title('sinit');
strlength = [strlength sinit]; 

%% convert to reduce order dynamics and integrate
gamma_hat = diag(gammas);
% Omg = calclagmat(N,Cb,B,Bd,Cs,gamma_hat,Minv,P,W,l_hat,m_hat);
% Omghist(:,:,end+1) = Omg

% lambda_c = diag(-biglambda*gammas-tau);
lambda_c = diag(diag(0.5*l_hat^-2*B'*(S*gamma_hat*Cs-W-Omg*P')*Cb'-(1/12)*l_hat^-2*(m*eye(nb))*Bd'*Bd));
K = Cs'*gamma_hat*Cs - Cb'*lambda_c*Cb;
% K = pinv(N)*(W*M^-1*U1 + Omg*P'*M^-1*U1)*pinv(U1);
M_til = U2'*M*U2;
K_til = U2'*K*U2;
W_til = W*U2-X1*U1'*K*U2;
X2dd = (W_til-X2*K_til)*pinv(M_til);
dy = [reshape(X2d,numel(X2d),1);reshape(X2dd,numel(X2dd),1)];

end


function gammas_hat = find_gamma(strlength,S,Sd,damp,k,ns)
gammas_hat = diag(zeros(ns,1));
for i=1:ns
    
    if(norm(S(:,i))>strlength(i))
        gammas_hat(i,i) = k*(1-strlength(i)/norm(S(:,i)))+damp*(S(:,i)'*Sd(:,i))/norm(S(:,i))^2;
    else
        gammas_hat(i,i) = 0; 
    end
end
end