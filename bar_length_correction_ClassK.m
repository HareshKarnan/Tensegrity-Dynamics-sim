function [N_Cor,Ndot_Cor]=bar_length_correction_ClassK(N,Ndot,C_b,P,D,bar_len_const_hat)

C_r = 1/2*abs(C_b);
q=0.5;

B=N*C_b';Bdot=Ndot*C_b';
R=N*C_r';Rdot=Ndot*C_r';

for j=1:size(B,2)
    
    a0=(q*bar_len_const_hat(j,j))^2*norm(Bdot(:,j))^2*((B(:,j)'*Bdot(:,j))^2-norm(B(:,j))^2*norm(Bdot(:,j))^2);
    a1=2*(q*bar_len_const_hat(j,j))^2*((B(:,j)'*Bdot(:,j))^2-norm(B(:,j))^2*norm(Bdot(:,j))^2);
    a2=norm(Bdot(:,j))^4-(q*bar_len_const_hat(j,j))^2*norm(B(:,j))^2;
    a3=2*norm(Bdot(:,j))^2;
    x=roots([1 a3 a2 a1 a0]);
    x= x(imag(x)==0);
    J=zeros(1,size(x,1));
    for i=1:size(x,1)
        v=q*bar_len_const_hat(j,j)*((x(i)*eye(3)+Bdot(:,j)*Bdot(:,j)')\B(:,j));
        p=bar_len_const_hat(j,j)*v-B(:,j);
        r=-v*v'*Bdot(:,j);
        J(i)=q*norm(p)^2+norm(r)^2;
    end
    [~,min_index]=min(J);
    v=q*bar_len_const_hat(j,j)*((x(min_index)*eye(3)+Bdot(:,j)*Bdot(:,j)')\B(:,j));
    p=bar_len_const_hat(j,j)*v-B(:,j);
    r=-v*v'*Bdot(:,j);
    B(:,j)=B(:,j)+p;
    Bdot(:,j)=Bdot(:,j)+r;
end

    [U,SIGMA,V] = svd(C_r*P);
    rank_SIGMA = rank(SIGMA);
    U1 = U(:,1:rank_SIGMA);
    U2 = U(:,rank_SIGMA+1:end);
    V1 = V(:,1:rank_SIGMA);
    V2 = V(:,rank_SIGMA+1:end);
    Sigma1 = SIGMA(1:rank_SIGMA,1:rank_SIGMA);
    R = (.5*D-.25*B*C_b*P)*V1*(Sigma1^-1)*U1'+R*(U2*U2');
    Rdot = -.25*Bdot*C_b*P*V1*(Sigma1^-1)*U1'+Rdot*(U2*U2');
    N_Cor = [B R]/[C_b' C_r'];
    Ndot_Cor = [Bdot Rdot]/[C_b' C_r'];

if norm(C_b(:,size(C_b,2)))==0
    for k=1:size(C_b,2)-2*size(B,2)
        N_Cor(:,2*size(B,2)+k)=N(:,2*size(B,2)+k);
        Ndot_Cor(:,2*size(B,2)+k)=Ndot(:,2*size(B,2)+k);
    end
end