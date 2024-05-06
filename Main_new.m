clc
clear all
close all


%% initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indim=5;  K=500;  N=2000;  p=2;  t=1;  rho=5;  eta=0.00002;  ASNR=8;  iteration=2000;
sigma1=sqrt(5);  sigma2=sqrt(5);  sigmat=sqrt(5);  alpha=0; epsilonr=0.00;

%%
H=zeros(indim,p);
freq=[0.1 0.15 0.2];
qreq=[0 -0.025 -0.05];
i=sqrt(-1);
for k=1:p
    for j=1:indim
        H(j,k)=1/sqrt(indim)*exp(-1*i*2*pi*freq(k)*(j-1));
    end
end
for k=1:t
    for j=1:indim
        B(j,k)=1/sqrt(indim)*exp(-1*i*2*pi*qreq(k)*(j-1));
    end
end

C=[H B];
%%

covs=zeros(indim,indim);
for i=1:indim
    for j=1:indim
         covs(i,j)=0.0^(abs(i-j));
    end
end
base_covs=covs/norm(covs,'fro')

dev=zeros(indim,indim);
for i=1:indim
    for j=1:indim
         dev(i,j)=(-0.95)^(abs(i-j));
    end
end

base_covt=base_covs+alpha*dev;
base_covt=base_covt/norm(base_covt,'fro')
covt=sigmat^2*base_covt
a=sigmat*random('normal',0,1,indim,N);
b1=sigma1*random('normal',0,1,indim,K/2);
b2=sigma2*random('normal',0,1,indim,K/2);
b=[b1 b2];
noiset=(base_covt)^(0.5)*(a);
noises=(base_covs)^(0.5)*(b);
%%
[Rs sigmas]=covestsw(noises,K,indim,2)

SS=0;
for i=1:K
   SS=SS+noises(:,i)*noises(:,i)';
end
SC=SS/K
norm(base_covt-Rs,'fro')^2
puse=input('Press enter to Proceed:')
%%%%%%%%%%%%%%%%% Making observations%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu1=[zeros(1,N/2) ones(1,N/2)];
TETA1=random('uniform',0.3,0.35,p,1);
scale=sqrt(((10^(ASNR/10))/((H*TETA1)'*(covt)^(-1)*(H*TETA1))));

for k=1:N  
    TETA1=random('uniform',0.3,0.35,p,1);
    PHI1=random('uniform',0.3,0.35,t,1);
    
    X(:,k)=mu1(k)*H*TETA1;
    Y(:,k)=scale*X(:,k)+B*PHI1+noiset(:,k);
end
SNR=10*log10((scale*H*TETA1)'*(covt)^(-1)*(scale*H*TETA1))
%% Diffrent Tests with plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Htilda=(SC^(-0.5))*H;    
PHtilda=((Htilda)*(((Htilda)')*(Htilda))^(-1)*((Htilda)'));
     
Btilda=(SC^(-0.5))*B;    
PBtildaO=eye(indim)-((Btilda)*(((Btilda)')*(Btilda))^(-1)*((Btilda)'));
     
Ctilda=(SC^(-0.5))*C;    
PCtildaO=eye(indim)-((Ctilda)*(((Ctilda)')*(Ctilda))^(-1)*((Ctilda)'));
for i=1:N
clc; i  

Ytilda=(SC^(-0.5))*Y(:,i);
TS1(i)=Ytilda'*PBtildaO*Ytilda;   TS2(i)=Ytilda'*PCtildaO*Ytilda;  
AMF(i)=TS1(i)/TS2(i);


TS1(i)=Ytilda'*PBtildaO*PHtilda*PBtildaO*Ytilda;   TS2(i)=Ytilda'*PBtildaO*Ytilda;  
ASD(i)=TS1(i)/TS2(i);

ybar=(covt^(-0.5))*Y(:,i);
Bbar=(covt^(-0.5))*B;    PBbarO=eye(indim)-((Bbar)*(((Bbar)')*(Bbar))^(-1)*((Bbar)'));
TS1(i)=ybar'*PBbarO*ybar;
Cbar=(covt^(-0.5))*C;    PCbarO=eye(indim)-((Cbar)*(((Cbar)')*(Cbar))^(-1)*((Cbar)'));
TS2(i)=(ybar'*PCbarO*ybar);
optAMF(i)=TS1(i)/TS2(i);

Hbar=(covt^(-0.5))*H;    PHbar=((Hbar)*(((Hbar)')*(Hbar))^(-1)*((Hbar)'));
TS1(i)=ybar'*PBbarO*PHbar*PBbarO*ybar;   TS2(i)=ybar'*PBbarO*ybar;  
optASD(i)=TS1(i)/TS2(i);



count=1;
R1=Rs;
U1=0.1*ones(indim,indim);
lambda1=3.0;
gamma1=3.0;
Z1=R1;
 while count<iteration
        count=count+1;
        ybar=(R1^(-0.5))*Y(:,i);
        Bbar=(R1^(-0.5))*B;
        Cbar=(R1^(-0.5))*C;
        PCbarO=eye(indim)-((Cbar)*(((Cbar)')*(Cbar))^(-1)*((Cbar)'));
        beta=(Cbar'*Cbar)^(-1)*Cbar'*ybar;
        sigma2=real((1/indim)*ybar'*PCbarO*ybar);
        R1=R1-eta*real((R1^(-1))+4*rho*double((norm(R1-Rs,'fro')^2-epsilonr)>0)*max([0,norm(R1-Rs,'fro')^2-epsilonr])^3*(R1-Rs)...
            +4*lambda1*double((norm(R1-Rs,'fro')^2-epsilonr)>0)*max([0,norm(R1-Rs,'fro')^2-epsilonr])*(R1-Rs)...
            +4*rho*double((norm(R1,'fro')^2-1)>0)*max([0,norm(R1,'fro')^2-1])^3*(R1)...
            +4*gamma1*double((norm(R1,'fro')^2-1)>0)*max([0,norm(R1,'fro')^2-1])*(R1)...
            +U1'+rho*(R1-Z1));
        Z1=Z1-eta*real(-(1/sigma2)*(Z1^(-1)*(Y(:,i)-C*beta)*(Y(:,i)-C*beta)'*Z1^(-1))+rho*(Z1-R1)-U1');
        U1=U1+rho*(R1-Z1);
        %U1=max(0.05,U1);
        gamma1=gamma1+rho*max([0,norm(R1,'fro')^2-1])^2;
        gamma1=min(1000,gamma1);
        lambda1=lambda1+rho*max([0,norm(R1-Rs,'fro')^2-epsilonr])^2; 
        lambda1=min(1000,lambda1);
 end

%      norm(R1,'fro')
Cbar=(R1^(-0.5))*C;    
PCbarO=eye(indim)-((Cbar)*(((Cbar)')*(Cbar))^(-1)*((Cbar)'));
TS2(i)=(ybar'*PCbarO*ybar);
     
count=1;
R0=Rs;
U0=ones(indim,indim);
lambda0=3.0;
gamma0=3.0;
Z0=R0;
loss=[];
normf=[];
while count<iteration
        count=count+1;
        ybar=(R0^(-0.5))*Y(:,i);
        Bbar=(R0^(-0.5))*B;
        Cbar=(R0^(-0.5))*C;
        PBbarO=eye(indim)-((Bbar)*(((Bbar)')*(Bbar))^(-1)*((Bbar)'));
        phi=(Bbar'*Bbar)^(-1)*Bbar'*ybar;
        sigma2=real((1/indim)*ybar'*PBbarO*ybar);
        
        loss=[loss, norm(R0-base_covt,'fro')];
        normf=[normf, norm(R0-Rs,'fro')^2];
        R0=R0-eta*real((R0^(-1))+4*rho*double((norm(R0-Rs,'fro')^2-epsilonr)>0)*max([0,norm(R0-Rs,'fro')^2-epsilonr])^3*(R0-Rs)...
            +4*lambda0*double((norm(R0-Rs,'fro')^2-epsilonr)>0)*max([0,norm(R0-Rs,'fro')^2-epsilonr])*(R0-Rs)...
            +4*rho*double((norm(R0,'fro')^2-1)>0)*max([0,norm(R0,'fro')^2-1])^3*(R0)...
            +4*gamma0*double((norm(R0,'fro')^2-1)>0)*max([0,norm(R0,'fro')^2-1])*(R0)...
            +U0'+rho*(R0-Z0));
        Z0=Z0-eta*real(-(1/sigma2)*(Z0^(-1)*(Y(:,i)-B*phi)*(Y(:,i)-B*phi)'*Z0^(-1))+rho*(Z0-R0)-U0');
        U0=U0+rho*(R0-Z0);
        %U0=max(0.05,U0)
        gamma0=gamma0+rho*max([0,norm(R0,'fro')^2-1])^2;
        gamma0=min(1000,gamma0);
        lambda0=lambda0+rho*max([0,norm(R0-Rs,'fro')^2-epsilonr])^2;
        lambda0=min(1000,lambda0);    
 end
 if i==1
    figure
    subplot(2,1,1)
    plot(real(loss))
    xlabel('Iteration')
    ylabel('Loss')
        subplot(2,1,2)
    plot(real(normf))
    xlabel('Iteration')
    ylabel('Norm')
 end
%      norm(R0,'fro')
Bbar=(R0^(-0.5))*B;    
PBbarO=eye(indim)-((Bbar)*(((Bbar)')*(Bbar))^(-1)*((Bbar)'));
TS1(i)=ybar'*PBbarO*ybar;
Proposed(i)=log(det(R0)/det(R1))+(2*indim)*log(TS1(i)/TS2(i));
end
subplot(1,5,1)
histogram(real(optAMF(1:N/2)),20);hold on; histogram(real(optAMF(N/2+1:end)),50)   
subplot(1,5,2)
histogram(real(AMF(1:N/2)),20);hold on; histogram(real(AMF(N/2+1:end)),50)
subplot(1,5,3)
histogram(real(ASD(1:N/2)),20);hold on; histogram(real(ASD(N/2+1:end)),50)  
subplot(1,5,4)
histogram(real(optASD(1:N/2)),20);hold on; histogram(real(optASD(N/2+1:end)),50) 
subplot(1,5,5)
histogram(real(Proposed(1:N/2)),20);hold on; histogram(real(Proposed(N/2+1:end)),50) 


optAMF=real(optAMF);
AMF=real(AMF);
ASD=real(ASD);
optASD=real(optASD);
Proposed=real(Proposed);
%% ROC calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
last=1000;
starting=-1000;
step=(last-starting)/1000000;
for th=starting:step:last    
    pd1(i)=(200*(sum(AMF(N/2+1:N)>th)))/N;
    pf1(i)=(200*(sum(AMF(1:N/2)>th)))/N;  
    
    pd2(i)=(200*(sum(optAMF(N/2+1:N)>th)))/N;
    pf2(i)=(200*(sum(optAMF(1:N/2)>th)))/N; 
    
    pd3(i)=(200*(sum(ASD(N/2+1:N)>th)))/N;
    pf3(i)=(200*(sum(ASD(1:N/2)>th)))/N;  
    
    pd4(i)=(200*(sum(optASD(N/2+1:N)>th)))/N;
    pf4(i)=(200*(sum(optASD(1:N/2)>th)))/N;  
    
    pd5(i)=(100*(sum(Proposed(N/2+1:N)>th)))/sum(~isnan(Proposed(N/2+1:N)));
    pf5(i)=(100*(sum(Proposed(1:N/2)>th)))/sum(~isnan(Proposed(1:N/2)));   
    
    i=i+1;
end
%% ROC plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(pf1,pd1,'b--','LineWidth',2)
hold on
plot(pf2,pd2,'b','LineWidth',2)
hold on
plot(pf3,pd3,'r--','LineWidth',2)
hold on
plot(pf4,pd4,'r','LineWidth',2)
hold on
plot(pf5,pd5,'k--','LineWidth',2)
grid on
    legend({'AMF','AMF, known R','ASD','ASD, known R','Proposed'}, ...
        'Interpreter', 'LaTeX')
    xlabel('Probability of False Alarm (\%)', 'Interpreter', 'LaTeX')
    ylabel('Probability of Detection (\%)', 'Interpreter', 'LaTeX')
    title('', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')
%  print('Performance-het','-depsc')
SNR=10*log10((scale*H*TETA1)'*(covt)^(-1)*(scale*H*TETA1))
     