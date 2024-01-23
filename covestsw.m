function [Rhat sigma] =covest(X,K,indim,sign)
sigma=ones(sign,1);
Rhat=eye(indim);
flag=true;
counter=1;
while flag==true
counter=counter+1;
temp=0;
s=[ones(1,K/2)*sigma(1) ones(1,K/2)*sigma(2)];
for k=1:K
    temp=temp+(1/indim)*(1/s(k))*X(:,k)*X(:,k)';
end
Rhat=temp;
Rhat=Rhat/norm(Rhat,'fro');
temp=Rhat^(-0.5)*X;
for i=1:sign
    if i==1
    sigma(i)=(2/(indim*K))*trace(temp(:,1:K/2)'*temp(:,1:K/2));
    end
    if i==2
    sigma(i)=(2/(indim*K))*trace(temp(:,K/2+1:end)'*temp(:,K/2+1:end));
    end
end
if counter==10
    flag=false
end
end

end