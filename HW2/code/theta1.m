function [Pi , Mu , Sigma]=theta1(X,tau)
    Pi=[mean(tau(:,1)),mean(tau(:,2)),mean(tau(:,3)),mean(tau(:,4))];
    
    Mu=zeros(4,2);
    for j=1:4
        Mu(j,:)=sum([tau(:,j).*X(:,1),tau(:,j).*X(:,2)])./sum(tau(:,j));
    end
        
    sigma2=zeros(4,1);
    n=length(X);
    for j=1:4
        for i=1:n
            sigma2(j,1)=sigma2(j,1)+tau(i,j)*((X(i,:)-Mu(j,:))*((X(i,:)-Mu(j,:))'));
        end
        sigma2(j)=sigma2(j)/(2*sum(tau(:,j)));
    end
    Sigma=zeros(2,2,4);
    for k=1:4
        Sigma(:,:,k)=sigma2(k)*eye(2,2);
    end
end