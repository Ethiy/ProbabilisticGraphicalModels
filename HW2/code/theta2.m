function [Pi , Mu , Sigma]=theta2(X,tau)
    Pi=[mean(tau(:,1)),mean(tau(:,2)),mean(tau(:,3)),mean(tau(:,4))];
    
    Mu=zeros(4,2);
    for j=1:4
        Mu(j,:)=sum([tau(:,j).*X(:,1),tau(:,j).*X(:,2)])./sum(tau(:,j));
    end
        
    Sigma=zeros(2,2,4);
    n=length(X);
    for j=1:4
        for i=1:n
            Sigma(:,:,j)=Sigma(:,:,j)+tau(i,j)*((X(i,:)-Mu(j,:))'*((X(i,:)-Mu(j,:))));
        end
         Sigma(:,:,j)= Sigma(:,:,j)/sum(tau(:,j));
    end
end