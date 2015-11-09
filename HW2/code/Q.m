function y=Q(X,Pi , Mu , Sigma)
    n=length(X);
    y=zeros(n,4);
    for i=1:n
        for j=1:4
            y(i,j)=Pi(j)*normal(X(i,:),Mu(j,:),Sigma(:,:,j));
        end
        y(i,:)=y(i,:)/sum(y(i,:));
    end
end