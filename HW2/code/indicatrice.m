function Y=indicatrice(X,mu)
    n=length(X);
    D=distance(X,mu);
    Y=zeros(n,4);
    for i=1:n
        for j=1:4
            if D(i,j)==min(D(i,:))
                Y(i,j)=1;
            end
        end
    end
end