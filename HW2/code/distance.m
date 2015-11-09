function y=distance(X,mu)
    n=length(X);
    Z=zeros(n,4);
    for i=1:n
        for j=1:4
            Z(i,j)=sqrt((X(i,1)-mu(j,1))^2+(X(i,2)-mu(j,2))^2);
        end
    end
    y=Z;
end
