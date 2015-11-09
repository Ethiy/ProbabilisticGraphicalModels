function Y=cluster(j,X,mu)
    n=length(X);
    D=distance(X,mu);
    Y=[];
    for i=1:n
        if D(i,j)==min(D(i,:))
            Y=cat(1,Y,X(i,:));
        end
    end
end