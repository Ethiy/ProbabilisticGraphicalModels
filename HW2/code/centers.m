function mu1=centers(X,mu0)
    mu1=[];
    for i=1:4
        M=cluster(i,X,mu0);
        mu1=cat(1,mu1,[mean(M(:,1)),mean(M(:,2))]);
    end
end