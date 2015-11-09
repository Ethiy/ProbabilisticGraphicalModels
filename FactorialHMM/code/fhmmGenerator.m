function X=fhmmGenerator(T,N)

    Cov=.0025; 
    M=3;
    K=2;
    
    N=T*N;
    X=zeros(N,1);
    S=zeros(N,3);

    W=rand(K,M);

    P1=rand(K,K);
    P1=rowDivision(P1,rowSum(P1));   
    
    P2=rand(K,K);
    P2=rowDivision(P2,rowSum(P2)); 
    
    P3=rand(K,K);
    P3=rowDivision(P3,rowSum(P3));
    
    P={P1,P2,P3};
    
    u=rand(1,3);
    
    for m=1:3
        if u(m)<.5
            S(1,m)=1;
        else
            S(1,m)=2;
        end
    end
    
    state=eye(2);
    
    mu=0;
    for m=1:3
        mu=W(:,m)'*state(:,S(1,m));
    end
    
    X(1)=mu+Cov*randn;
    
    for n=2:N
        u=rand(1,3);
        for m=1:3
            if u(m)<P{m}(S(n-1,m),1)
                S(n,m)=1;
            else
                S(n,m)=2;
            end
        end
        mu=0;
        for m=1:3
            mu=W(:,m)'*state(:,S(n,m));
        end
        X(n)=mu+Cov*randn;
    end
    
end