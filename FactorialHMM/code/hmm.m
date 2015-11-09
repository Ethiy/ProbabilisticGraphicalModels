function [means,Sigmas,A,Pi,LogLik]=hmm(K,U,Cycles)
    %% Initialization
       
    [T,p]=size(U);
    
    C=diag(rand(1,p));
    Sigmas={C};
    for k=2:K
        Sigmas{k}=C;
    end
    
    means=rand(K,p);
    
    Pi=rand(1,K);
    Pi=Pi/sum(Pi);

    A=rand(K);
    A=rowDivision(A,rowSum(A));

    LogLik=[];
    
    for c=1:Cycles
    %% Forward Backward
    
        [Alpha Beta scale NormDist]=forwardBackward(U,K,A,Pi,means,Sigmas);
    
        Gamma=(Alpha.*Beta); 
        Gamma=rowDivision(Gamma,rowSum(Gamma));
        Gammasum=sum(Gamma);
        
        Xi=zeros(T-1,K*K);
        for i=1:T-1
            t=A.*( Alpha(i,:)' * (Beta(i+1,:).*NormDist(i+1,:)));
            Xi(i,:)=t(:)'/sum(t(:));
        end
        
    %% M Step
    
        means=zeros(K,p);
        means=Gamma'*U;
        means=rowDivision(means,Gammasum');
        
        sxi=rowSum(Xi')';
        sxi=reshape(sxi,K,K);
        A=rowDivision(sxi,rowSum(sxi));
        
        Pi=Gamma(1,:);
        
        for i=1:K
            d=(U-ones(T,1)*means(i,:));
            Sigmas{i}=rowProduct(d,Gamma(:,i))'*d;
            Sigmas{i}=Sigmas{i}/sum(Gamma(:,i));
        end
        
    %% Log Likelihood
        scale=log(scale);
        LogLik=[LogLik sum(scale)];
    end

end