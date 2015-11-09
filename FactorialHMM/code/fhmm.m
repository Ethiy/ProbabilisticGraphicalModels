function [W,C,P,Pi,LogLik]=fhmm(X,K,M,Cycles)
    %% Initialisation
    [T, p]=size(X);
    C=diag(rand(1,p));
    W=randn(K*M,p);
    
    Pi=rand(K,M);
    Pi=columnDivision(Pi,columnSum(Pi));
    
    P=rand(K*M,K);
    P=rowDivision(P,rowSum(P));
    
    LogLik=[];
    
    Ranker=zeros(K^M,M);
    for k=1:K^M
        Ranker(k,:)=base(k-1,K,M);
    end
    
    alpha=zeros(T,K^M);
    B=zeros(T,K^M);  %emission proba      
    beta=zeros(T,K^M);
    gamma=zeros(T,K^M);    
    eta =zeros(T*K*M,K*M);
    GammaX=zeros(K*M,p);
    Eta=zeros(K*M,K*M); 
    xi=zeros(M*K,K);
    
    
    %% Rearangement  variables for the FB
    Wmk=zeros(K^M,p);
    Pm=ones(K^M,K^M);
    Pis=ones(K^M,1);
    collapse=zeros(K^M,M*K);
    collapseb=zeros(K^M,M*K*M*K);
    for i=K^M
            for m=1:M
                collapse(i,(m-1)*K+Ranker(i,m))=1;
                for l=1:M
                    collapseb(i,((m-1)*K+Ranker(i,m)-1)*M*K+(l-1)*K+Ranker(i,l))=1;
                end
            end
    end
    
    
    
    normconstant=(2*pi)^(-p/2);
    
    for c=1:Cycles
        %% Forward Backward
        
        %% rearanging
        for i=K^M
            for m=1:M
                Wmk(i,:)=Wmk(i,:)+W((m-1)*K+Ranker(i,m),:);
                Pis(i,:)=Pis(i,:)*Pi(Ranker(i,m),m);
            end
        
            for j=1:K^M
                for l=1:M
                    Pm(i,j)=Pm(i,j)*P((l-1)*K+Ranker(i,l),Ranker(j,l));
                end
            end
        end
        
        %% emission probabilities
        
        invers=inv(C);      
        cst=normconstant/sqrt(det(C));
        
        for l=1:(K^M),
            d=ones(T,1)*Wmk(l,:)-X(1:T,:);
            B(:,l)=cst*exp(-0.5*rowSum((d*invers).*d));
        end
        
        %% Forward
        
        scale=zeros(T,1);
        alpha(1,:)=Pis'.*B(1,:);
        scale(1)=sum(alpha(1,:)); 
        alpha(1,:)=alpha(1,:)/scale(1);
        for i=2:T
            alpha(i,:)=(alpha(i-1,:)*Pm).*B(i,:); 
            scale(i)=sum(alpha(i,:));
            alpha(i,:)=alpha(i,:)/scale(i);
        end
        
        %% Backward
        
        beta(T,:)=ones(1,K^M)/scale(T);
        for i=T-1:-1:1
            beta(i,:)=(beta(i+1,:).*B(i+1,:))*(Pm')/scale(i); 
        end
        
        %% Inference
        
        gamma=(alpha.*beta); 
        gamma=rowDivision(gamma,rowSum(gamma));
                       
        gamma1=gamma*collapse;
        for i=1:T
            for j=1:M
                ind1=(j-1)*K+1:j*K;
                gamma1(i,ind1)=gamma1(i,ind1)/sum(gamma1(i,ind1));
            end
        end
        
        for i=1:T-1
            tmp=(alpha(i,:)*collapse)'*((beta(i+1,:).*B(i+1,:))*collapse);
            for j=1:M
                ind1=(j-1)*K+1:j*K;
                tmp2=P(ind1,:).*tmp(ind1,ind1);
                xi(ind1,:)=xi(ind1,:)+tmp2/sum(tmp2(:));
            end    
        end
        
        tmp=gamma*collapseb;
        
        for i=1:T
            ind1=(i-1)*K*M+1:i*K*M;
            eta(ind1,:)=reshape(tmp(i,:),[],K*M);
            for j=1:M
                ind2=(i-1)*K*M+(j-1)*K+1:(i-1)*K*M+j*K;
                for l=1:M
                    if (j==l)
                        eta(ind2,(j-1)*K+1:j*K)=diag(gamma1(i,(j-1)*K+1:j*K)); 
                    else
                        tmp3=sum(sum(eta(ind2,(l-1)*K+1:l*K)));
                        eta(ind2,(l-1)*K+1:l*K)=eta(ind2,(l-1)*K+1:l*K)/tmp3;
                    end
                end
            end
            Eta=Eta+eta(ind1,:);
            GammaX=GammaX+gamma1(i,:)'*X(i,:);
        end
        
        Scale=log(scale);
        Eta=(Eta+Eta')/2;
            
        LogLik=[LogLik sum(Scale)];
        
        
        %% M Step
        
        W=pinv(Eta)*GammaX;
        
        C=X'*X/T-GammaX'*W/T;
        
        for i=1:K*M
            ind1=sum(xi(i,:));
            if(ind1==0)
                P(i,:)=ones(1,K)/K;
            else
                P(i,:)=xi(i,:)/ind1;
            end
        end
        
        Pi=reshape(gamma1(1,:),K,M);
    end


end