function [W,C,P,Pi,LogLik]=phmm(X,K,M,Cycles,iter)
    %% Initialisation
    [T, p]=size(X);
    epsilon=exp(-700);
    
    C=diag(rand(1,p));
    W=randn(K*M,p);
    
    Pi=rand(K,M);
    Pi=columnDivision(Pi,columnSum(Pi));
    
    P=rand(K*M,K);
    P=rowDivision(P,rowSum(P));
    
    LogLik=[];
    
    normconstant=(2*pi)^(-p/2);
    
    mS=ones(T,M*K)/K; 
    h=ones(T,M*K)/K;
    logmS=log(mS);
    exph=exp(h);
    
    alpha=zeros(T,M*K);
    beta=zeros(T,M*K);
    gamma=zeros(T,K*M); 
    
    for c=1:Cycles
        
        %% E step
        
        invers=inv(C);      
        cst=normconstant/sqrt(det(C));
        for l=1:iter
            mSold=mS;
            logmSold=logmS;
            
            %% Finding the new h
            
            for i=1:T
               Yhat=mS(i,:)*W;
               for j=1:M
                   ind1=(j-1)*K+1:j*K;
                   Wj=W(ind1,:);
                   mSold=mS(i,ind1);
                   logP=log(P(ind1,:)+epsilon);
                   logPi=log(Pi(:,j)+epsilon);
                   h(i,ind1)=(Wj*invers*(X(i,:)-Yhat)'+Wj*invers*Wj'*mS(i,ind1)'- 0.5*diag(Wj*invers*Wj'))';
                   h(i,ind1)=h(i,ind1)-max(h(i,ind1)')'*ones(1,K);
               end
            end
            exph=exp(h);
            
            
            %% Computing mean S for next time
            
            scale=zeros(T,M);
            for j=1:M
                ind1=(j-1)*K+1:j*K;
                alpha(1,ind1)=exph(1,ind1).*(Pi(:,j)');
                scale(1,j)=rowSum(alpha(1,ind1))+epsilon;
                alpha(1,ind1)=rowDivision(alpha(1,ind1),scale(1,j));
                for i=2:T
                    alpha(i,ind1)=(alpha(i-1,ind1)*P(ind1,:)).*exph(i,ind1);
                    scale(i,j)=rowSum(alpha(i,ind1))+epsilon;
                    alpha(i,ind1)=rowDivision(alpha(i,ind1),scale(i,j));
                end
                beta(T,ind1)=rowDivision(ones(1,K),scale(T,j));
                for i=T-1:-1:1
                    beta(i,ind1)=(beta(i+1,ind1).*exph(i+1,ind1))*(P(ind1,:)');
                    beta(i,ind1)=rowDivision(beta(i,ind1),scale(i,j));
                end
    
                mS(:,ind1)=(alpha(:,ind1).*beta(:,ind1));
                mS(:,ind1)=rowDivision(mS(:,ind1),rowSum(mS(:,ind1))+epsilon);
            
            end
            
            logmS=log(mS+(mS==0).*epsilon);
        end
        
        %%  Computing Likelihood 
        
        Ranker=zeros(K^M,M);
        for k=1:K^M
            Ranker(k,:)=base(k-1,K,M);
        end
        
        MS=ones(T,K^M);
        Wkm=zeros(K^M,p);
        for i=1:K^M
            for j=1:M;
                Wkm(i,:)=Wkm(i,:)+W((j-1)*K+Ranker(i,j),:);
                MS(:,i)=MS(:,i).*mS(:,(j-1)*K+Ranker(i,j));
            end
        end
        
        logPi=log(Pi+(Pi==0)*epsilon);
        logP=log(P+(P==0)*epsilon);
        logmS=log(mS+(mS==0)*epsilon);
        
        ll=0;
        
        for l=1:(K^M)
            d= ones(T,1)*Wkm(l,:)-X;
            ll=ll - 0.5*sum(MS(:,l).*rowSum((d*invers).*d));
        end
        
        ll=ll+T*log(cst);
        
        ll=ll+sum(mS(1,:)*logPi(:))-sum(sum(mS(1,:).*logmS(1,:)));
        
        for i=2:T
            for j=1:M
                ind2=(j-1)*K+1:j*K; 
                ll=ll+sum(sum(mS(i-1,ind2).*(mS(i,ind2)*logP(ind2,:)')))-sum(sum(mS(i,ind2).*logmS(i,ind2)));
            end
        end
        
        LogLik=[LogLik, ll];
        
        %% Inference
        
        gamma=mS;
        Eta=gamma'*gamma;
        gammasum=sum(gamma);
        for j=1:M
            ind2=(j-1)*K+1:j*K;
            Eta(ind2,ind2)=diag(gammasum(ind2)); 
        end
        
        GammaX=gamma'*X;
        Eta=(Eta+Eta')/2;
        
        xi=zeros(M*K,K); 
        for i=1:T-1
            for j=1:M
                ind=(j-1)*K+1:j*K;
                tmp = P(ind,:).*(alpha(i,ind)'*(beta(i+1,ind).*exph(i+1,ind)));
                xi(ind,:)=xi(ind,:)+tmp/sum(tmp(:));
            end
        end
        
        %% M Step
        
         W=pinv(Eta)*GammaX;
         
        C=X'*X/T-GammaX'*W/T;
        
        Pi=reshape(gamma(1,:),K,M);
        
        for i=1:K*M
            ind1=sum(xi(i,:));
            if(ind1==0)
                P(i,:)=ones(1,K)/K;
            else
                P(i,:)=xi(i,:)/ind1;
            end
        end

    end
end