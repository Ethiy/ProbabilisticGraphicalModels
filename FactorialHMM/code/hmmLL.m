function LogLik=hmmLL(U,K,means,Sigmas,A,Pi)
    [T,p]=size(U);
    
    LogLik=[];
    
    %% The normal distribution functions
    NormDist=zeros(T,K);
    for i=1:K
        for t=1:T
            normconstant=(2*pi)^(-p/2)/sqrt(det(Sigmas{i}));
            invers=inv(Sigmas{i});
            d=U(t,:)-means(i,:);
            NormDist(t,i)=normconstant*exp(-0.5*d*invers*d');
        end
    end
    
    %% Initilization of Alpha
    
    % the scale is introduced so as to normalize the Alphas and Betas
    % because otherwise they will go to zero
    scale=zeros(T,1);
    
    Alpha=zeros(T,K);
    
    
    Alpha(1,:)=Pi.*NormDist(1,:);
    scale(1)=sum(Alpha(1,:));
    Alpha(1,:)=Alpha(1,:)/scale(1);
    
    %% Forward
    
    for t=2:T
      Alpha(t,:)=(Alpha(t-1,:)*A).*NormDist(t,:);
      scale(t)=sum(Alpha(t,:));
      Alpha(t,:)=Alpha(t,:)/scale(t);
    end
    
    %% Log likekihood
    
    scale=log(scale);
    LogLik=cumsum(scale);
end