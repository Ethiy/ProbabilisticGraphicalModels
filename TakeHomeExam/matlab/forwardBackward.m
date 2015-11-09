function [Alpha Beta scale NormDist]=forwardBackward(U,K,A,Pi,means,Sigmas)
%K number of states
%A is the Kernel 
%Pi is the a priori distribution
%means=(means(1),...,means(4))
%Sigmas=(Sigmas(1),...,Sigma(4))
    [T,p]=size(U);
    
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
    % because otherwise it will go to zero
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
    
    %% Initialization of Beta
    
    Beta=zeros(T,K);
    Beta(T,:)=ones(1,K)/scale(T);
    
    %% Backward
    
    for t=T-1:-1:1
      Beta(t,:)=(Beta(t+1,:).*NormDist(t+1,:))*(A')/scale(t); 
    end
end