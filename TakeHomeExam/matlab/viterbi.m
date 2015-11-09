function q=viterbi(K,Pi,U,A,means,Sigmas)
    [T,p]=size(U);
    q=zeros(T,1);
    NormDist=zeros(K,T);
    for i=1:K
        for t=1:T
            normconstant=(2*pi)^(-p/2)/sqrt(det(Sigmas{i}));
            invers=inv(Sigmas{i});
            d=U(t,:)-means(i,:);
            NormDist(i,t)=normconstant*exp(-0.5*d*invers*d');
        end
    end
    
    E=log(NormDist);
    P=log(A);
    
    bp = zeros(K,T);
    bp(:,1) = 1:K;
    vit = log(Pi(:))+E(:,1);
    
    for t = 2:T
        [vit,idx] = max(bsxfun(@plus,vit,P),[],1);
        vit = vit(:)+E(:,t);
        bp = bp(idx,:);
        bp(:,t) = 1:K;
    end
    
    [vit,idx] = max(vit);
    
    q=bp(idx,:);
end