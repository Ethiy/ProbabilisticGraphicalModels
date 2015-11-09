function y=densQ(x,pi,mu1,mu0,sigma1,sigma2)
    x1=[x(:,1)-mu1(1),x(:,2)-mu1(2)];
    x0=[x(:,1)-mu0(1),x(:,2)-mu0(2)];
    y1=x1*sigma1^-1*x1';
    y2=x0*sigma2^-1*x0';
    z1=zeros(length(x),1);
    z2=zeros(length(x),1);
    for i=1:length(x)
      z1(i)=y1(i,i);
      z2(i)=y2(i,i);
    end
    y=(pi*exp(-1/2*z1))./(pi*exp(-1/2*z1)+(1-pi)*exp(-1/2*z2));
end