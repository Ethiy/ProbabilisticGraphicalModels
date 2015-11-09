function y=normal(x,mu,Sigma)
    y=1/(2*pi)*(1/det(Sigma))^(1/2)*exp(-1/2*(x-mu)*Sigma^-1*(x-mu)');
end