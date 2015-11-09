function z=funky1(x,y,Sigma)
    z=Sigma(1,1,1)*(x-Mu(1,1))^2+Sigma(2,2,1)*(y-Mu(1,2))^2+2*Sigma(2,1,1)*(x-Mu(1,1))*(y-Mu(1,2))-3;
end