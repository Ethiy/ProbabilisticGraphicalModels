clear;

%%---------------------Read the data---------------------------------------
Z=dlmread('classification_data_HWK1/classificationC.train');

X=Z(:,1:2);
Y=Z(:,3);
%We keep the same notations for X and Y

n=length(Y);% is the number of couples (X_i,Y_i)

%%--------------------------Estimators-------------------------------------

Pi=mean(Y);%pi estimator

mu1=zeros(1,2);
mu11=sum(X(:,1).*Y)/sum(Y);
mu12=sum(X(:,2).*Y)/sum(Y);

mu1(1)=mu11;
mu1(2)=mu12;%mu1 estimator


mu0=zeros(1,2);
mu01=sum(X(:,1).*(1.-Y))/sum(1.-Y);
mu02=sum(X(:,2).*1.-Y)/sum(1.-Y);

mu0(1)=mu01;
mu0(2)=mu02;%mu0 estimator

sigma1=zeros(2,2);

for i=1:n
    sigma1=Y(i)*(X(i,:)-mu1)'*(X(i,:)-mu1)+sigma1;
end

sigma1=sigma1/sum(Y);

sigma0=zeros(2,2);
for i=1:n
    sigma0=(1-Y(i))*(X(i,:)-mu0)'*(X(i,:)-mu0)+sigma0;
end

sigma0=sigma0/(n-sum(Y));


%%---------------------Plotting--------------------------------------------
figure(1),
plot(X(:,1),X(:,2),'+');

Mx=ceil(max(X(:,1)));
mx=floor(min(X(:,1)));
My=ceil(max(X(:,2)));
my=floor(min(X(:,1)));

x=linspace(mx,Mx);
y=linspace(my,My);
[X3,Y3] = meshgrid(x,y);



%%------------------Misclassification Error-----------------

r=zeros(n,1);
Z=densQ(X,Pi,mu1,mu0,sigma1,sigma0);
for i=1:n
    r(i)=min(Z(i),1-Z(i));
end

R1=mean(r);

%%-------------------Misclassification Error test data--------------------
Z=dlmread('classification_data_HWK1/classificationC.test');

Y=Z(:,3);
n=length(Y);

X=Z(:,1:2);

r=zeros(n,1);
Z=densQ(X,Pi,mu1,mu0,sigma1,sigma0);
for i=1:n
    r(i)=min(Z(i),1-Z(i));
end

R2=mean(r);
