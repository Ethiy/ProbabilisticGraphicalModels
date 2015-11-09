clear;

%%-----------------File Reading------------------------
Z=dlmread('classification_data_HWK1/classificationA.train');

Y=Z(:,3);
n=length(Y);

X=[Z(:,1:2),ones(n,1)];

%%-----------------Estimator--------------------------------
W=(X'*X)^-1*X'*Y;

sigma=sqrt(mean((Y-X*W).^2));

%%-----------------Plotting--------------------------------
Mx=ceil(max(X(:,1)));
mx=floor(min(X(:,1)));
My=ceil(max(X(:,2)));
my=floor(min(X(:,1)));

figure(1),
x=linspace(mx,Mx);
y=-W(1)/W(2)*x-W(3)/W(2)+1/2*1/W(2);
plot(x,y,X(:,1),X(:,2),'+');

%%------------------Misclassification Error training data-----------------

r=zeros(n,1);
Z=X*W;
for i=1:n
    r(i)=min(Z(i),1-Z(i));
end

R1=mean(r);

%%-------------------Misclassification Error test data--------------------
Z=dlmread('classification_data_HWK1/classificationA.test');

Y=Z(:,3);
n=length(Y);

X=[Z(:,1:2),ones(n,1)];

r=zeros(n,1);
Z=X*W;
for i=1:n
    r(i)=min(Z(i),1-Z(i));
end

R2=mean(r);