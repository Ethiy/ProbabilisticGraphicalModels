clear;
%%Read the data------------------------------------------------------------
X=dlmread('C:\Users\Oussama\Desktop\DM2\classification_data_HWK2\EMGaussian.data');
n=length(X);

%%a-Kmeans-----------------------------------------------------------------

Mu0=[X(50,:);X(143,:);X(299,:);X(450,:)];

TK=20;

for t=1:TK
    Mu0=centers(X,Mu0);
end

X1=cluster(1,X,Mu0);
X2=cluster(2,X,Mu0);
X3=cluster(3,X,Mu0);
X4=cluster(4,X,Mu0);

figure(1),
plot(X1(:,1),X1(:,2),'.r',X2(:,1),X2(:,2),'.g',X3(:,1),X3(:,2),'.y',X4(:,1),X4(:,2),'.b',Mu0(:,1),Mu0(:,2),'+k');

%%b-EM-gaussian-mixture-prop-identity--------------------------------------

tau=indicatrice(X,Mu0);
[Pi , Mu , Sigma]=theta1(X,tau);

TEM1=150;

for t=1:TEM1
    tau=Q(X,Pi , Mu , Sigma);
    [Pi , Mu , Sigma]=theta1(X,tau);
end

X1=cluster(1,X,Mu);
X2=cluster(2,X,Mu);
X3=cluster(3,X,Mu);
X4=cluster(4,X,Mu);

u = -10:0.125:10;
v = -10:0.125:10;
[U,V] = meshgrid(u,v);
b = [3,3];

s1=Sigma(:,:,1)^-1;
W1=s1(1,1).*(U-Mu(1,1)).^2+s1(2,2).*(V-Mu(1,2)).^2+2*s1(2,1).*(U-Mu(1,1)).*(V-Mu(1,2));

s2=Sigma(:,:,2)^-1;
W2=s2(1,1).*(U-Mu(2,1)).^2+s2(2,2).*(V-Mu(2,2)).^2+2*s2(2,1).*(U-Mu(2,1)).*(V-Mu(2,2));

s3=Sigma(:,:,1)^-1;
W3=s3(1,1).*(U-Mu(3,1)).^2+s3(2,2).*(V-Mu(3,2)).^2+2*s3(2,1).*(U-Mu(3,1)).*(V-Mu(3,2));

s4=Sigma(:,:,1)^-1;
W4=s4(1,1).*(U-Mu(4,1)).^2+s4(2,2).*(V-Mu(4,2)).^2+2*s4(2,1).*(U-Mu(4,1)).*(V-Mu(4,2));

figure(2),
plot(X1(:,1),X1(:,2),'.r',X2(:,1),X2(:,2),'.g',X3(:,1),X3(:,2),'.y',X4(:,1),X4(:,2),'.b',Mu(:,1),Mu(:,2),'+k');
hold on
contour(U,V,W1,b,'r');
contour(U,V,W2,b,'g');
contour(U,V,W3,b,'y');
contour(U,V,W4,b,'b');

%%c-EM-gaussian-mixture-general--------------------------------------------

tau=indicatrice(X,Mu0);
[Pi , Mu , Sigma]=theta2(X,tau);

TEM1=150;

for t=1:TEM1
    tau=Q(X,Pi , Mu , Sigma);
    [Pi , Mu , Sigma]=theta2(X,tau);
end

X1=cluster(1,X,Mu);
X2=cluster(2,X,Mu);
X3=cluster(3,X,Mu);
X4=cluster(4,X,Mu);

u = -10:0.125:10;
v = -10:0.125:10;
[U,V] = meshgrid(u,v);
b = [3,3];

s1=Sigma(:,:,1)^-1;
W1=s1(1,1).*(U-Mu(1,1)).^2+s1(2,2).*(V-Mu(1,2)).^2+2*s1(2,1).*(U-Mu(1,1)).*(V-Mu(1,2));

s2=Sigma(:,:,2)^-1;
W2=s2(1,1).*(U-Mu(2,1)).^2+s2(2,2).*(V-Mu(2,2)).^2+2*s2(2,1).*(U-Mu(2,1)).*(V-Mu(2,2));

s3=Sigma(:,:,1)^-1;
W3=s3(1,1).*(U-Mu(3,1)).^2+s3(2,2).*(V-Mu(3,2)).^2+2*s3(2,1).*(U-Mu(3,1)).*(V-Mu(3,2));

s4=Sigma(:,:,1)^-1;
W4=s4(1,1).*(U-Mu(4,1)).^2+s4(2,2).*(V-Mu(4,2)).^2+2*s4(2,1).*(U-Mu(4,1)).*(V-Mu(4,2));


figure(3),
plot(X1(:,1),X1(:,2),'.r',X2(:,1),X2(:,2),'.g',X3(:,1),X3(:,2),'.y',X4(:,1),X4(:,2),'.b',Mu(:,1),Mu(:,2),'+k');

hold on
contour(U,V,W1,b,'r');
contour(U,V,W2,b,'g');
contour(U,V,W3,b,'y');
contour(U,V,W4,b,'b');
