clc
clear all

%% Variable Initiation

K=4;

solution_05nov2014;
clc
clf;
Sigmas=sigmas;
means=mu;
U1 = x;
U2=xtest;
Pi=ones(1,K)./K;
A=.5.*eye(K)+(1/6.*(ones(K)-eye(K)));

%% Forward Backward

[Alpha Beta scale NormDist]=forwardBackward(U1,K,A,Pi,means,Sigmas);

Gamma=(Alpha.*Beta); 
Gamma=rowDivision(Gamma,rowSum(Gamma));

T=length(U1);
Xi=zeros(T-1,K*K);
for i=1:T-1
    t=A.*( Alpha(i,:)' * (Beta(i+1,:).*NormDist(i+1,:)));
    Xi(i,:)=t(:)'/sum(t(:));
end;

%% Question 2

time=1:100;

figure(1),
subplot(2,2,1)
plot(time,Gamma(time,1))
title('state=1')

subplot(2,2,2)
plot(time,Gamma(time,2),'r')
title('state=2')

subplot(2,2,3)
plot(time,Gamma(time,3),'g')
title('state=3')

subplot(2,2,4)
plot(time,Gamma(time,4),'y')
title('state=4')

%% Question 4

[means1,Sigmas1,A1,Pi1,LogLik1]=hmm(K,U1,50,mu,sigmas);

[means2,Sigmas2,A2,Pi2,LogLik2]=hmm(K,U2,50,mu,sigmas);

%% Question 5

figure(2),
plot(1:50,LogLik1,1:50,LogLik2);

%% Question 7

q=viterbi(K,Pi1,U1,A1,means1,Sigmas1);

I1=[];
I2=[];
I3=[];
I4=[];

for i=1:length(q)
    if q(i)==1
        I1=[I1,i];
    elseif q(i)==2
        I2=[I2,i];
    elseif q(i)==3
        I3=[I3,i];
    else
        I4=[I4,i];
    end
end

figure(3),
scatter(U1(I1,1),U1(I1,2));
hold on
scatter(U1(I2,1),U1(I2,2),'r');
scatter(U1(I3,1),U1(I3,2),'g');
scatter(U1(I4,1),U1(I4,2),'k');

hold off;

%% Question 8

[Alpha Beta scale NormDist]=forwardBackward(U2,K,A1,Pi1,means1,Sigmas1);

Gamma=(Alpha.*Beta); 
Gamma=rowDivision(Gamma,rowSum(Gamma));

T=length(U2);
Xi=zeros(T-1,K*K);
for i=1:T-1
    t=A.*( Alpha(i,:)' * (Beta(i+1,:).*NormDist(i+1,:)));
    Xi(i,:)=t(:)'/sum(t(:));
end

time=1:100;

figure(4),
subplot(2,2,1)
plot(time,Gamma(time,1))
title('state=1')

subplot(2,2,2)
plot(time,Gamma(time,2),'r')
title('state=2')

subplot(2,2,3)
plot(time,Gamma(time,3),'g')
title('state=3')

subplot(2,2,4)
plot(time,Gamma(time,4),'y')
title('state=4')

%% Question 9

P=Gamma(time,:);

[l states]=max(P,[],2);

figure(5),
plot(states,'*r');
hold on
%% Question 10

q2=viterbi(K,Pi1,U2,A1,means1,Sigmas1);

plot(q2(1:100),'-');
hold off