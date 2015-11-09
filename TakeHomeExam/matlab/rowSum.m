function Y=rowSum(X)

Y=zeros(size(X(:,1)));

for i=1:length(X(1,:))
  Y=Y+X(:,i);
end
