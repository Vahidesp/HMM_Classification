function [ la,lb ] = alfabeta(data,m,mu,sigma,trpro,delta)
% This function computes the logarithms of the forward and backward probabilities
% Please refer to Zucchini et al., (2017) Equations (4.1) and (4.2) on p. 60 
x=data';
n=length(x);
alpha=zeros(m,n);beta=zeros(m,n);
allprob=normpdf(repmat(x,1,m),repmat(mu,n,1),repmat(sigma,n,1));
foo=delta'.*allprob(1,:);
sumfoo=sum(foo);
lscale=log(sumfoo);
foo=foo/sumfoo;
alpha(:,1)= log(foo)+lscale';

for i=2:n
    foo=(foo*trpro).*allprob(i,:);
    sumfoo=sum(foo);
    lscale=lscale+log(sumfoo);
    foo=foo/sumfoo;
    alpha(:,i)=log(foo)+lscale;
end
beta(:,n)=zeros(1,m);
foo=repmat(1/m,1,m);
lscale=log(m);
for i=n-1:-1:1
    foo=trpro*(allprob(i+1,:).*foo)';
    beta(:,i)=log(foo)+lscale;
    sumfoo=sum(foo);
    foo=(foo/sumfoo)';
    lscale=lscale+log(sumfoo);
end

la=alpha;lb=beta;
end

