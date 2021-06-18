function [mu,sigma,trpro,delta] = EMHMM(data,m,mu,sigma,trpro,delta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%trpro=[0.7 0.3;0.1 0.9];
maxiter=10000;tol=1e-6;
if exist('delta')==0
   delta=double(stationary(trpro));
end
x=data';
n=length(x);
munew=mu;
trpronew=trpro;
deltanew=delta;
sigmanew=sigma;
for iter=1:maxiter
    lprob=log(normpdf(repmat(x,1,m),repmat(mu,n,1),repmat(sigma,n,1)));
    [la,lb]=alfabeta(data,m,mu,sigma,trpro,delta);
    c=max(la(:,n));
    llk=c+log(sum(exp(la(:,n)-c)));
    for j=1:m
        for k=1:m
            trpronew(j,k)=trpro(j,k).*sum(exp(la(j,1:(n-1))+(lprob(2:n,k))'+lb(k,2:n)-llk));
        end
        munew(j)=sum(exp(la(j,:)+lb(j,:)-llk).*x')/sum(exp(la(j,:)+lb(j,:)-llk));
        sigmanew(j)=sqrt(sum(exp(la(j,:)+lb(j,:)-llk).*((x-mu(j))').^2)/sum(exp(la(j,:)+lb(j,:)-llk)));
    end
    trpronew=bsxfun(@rdivide,trpronew,sum(trpronew,2));
    deltanew=exp(la(:,1)+lb(:,1)-llk);
    deltanew=deltanew/sum(deltanew);
    crit=sum(sum(abs(mu-munew)))+sum(sum(abs(trpro-trpronew)))+sum(sum(abs(delta-deltanew),2))...
         +sum(sum(abs(sigma-sigmanew),2));
     if crit<tol
         np=m.*m+m-1;
         AIC=-2.*(llk-np);
         BIC=-2.*llk+np*log(n);
         mllk=-llk;
     end
     mu=munew;
     sigma=sigmanew;
     trpro=trpronew;
     delta=deltanew;    
end
end

