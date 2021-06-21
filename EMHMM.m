function [mu,sigma,trpro,delta] = EMHMM(data,m,mu,sigma,trpro,delta)
% EM estimation of a HMM
%This function implements the EM algorithm as described in Sections 2.1
% With the initial values lambda, gamma and delta for the natural parameters 
% In each iteration the logs of the forward and backward probabilities are computed
% The log- likelihood l, is needed in both E and M steps

maxiter=10000; %total number of iterations
tol=1e-6;      % acceptable threshold
if exist('delta')==0
   delta=double(stationary(trpro));
end
x=data';       % time series of normalized monthly river discharges
n=length(x);
munew=mu;
trpronew=trpro;
deltanew=delta;
sigmanew=sigma;
for iter=1:maxiter
    lprob=log(normpdf(repmat(x,1,m),repmat(mu,n,1),repmat(sigma,n,1)));    % lognormal of liklihood function
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

