function [mu,sigma,delta,trpro,locstate] = HMM(data,m,mu,sigma,trpro,delta)
%********************************************************************************************
% This function performs Hidden Markov Model classification of                              *
% for the time series of monthly river discharges                                           *
% The HMM parameters are calculated with EMHMM function usig EM algorithm and               *
% determine the sequence of hidden climate states with Viterbi algorithm(Viterbi, 1967)     *
% i.e. determines for each time t the most likely state                                     *
% m       :is the number of climate states                                                  *
% data    :is the time series of monthly streamflows                                        *
% mu      :is the initial mean for each climate state (optinal)                             *
% sigma   :is the initial standard deviation for each climate state(optinal)                *
% trpro   :is the initial transition probability of HMM (optinal)                           *
% delta   :is the initial distribution of hidden climate states (optinal)                   *
%locstate :is the classification of monthly river discharges                                *
%********************************************************************************************

% initialization of the HMM parameters
%% 
innfac=quantile(data,m);  

if m==2
    
    [~,yi1]=find(data<innfac(2));
    [~,yi2]=find(data>=innfac(2));
    muu1=mean(data(yi1));muu2=mean(data(yi2));
    mu=[muu1 muu2];
    sigma1=std(data(yi1));sigma2=std(data(yi2));
    sigma=[sigma1 sigma2];
elseif m==3
       [~,yi1]=find(data<innfac(1));
       [~,yi2]=find(data>=innfac(1)& data<innfac(3));
       [~,yi3]=find(data>=innfac(3));
       muu1=mean(data(yi1));muu2=mean(data(yi2));muu3=mean(data(yi3));mu=[muu1 muu2 muu3];
       sigma1=std(data(yi1));sigma2=std(data(yi2));sigma3=std(data(yi3));sigma=[sigma1 sigma2 sigma3];
else
    for i=1:m
        mu(1,i)=min(data) + mean(data)*rand();
        sigma(1,i)=min(data)+ std(data)*rand();
    end
    mu=sort(mu);sigma=(sort(sigma));
end 

trpro=ones(m)*(1/m);
delta=double(stationary(trpro));
%% 

%% 
% EM estimation of a HMM
%This function implements the EM algorithm as described in Sections 2.1
% With the initial values lambda, gamma and delta for the natural parameters 
% In each iteration the logs of the forward and backward probabilities are computed
% The log- likelihood l, is needed in both E and M steps
[mu,sigma,trpro,delta] = EMHMM(data,m,mu,sigma,trpro,delta);
%% 

%%
% This section performs Viterbi algorithm (Viterbi, 1967)
% computes the most likely sequence of states, given the parameters and the observations
% Forward and Backward probabilities 
% This function computes the logarithms of the forward and backward probabilities
% Please refer to Zucchini et al., (2017) Equations (4.1) and (4.2) on p. 60 
[la,lb]=alfabeta(data,m,mu,sigma,trpro,delta);



x=data';
n=length(x);
co_state=zeros(n,m);
locstate=zeros(1,n);
y=zeros(1,n);
c=max(la(:,n));
llk=c+log(sum(exp(la(:,n)-c)));

for i=1:n
    co_state(i,:)=exp(la(:,i)+lb(:,i)-llk)';
    [~,locstate(i)]=max(co_state(i,:));
end

for i=1:m
    y(locstate==i)=mu(i);
end


amonth=[data' y' locstate'];
end

