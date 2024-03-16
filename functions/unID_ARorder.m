%% Finds optimal VAR model order (Akaike Information Criterion and Bayesian Information Criterion)

%%% INPUT
% data: original time series 
% pmax: maximum model order
%%% OUTPUT
% pottaic: model order selected with AIC
% pottbic: model order selected with BIC

function out=VARorder(data,pmax)

[N,M]=size(data);

% figures of merit
aic=NaN*ones(pmax,1); bic=aic; detSu=aic;
for p=1:pmax   
    ret=unID_LinReg(data,(1:M),(1:M),(1:p)); 
    detSu(p)=det(ret.es2u);
    %formula multivariate AIC 
    aic(p)=log(detSu(p))+2*M*M*p/N; % S covariance matrix
    %formula multivariate BIC
    bic(p)=log(detSu(p))+log(N)*M*M*p/N; % S covariance matrix
end
pottaic=find(aic == min(aic));
pottbic=find(bic == min(bic));

out.pottaic=pottaic;
out.pottbic=pottbic;
out.aic=aic;
out.bic=bic;
out.detSu=detSu;

end
