%% k-nearest neighbor Estimation of the Conditional Entropy
% Computes the terms of the information decomposition in the univariate framework
% by using the observation matrix B

% INPUT: embedding matrix B, number of neighbors k, distance metric:
% 'maximum' Chebyshev distance (default)
% OUTPUT: estimates of conditional entropy Hy_Y, information storage Sy_Y,
% entropy Hy, and dynamic entropies HY and HyY

function out=unID_knn(B,k,metric)

    if ~exist('metric','var'), metric='chebychev'; end

    A=B(:,2:end);
    N=size(B,1);
    m=size(A,2)-1;

    %% kNN analysis
    %%% neighbor search in space of higher dimension    
    [~, distances] =  knnsearch(B,B,'K',k+1,'Distance',metric);
    dd = distances(:,end);

    %%% range search in the subspace of lower dimension - M_Y
    [~, distance_Y] =  knnsearch(A,A,'K',N,'Distance',metric);
    count_Y = distance_Y(:,2:end) < dd;
    count_Y = max(k-1, sum(count_Y,2));

    %%% range search in the subspace of lower dimension - M_y
    [~, distance_y] =  knnsearch(B(:,1),B(:,1),'K',N,'Distance',metric);
    count_y = distance_y(:,2:end) < dd;
    count_y = max(k-1, sum(count_y,2));

    %% computes CE and H

    dd2=2*dd;
    dd2(dd2==0)=[]; % do not accept distance=0
    Hy= psi(N)-(1/N)*( sum(psi(count_y+1)))+mean(log(dd2));
    HY= psi(N)-(1/N)*( sum(psi(count_Y+1)))+m*mean(log(dd2));
    HyY= psi(N)-psi(k)+(m+1)*mean(log(dd2));
    Hy_Y= - psi(k) + (1/N)*( sum(psi(count_Y+1))) + mean(log(dd2));
    % Hy_Y = HyY - HY;

    %%% OUTPUT
    out.Hy_Y=Hy_Y;
    out.Hy=Hy;
    out.HY=HY;
    out.HyY=HyY;
    out.Sy_Y=Hy-Hy_Y;


end