%% k-nearest neighbor Estimation of the Conditional Entropy
% Computes the terms of the information decomposition in the univariate framework
% by using the observation matrix B

% INPUT: embedding matrix B, number of neighbors k, distance metric:
% 'maximum' Chebyshev distance (default)
% OUTPUT: estimates of conditional entropy Hy_Y, information storage Sy_Y,
% entropy Hy, and dynamic entropies HY and HyY

function out=unID_knn(B,k,metric)

    if ~exist('metric','var'), metric='maximum'; end

    A=B(:,2:end);
    N=size(B,1);
    m=size(A,2);

    %% kNN analysis
    %%% neighbor search in space of higher dimension
    atria_yY = nn_prepare(B, metric);
    [~, distances] = nn_search(B, atria_yY, (1:N)', k, 0);
    dd=distances(:,k);

    %%% range search in the subspace of lower dimension - M_Y
    if ~isempty(A)
        atria_Y = nn_prepare(A, metric);
        [count_Y, tmp] = range_search(A, atria_Y, (1:N)', dd, 0);
        tmp=tmp(:,2);
        for n=1:length(tmp)
            count_Y(n)=max(k-1,count_Y(n)-sum(tmp{n}==dd(n)));
        end
    else
        count_Y=(N-1)*ones(N,1);
    end
    %%% range search in the subspace of lower dimension - M_y
    if ~isempty(B(:,1))
        atria_y = nn_prepare(B(:,1), metric);
        [count_y, tmp] = range_search(B(:,1), atria_y, (1:N)', dd, 0);
        tmp=tmp(:,2);
        for n=1:length(tmp)
            count_y(n)=max(k-1,count_y(n)-sum(tmp{n}==dd(n)));
        end
    else
        count_y=(N-1)*ones(N,1);
    end


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
