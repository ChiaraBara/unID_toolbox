%% k-nearest neighbor Estimation of entropy
% Computes entropy by using the observation matrix B (which coincides with the signal)

% INPUT: embedding matrix B, number of neighbors k, distance metric:
% 'maximum' Chebyshev distance (default)
% OUTPUT: entropy Hy

function Hy=unID_Hknn(B,k,metric)

    if ~exist('metric','var'), metric='maximum'; end
    N=size(B,1);

    %% kNN analysis
    %%% neighbor search in space of higher dimension
    atria_yY = nn_prepare(B, metric);
    [~, distances] = nn_search(B, atria_yY, (1:N)', k, 0);
    dd=distances(:,k);

    %% computes H
    dd2=2*dd; dd2(dd2==0)=[]; % do not accept distance=0
    Hy= psi(N)-psi(k)+mean(log(dd2));

end
