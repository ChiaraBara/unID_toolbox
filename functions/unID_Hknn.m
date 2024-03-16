%% k-nearest neighbor Estimation of the Conditional Entropy
% Computes entropy by using the observation matrix B (which coincide with the signal)

% INPUT: embedding matrix B, number of neighbors k, distance metric:
% 'maximum' Chebyshev distance (default)
% OUTPUT: entropy Hy

function Hy=unID_Hknn(B,k,metric)

    if ~exist('metric','var'), metric='chebychev'; end
    N=size(B,1);

    %% kNN analysis
    %%% neighbor search in space of higher dimension    
    [~, distances] =  knnsearch(B,B,'K',k+1,'Distance',metric);
    dd = distances(:,end);

    %% computes H
    dd2=2*dd; dd2(dd2==0)=[]; % do not accept distance=0
    Hy= psi(N)-psi(k)+mean(log(dd2));

end