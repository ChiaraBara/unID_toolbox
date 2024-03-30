%% Entropy of a discrete multidimensional variable A

%%% INPUT
% A: Np*d matrix, realizations discrete d-dimensional vector variable, Np observations
% base: base of the logarithm (2, or not pass argument to measure in bits, 0 to measure in nats)
%%% OUTPUT
% H: entropy esimates by using the frequentistic approach

function H=unID_H(A,base)

    if nargin<2, base=2; end %default entropy in bits

    Np=size(A,1); %total number of patterns
    [C,~,ic] = unique(A,'rows');
    Np0=size(C,1); %number of patterns with nonzero probability
    for cnt=1:Np0
        p(cnt)=sum(ic==cnt)/Np;
    end

    switch base
        case 2
            H=sum(-p.*log2(p));
        case 0
            H=sum(-p.*log(p));
    end
    
end
