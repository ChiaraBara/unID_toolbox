%% Permutation-based Estimation of the Conditional Entropy
% Computes the terms of the information decomposition in the univariate framework
% by using the observation matrix B

% INPUT: embedding matrix B, log base 
% OUTPUT: estimates of conditional entropy Hy_Y and dynamic entropies HY and HyY

function out=unID_perm(B,base)

    A=B(:,2:end);
    
    %%% ranking procedure
    [~,Bp]=sort(B,2); [~,Bp]=sort(Bp,2); 
    [~,Ap]=sort(A,2); [~,Ap]=sort(Ap,2);
    
    HyY=unID_H(Bp,base);
    HY=unID_H(Ap,base);
    Hy_Y=HyY-HY;
    
    %%% OUTPUT
    out.Hy_Y=Hy_Y;
    out.HY = HY;
    out.HyY = HyY;

end
