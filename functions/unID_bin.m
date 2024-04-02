%% Binning Estimation of the Conditional Entropy
% Computes the terms of the information decomposition in the univariate framework
% by using the observation matrix B

% INPUT: embedding matrix B, number of quantization bins b, log base 
% OUTPUT: estimates of conditional entropy Hy_Y, information storage Sy_Y,
% entropy Hy, and dynamic entropies HY and HyY

function out=unID_bin(B,b,base)

    M=size(B,2);
    for im=1:M
        Bq(:,im)=unID_quantization(B(:,im),b);
    end
    Aq=Bq(:,2:end);
    
    HyY=unID_H(Bq,base);
    HY=unID_H(Aq,base);
    Hy=unID_H(Bq(:,1),base);
    Hy_Y=HyY-HY;
    
    %%% OUTPUT
    out.Hy_Y=Hy_Y;
    out.Hy=Hy;
    out.HY=HY;
    out.HyY=HyY;
    out.Sy_Y=Hy-Hy_Y;

end
