%% Binning Estimation of entropy
% Computes entropy by using the observation matrix B (which coincides with the signal)

% INPUT: embedding matrix B, number of quantization bins b, log base 
% OUTPUT: entropy Hy

function Hy=unID_Hbin(B,b,base)

    Bq=unID_quantization(B,b);
    Hy=unID_H(Bq,base);

end
