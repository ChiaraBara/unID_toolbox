%% Linear Gaussian Estimation of entropy
% Computes entropy by using the observation matrix B (which coincides with the signal)

% INPUT: embedding matrix B
% OUTPUT: entropy Hy

function Hy=unID_Hlin(B)


    %% compute Conditional entropy
    Yb=B'; % inversion works with data organized in rows

    %Entropy measures estimates for the Gaussian case (Barnett PRL 2009):
    Hy=0.5*size(Yb,1)*log(2*pi*exp(1)*cov(Yb));

end
