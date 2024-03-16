%% Linear Gaussian Estimation of the Conditional Entropy
% Computes the terms of the information decomposition in the univariate framework
% by using the observation matrix B

% INPUT: embedding matrix B
% OUTPUT: estimates of conditional entropy Hy_Y, information storage Sy_Y,
% entropy Hy, and dynamic entropies HY and HyY

function out=unID_lin(B)


    %% compute Conditional entropy
    Yb=B(:,1)'; % inversion works with data organized in rows
    A=B; A(:,1)=[]; Z=A';
    
    m=size(A,2);

    Am=Yb/Z; % least squares!

    Yp=Am*Z; 
    Up=Yb-Yp;
    S=cov(Up');

    %Entropy measures estimates for the Gaussian case (Barnett PRL 2009):
    Hy=0.5*size(Yb,1)*log(2*pi*exp(1)*cov(Yb));
    HY=0.5*log(det(cov(A)))+0.5*size(Yb,1)*log((2*pi*exp(1))^m); 
    HyY=0.5*log(det(cov(B)))+0.5*size(Yb,1)*log((2*pi*exp(1))^(m+1)); 
    Hy_Y=0.5*log(det(S))+0.5*size(Yb,1)*log(2*pi*exp(1)); 
    % Hy_Y = HyY - HY;
    
    %%% OUTPUT
    out.Hy_Y=Hy_Y;
    out.Hy=Hy;
    out.HY=HY;
    out.HyY=HyY;
    out.Sy_Y=Hy-Hy_Y;

end