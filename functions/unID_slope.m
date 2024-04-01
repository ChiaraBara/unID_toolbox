%% Slope-based Estimation of the Conditional Entropy
% Computes the terms of the information decomposition in the univariate framework
% by using the observation matrix B

% INPUT: embedding matrix B, discretization parameters delta and gamma, log base 
% OUTPUT: estimates of conditional entropy Hy_Y and dynamic entropies HY and HyY

function out = unID_slope(B,delta,gamma,base)

    Bdiff = -diff(B,1,2);  
    [N,M]=size(Bdiff);

    %%% slope vector
    Bs = zeros(size(Bdiff));
    %%% discretization procedure
    for in = 1:N
        for im = 1:M           
            if abs(Bdiff(in,im))<=delta
                Bs(in,im) = 0;
            elseif (Bdiff(in,im)>delta) && (Bdiff(in,im)<=gamma)
                Bs(in,im) = 1;
            elseif Bdiff(in,im)>gamma
                Bs(in,im) = 2;
            elseif (Bdiff(in,im)>=-gamma) && (Bdiff(in,im)<-delta)
                Bs(in,im) = -1;
            elseif Bdiff(in,im)<-gamma
                Bs(in,im) = -2;
            end   
        end
    end

    As = Bs; As(:,1)=[];

    HyY = unID_H(Bs,base);
    HY = unID_H(As,base);
    Hy_Y = HyY - HY;
    
    %%% OUTPUT
    out.Hy_Y=Hy_Y;
    out.HY = HY;
    out.HyY = HyY;

end
