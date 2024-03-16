%% Kernel Estimation of the Conditional Entropy
% Computes the terms of the information decomposition in the univariate framework
% by using the observation matrix B

% INPUT: embedding matrix B, distance threshold r, distance metric:
% 'maximum' Chebyshev distance ['c'](default), normalized distance ['r'],
% Euclidean distance ['e']
% OUTPUT: estimates of conditional entropy Hy_Y, information storage Sy_Y,
% entropy Hy, and dynamic entropies HY and HyY

function out=unID_ker(B,r,norma)

    if ~exist('norma','var') 
        norma='c'; %default Chebyshev
    end

    dataMat=flipdim(B',1);
    dim=size(B,2)-1;
    N=size(B,1)+dim;

    correl = zeros(1,3);
    vecm = [1,dim,dim+1];
    count = zeros(length(vecm),N-dim);
    for m = 1:length(vecm)
        if m == 1
            tempMat = dataMat(end,:);
        else
            tempMat = dataMat(1:vecm(m),:);
        end

        for i = 1:N-dim

            tmp = tempMat; tmp(:,i)=[]; % excluding self-matching
            if norma=='c' % Chebyshev distance
                dist = max(abs(tmp - repmat(tempMat(:,i),1,N-dim-1)),[],1);
            elseif norma=='r' % normalized distance (Range Entropy)
                tmp1 = max(abs(tmp - repmat(tempMat(:,i),1,N-dim-1)),[],1); 
                tmp2 = min(abs(tmp - repmat(tempMat(:,i),1,N-dim-1)),[],1); 
                dist = (tmp1-tmp2)./(tmp1+tmp2);
            elseif norma=='e' % Euclidean distance
                dist=nan*ones(1,size(tmp,2));
                for q=1:size(tmp,2)
                    dist(q)=norm(tempMat(:,i)-tmp(:,q),2);
                end
            end

            % calculate Heaviside function of the distance
            D = (dist < r);

            count(m,i) = sum(D)/(N-dim-1); %computes probability
        end

        correl(m) = sum(count(m,:))/(N-dim); %average of probability
    end

    Hy = -log(correl(1));
    HY = -log(correl(2));
    HyY = -log(correl(3));
    Hy_Y = -log(correl(3)/correl(2));

    %%% OUTPUT
    out.Hy_Y=Hy_Y;
    out.Hy=Hy;
    out.HY=HY;
    out.HyY=HyY;
    out.Sy_Y=Hy-Hy_Y;

end