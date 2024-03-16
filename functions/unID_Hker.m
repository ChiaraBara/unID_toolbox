%% Kernel Estimation of the Conditional Entropy
% Computes entropy by using the observation matrix B (which coincide with the signal)

% INPUT: embedding matrix B, distance threshold r, distance metric:
% 'maximum' Chebyshev distance ['c'](default), normalized distance ['r'],
% Euclidean distance ['e']
% OUTPUT: entropy Hy

function Hy=unID_Hker(B,r,norma)

    if ~exist('norma','var') 
        norma='c'; %default Chebyshev
    end

    dataMat=flipdim(B',1);
    N=size(B,1);
    count = zeros(1,N);
    
    for i = 1:N

        tmp = dataMat; tmp(:,i)=[]; % excluding self-matching
        if norma=='c' % Chebyshev distance
            dist = max(abs(tmp - repmat(dataMat(:,i),1,N-1)),[],1);
        elseif norma=='r' % normalized distance (Range Entropy)
            tmp1 = max(abs(tmp - repmat(dataMat(:,i),1,N-1)),[],1); 
            tmp2 = min(abs(tmp - repmat(dataMat(:,i),1,N-1)),[],1); 
            dist = (tmp1-tmp2)./(tmp1+tmp2);
        elseif norma=='e' % Euclidean distance
            dist=nan*ones(1,size(tmp,2));
            for q=1:size(tmp,2)
                dist(q)=norm(dataMat(:,i)-tmp(:,q),2);
            end
        end

        % calculate Heaviside function of the distance
        D = (dist < r);

        count(i) = sum(D)/(N-1); %computes probability
    end

    correl = sum(count)/(N); %average of probability
   
    Hy = -log(correl);

end