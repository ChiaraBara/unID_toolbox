%% Complexity analysis for a simulated univariate system
clear all; close all; clc;

addpath([pwd,'/functions/']);

%% Parameters of the simulation
N = 300; % time series length
vecr = [0:0.1:0.9]; % strength of the stochastic oscillation
f = 0.25; % frequency of the stochastic oscillation
M = 1; %n. of time series
p = 3; % maximum lag

%% Parameter of the estimators
base = 2; % 0: natural, 2: bits
m = 3; %memory of the process 
k = 10; %nearest neighbor: number of neighbors
r = 0.3; %kernel: threshold distance
b = 4; %binning: number of bins
delta = 1e-3; %slope: 1st threshold
gamma = 1; %slope: 2nd threshold

%% computation
for ir = 1:length(vecr)
    
    %%% Simulation setup
    par.poles=([vecr(ir) f]); % Oscillation
    par.coup=[]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
    par.Su=1; %variance of innovation processes
    
    %%% Theorical entropy rate from AR parameters
    [Am,Su] = var_simulations(M,par); % parameters
    ret = its_CElinVAR1(Am,Su,p); % exact values information dynamics
    varY = ret.Sy;
    CETh(ir) = 0.5*log(2*pi*exp(1)/varY);

    %%% Simulation
    Un = mvnrnd(zeros(1,M),Su,N); %white gaussian innovations
    Yn = var_filter(Am,Un); % realization (filters the noise)
    Yn = zscore(Yn); %normalization to zero mean and unit variance

    %%% Estimated measure
    V = [ones(m,1),(1:m)'];
    B = unID_buildvectors(Yn,1,V); %observation matrix (present and past)
    
    % Linear Entropy rate estimation
    outlin = unID_lin(B); 
    CELin(ir) = outlin.Hy_Y; 
    HLin(ir) = outlin.Hy; 
    HLin2(ir) = unID_Hlin(Yn);
    
    % knn Entropy rate estimation
    outknn = unID_knn(B,k);
    CEKnn(ir) = outknn.Hy_Y; 
    HKnn(ir) = outknn.Hy; 
    HKnn2(ir) = unID_Hknn(Yn,k);
    
    % Kernel Entropy rate estimation
    outker = unID_ker(B,r,'c'); %NOTE: in the case that Yn is not normalized to unit variance, impose r=r*std(Yn)
    CEKer(ir) = outker.Hy_Y;
    HKer(ir) = outker.Hy;
    HKer2(ir) = unID_Hker(Yn,r,'c');
    
    % Binning Entropy rate estimation
    outbin = unID_bin(B,b,base);
    CEBin(ir)=outbin.Hy_Y;
    HBin(ir)=outbin.Hy;
    HBin2(ir) = unID_Hbin(Yn,b,base);
    
    % Permutation Entropy rate estimation
    outperm = unID_perm(B,base);
    CEPerm(ir) = outperm.Hy_Y;
    
    % Slope Entropy rate estimation
    outslope = unID_slope(B,delta,gamma,base);
    CESlope(ir) = outslope.Hy_Y;
    
end

%% plot

figure(1);
% lin
a1 = subplot(2,3,1);
plot(vecr,CELin,'.-');
hold on; plot(vecr,CETh,'-k');
ylabel('[nats]');
xlabel('\rho');
title('CE_{lin}');
% knn
a2 = subplot(2,3,2);
plot(vecr,CEKnn,'.-');
hold on; plot(vecr,CETh,'-k');
ylabel('[nats]');
xlabel('\rho');
title('CE_{knn}');
% ker
a3 = subplot(2,3,3);
plot(vecr,CEKer,'.-');
ylabel('[nats]');
xlabel('\rho');
title('CE_{ker}');
% bin
a4 = subplot(2,3,4);
plot(vecr,CEBin,'.-');
ylabel('[bits]');
xlabel('\rho');
title('CE_{bin}');
% perm
a5 = subplot(2,3,5);
plot(vecr,CEPerm,'.-');
ylabel('[bits]');
xlabel('\rho');
title('CE_{perm}');
% slope
a6 = subplot(2,3,6);
plot(vecr,CESlope,'.-');
ylabel('[bits]');
xlabel('\rho');
title('CE_{slope}');
linkaxes([a1 a2 a3 a4 a5 a6],'y');

figure(2);
% lin
a1 = subplot(2,2,1);
plot(vecr,HLin,'.-');
hold on; plot(vecr,HLin2,'.-');
ylabel('[nats]');
xlabel('\rho');
title('H_{lin}');
% knn
a2 = subplot(2,2,2);
plot(vecr,HKnn,'.-');
hold on; plot(vecr,HKnn2,'.-');
ylabel('[nats]');
xlabel('\rho');
title('H_{knn}');
% ker
a3 = subplot(2,2,3);
plot(vecr,HKer,'.-');
hold on; plot(vecr,HKer2,'.-');
ylabel('[nats]');
xlabel('\rho');
title('H_{ker}');
% bin
a4 = subplot(2,2,4);
plot(vecr,HBin,'.-');
hold on; plot(vecr,HBin2,'.-');
ylabel('[bits]');
xlabel('\rho');
title('H_{bin}');
linkaxes([a1 a2 a3 a4],'y');
