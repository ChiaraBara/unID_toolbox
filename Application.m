clear all; close all; clc;

addpath([pwd,'/functions/']);

%%% load data
load('data_RR.mat');
pfilter=0.94; %filter parameter

%%% Parameter of estimators
base = 2; % 0: nats, 2: bits
pmax = 10;
m_knn = 3; m_ker = 3; m_bin = 3; m_perm = 4; m_slope = 3; %memory of the process 
k = 10; %nearest neighbor: number of neighbors
r = 0.3; %kernel: threshold distance
b = 6; %binning: number of bins
delta = 1e-3; %slope: 1st threshold
gamma = 1; %slope: 2nd threshold


%% Analysis
% pre-processing: AR highpass filter
Sf = detrend_AR_filter(data,1,pfilter); % AR highpass filtered series       
S = zscore(Sf); % normalization to zero mean and unit variance

% lin
out_p = unID_ARorder(S,pmax); %selection of optimal model order
V_lin = [ones(out_p.pottbic,1),(1:out_p.pottbic)']; 
B_lin = unID_buildvectors(S,1,V_lin);
outlin = unID_lin(B_lin);
CElin = outlin.Hy_Y;
Hlin1 = outlin.Hy;
Hlin2 = unID_Hlin(S);

% knn
V_knn = [ones(m_knn,1),(1:m_knn)']; 
B_knn = unID_buildvectors(S,1,V_knn);
outknn = unID_knn(B_knn,k);
CEknn = outknn.Hy_Y;
Hknn1 = outknn.Hy;
Hknn2 = unID_Hknn(S,k);

% ker
V_ker = [ones(m_ker,1),(1:m_ker)']; 
B_ker = unID_buildvectors(S,1,V_ker);
outker = unID_ker(B_ker,r,'c');
CEker = outker.Hy_Y; 
Hker1 = outker.Hy;
Hker2 = unID_Hker(S,r,'c');

% bin
V_bin = [ones(m_bin,1),(1:m_bin)']; 
B_bin = unID_buildvectors(S,1,V_bin);
outbin = unID_bin(B_bin,b,base);
CEbin = outbin.Hy_Y;
Hbin1 = outbin.Hy;
Hbin2 = unID_Hbin(S,b,base);

% perm
V_perm = [ones(m_perm,1),(1:m_perm)']; 
B_perm = unID_buildvectors(S,1,V_perm);
outperm = unID_perm(B_perm,base);
CEperm = outperm.Hy_Y;

% slope
V_slope = [ones(m_slope,1),(1:m_slope)']; 
B_slope = unID_buildvectors(S,1,V_slope);
outslope = unID_slope(B_slope,delta,gamma,base); 
CEslope = outslope.Hy_Y;

%% display
disp('Estimated values of CE:');
disp(['Lin: ', num2str(CElin),' nats']);
disp(['Knn: ', num2str(CEknn),' nats']);
disp(['Ker: ', num2str(CEker),' nats']);
disp(['Bin: ', num2str(CEbin),' bits']);
disp(['Perm: ', num2str(CEperm),' bits']);
disp(['Slope: ', num2str(CEslope),' bits']);

disp('Estimated values of H:');
disp(['Lin: ', num2str(Hlin1),' nats,', num2str(Hlin2),' nats']);
disp(['Knn: ', num2str(Hknn1),' nats,', num2str(Hknn2),' nats']);
disp(['Ker: ', num2str(Hker1),' nats,', num2str(Hker2),' nats']);
disp(['Bin: ', num2str(Hbin1),' bits,', num2str(Hbin2),' bits']);
