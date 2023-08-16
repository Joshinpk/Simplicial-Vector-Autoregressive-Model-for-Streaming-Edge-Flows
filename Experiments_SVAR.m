close all;clear all;
%%%PARAMETERS
mu1=0.01;%L1_lwr Hyperparameter
mu2=0.001;%L1_upr Hyperparameter 
lamda=0.01;%Group-Lasso Hyperparameter
GLasso_en=1;
DisplayGraph=1;
IncidDisplay=1;
P=2; %Filter Order
T=2000; %Number of Samples
Tstep=1;
gamma=0.98;
delta=1-gamma;
HodgeNormlzn=1;

%Register Parameter
param.GLasso_en=GLasso_en;
param.DisplayGraph=DisplayGraph;
param.IncidDisplay=IncidDisplay;
param.gamma=gamma; %forgetting factor
param.T=T;
param.Tstep=Tstep;
param.P=P;
param.P=P;
param.mu1=mu1; 
param.mu2=mu2;  
param.lamda=lamda;
param.HodgeNormlzn=HodgeNormlzn;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose Experiment
% ExptNum=1: Fig.3a (Sioux Fall)
% ExptNum=1: Fig.3b (Sioux Fall)
% ExptNum=1: Fig.3c (Cherry Hills)
%%%%%%%%%%%%%%%%%%%%%%%%%%

ExptNum=3

switch ExptNum
    case 1
        ModelNum=1;
        param.ModelNum=ModelNum;
        SiouxFall    % Loading the Synthetic Data
    case 2
        ModelNum=2;
        param.ModelNum=ModelNum;
        SiouxFall    % Loading the Synthetic Data
    case 3
        CherryHills
    otherwise
        disp('Enter a Valid Experiment number (1,2, or 3)')
end


%%%%%RUN SVAR%%%%%%
KV=[0,1,2,3,4,5]; % Run for different graph filter orders.
for k=1:length(KV)
    param.mu1=0.001;
    param.mu2=0.01;
    K_lwr=KV(k);
    K_upr=KV(k);
    K=K_lwr+K_upr+1;
    param.K=K;
    param.K_lwr=K_lwr;
    param.K_upr=K_upr;
    param.lamda=0.001;
    [fp,nmse]=SVAR(f_time_series,Hodge,param); % MAIN FUNCTION FOR THE S-VAR ALGORITHM
    f_pred_buffer(:,:,k)=fp;
    nmse_buffer(:,:,k)=nmse;
end
nmse_SVAR=squeeze(mean(nmse_buffer,1));
[f_pred_TIRSO,nmse_TIRSO]=TIRSO_Predict(f_time_series,param); % Comparison: TIRSO Algorithm
[f_pr_MA,nmse_MA]=MA_Predict(f_time_series,param); % Comparison: Moving Average Algorithm 

DO_Plot



