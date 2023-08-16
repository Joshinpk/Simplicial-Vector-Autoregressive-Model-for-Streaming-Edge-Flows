
function [f_pred_time_series,nmse_SVAR]=SAVR(f_time_series,Hodge,param)

P=param.P;
K_lwr=param.K_lwr;
K_upr=param.K_upr;
K=param.K;
mu1=param.mu1; 
mu2=param.mu2;  
lamda=param.lamda;
GLasso_en=param.GLasso_en;
DisplayGraph=param.DisplayGraph;
IncidDisplay=param.IncidDisplay;
gamma=param.gamma;
T=param.T;
Tstep=param.Tstep;
L1_lwr=Hodge.L1_lwr;
L1_upr=Hodge.L1_upr;

NumEdg=size(f_time_series,1);
param.NumEdg=NumEdg;
KP=K*P+1;
PHI=1*eye(KP);
r=0*rand(KP,1);

Beta=0.00*ones(KP,1);
BV=[];
BR=[];
etaV=[];
f_pred_time_series=zeros(size(f_time_series,1),size(f_time_series,2)+1);
f_ONE=ones(NumEdg,1);
for t=P+1:T
    ft=f_time_series(:,t);
    Ft=[];
    for p=1:P
        f_t_minus_p=f_time_series(:,t-p);
        Fs_lwr=Collect_Shifted_Mx(f_t_minus_p,L1_lwr,K_lwr,'lwr');
        Fs_upr=Collect_Shifted_Mx(f_t_minus_p,L1_upr,K_upr,'upr');
        Fs =[Fs_lwr,Fs_upr];
        Ft=[Ft,Fs];
        test1(t,:)=Fs_lwr(2,:);
    end
     Ft=[Ft,f_ONE];
     param.t=1;
     
    [Beta,PHI,r,eta]=COMID(Beta,Ft,ft,Hodge,PHI,r,param);
  
    etaV=[etaV;eta];
    BV=[BV;Beta'];
    BR=[BR;r'];
    Tstep_prediction=SVAR_Prediction(Beta,Ft,ft,Hodge,PHI,r,param);
    f_pred_time_series(:,t+Tstep)=Tstep_prediction(:,Tstep);
end
nmse_SVAR=CompNMSE(f_time_series,f_pred_time_series(:,1:end-Tstep));
end




