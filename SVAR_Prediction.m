 function Tstep_prediction=SVAR_Prediction(Beta,Ft,ft,Hodge,PHI,r,param)
K_lwr=param.K_lwr;
K_upr=param.K_upr;
K=param.K;
NumEdg=param.NumEdg;
f_ONE=ones(NumEdg,1);
Tstep=param.Tstep;
L1_lwr=Hodge.L1_lwr;
L1_upr=Hodge.L1_upr;
 
 for n=1:Tstep
    Fs_lwr=Collect_Shifted_Mx(ft,L1_lwr,K_lwr,'lwr');
    Fs_upr=Collect_Shifted_Mx(ft,L1_upr,K_upr,'upr');
  
    Fs =[Fs_lwr,Fs_upr]; 
    F_tp1=[Fs,Ft(:,1:end-K-1)];
    F_tp1=[F_tp1,f_ONE];
    
    Tstep_prediction(:,n)=F_tp1*Beta;
    ft=Tstep_prediction(:,n);
    [Beta,PHI,r,eta]=COMID(Beta,F_tp1,ft,Hodge,PHI,r,param);
 end
 end