function [Beta,PHI,r,eta]=COMID(Beta,Ft,ft,Hodge,PHI,r,param)
mu1=param.mu1; %L1_lwr Hyperparameter
mu2=param.mu2; %L1_upr Hyperparameter  
lamda=param.lamda; %Group-Lasso Hyperparameter
gamma=param.gamma; %forgetting factor
GLasso_en=param.GLasso_en;
P=param.P;
tp=(param.t);

L1_upr=Hodge.L1_upr;
L1_lwr_jn=Hodge.L1_lwr_jn;




K_lwr=param.K_lwr;
K_upr=param.K_upr;


delta=1-gamma;
K=param.K;
NumEdg=size(L1_lwr_jn,1);
PHI=gamma*PHI+delta*Ft'*(eye(NumEdg)+mu1*L1_lwr_jn+mu2*L1_upr)*Ft;
eta=1/(max(eigs(PHI))*tp);
r=gamma*r+delta*Ft'*ft;
Grad=PHI*Beta-r;

if GLasso_en==1
   for p=1:P
      nL=(p-1)*K+1:(p-1)*K+1+K_lwr;
      nU=(p-1)*K+K_lwr+2:p*K;
   
        
        Beta_L=Beta(nL);
        Grad_L=Grad(nL);
        GD_update_L=Beta_L-eta*Grad_L;
        GDUnormL=norm(GD_update_L);
        Beta_L=GD_update_L*max(0,1-eta*lamda/GDUnormL);
        
        Beta_U=Beta(nU);
        Grad_U=Grad(nU);
        GD_update_U=Beta_U-eta*Grad_U;
        GDUnormU=norm(GD_update_U);
        Beta_U=GD_update_U*max(0,1-eta*lamda/GDUnormU);
        
        Beta(nL)=Beta_L;
        Beta(nU)=Beta_U;
   end

    Beta_bias=Beta(end);
    Grad_bias=Grad(end);
    GD_update_bias=Beta_bias-eta*Grad_bias;
    Beta(end)=GD_update_bias;
    
else
    GD_update=Beta-eta*Grad;
    GDUnorm=norm(GD_update);
    Beta=GD_update*max(0,1-eta*lamda/GDUnorm);
end

end