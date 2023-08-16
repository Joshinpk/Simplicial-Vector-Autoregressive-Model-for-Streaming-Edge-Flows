 function [f_pr,nmse_MA]=MA_Predict(ft,param)
 P=param.P;
 T=param.T
 for t=P+1:T
     f_pr(:,t+1)=mean(ft(:,t-P+1:t),2);
 end
 nmse_MA=CompNMSE(ft,f_pr(:,1:end-1));
 end