function [m_prediction,nmse_TIRSO]=TIRSO_Predict(f_time_series,param)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIRSO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n_row_mx,nTimeInstants]=size(f_time_series);
tirsoObj = Tirso; % set tirso object up
tirsoObj.noOfNodes = n_row_mx;
tirsoObj.order     = param.P; % we can try a higher order later
tirsoObj.regPar    = 1e-1;
tirsoObj.b_shrinkSelfLoops  = 0; % Bolstad
tirsoObj.forgettingFactor   = param.gamma;
tirsoObj.h_stepsize         = @(ts)1/eigs(ts.m_Phi,1);
Tstep=param.Tstep;
% initialize
tState_in = tirsoObj.initialize(0, f_time_series( :,1:tirsoObj.order)');
for t = tirsoObj.order+1:nTimeInstants
mtemp= f_time_series( :,t);
    tState_in = tirsoObj.update(tState_in, mtemp);
    m_predic(:,:)=tState_in.predictManyFromBuffer(Tstep)';%10 step ahead prediction
     m_prediction(1: n_row_mx,t+Tstep)= m_predic(:,Tstep); %store 1 step ahead prediction
    %CHECK!!!!!!!!!!!
    %%%%%
end
%%
nmse_TIRSO=CompNMSE(f_time_series,m_prediction(:,1:end-Tstep));


end