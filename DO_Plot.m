%%%%%%%PLOTS%%%%%%%%%%%
x=1:T;
pix=600;
figure('Position',[1000,-1000,pix,500])
clr=[1 0 1; 0 0 0; 0.6350 0.0780 0.1840; 0 0 1;0.3010 0.7450 0.9330; 1 0 0; 0 1 0;0.9290 0.6940 0.1250];

set(gca, 'ColorOrder',clr , 'NextPlot', 'replacechildren');
plt=plot(x,nmse_SVAR,x,(mean(nmse_TIRSO',2)),x,(mean(nmse_MA',2)),'MarkerIndices', 1:100:length(x));
set(plt,{'Marker'},{'^';'o';'<';'diamond';'*';'+';'o';'square'})
set(gca,'YScale','log');
set(plt,'LineWidth',2)

grid on
box on
for i=1: size(clr,1)
plt(i).MarkerFaceColor=clr(i,:);
plt(i).MarkerSize=7;
end

switch ExptNum
   case 1
        YLim=2;
    case 2
        YLim=5;
    case 3
        YLim=0.25;
    otherwise
        disp('Enter a Valid Experiment number (1,2, or 3)')
end
ylim([0,YLim])


legend('S-VAR (K=1)','S-VAR (K=3)','S-VAR (K=5)','S-VAR (K=7)','S-VAR (K=9)','S-VAR (K=11)','TIRSO','MA')
xlabel('T')
ylabel('NMSE (log scale)')
set(gca,'fontsize',25)

