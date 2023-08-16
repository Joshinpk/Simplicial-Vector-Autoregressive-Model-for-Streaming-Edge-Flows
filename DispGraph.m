
function DispGraph(edg,f_val,ttag)
% f_test_str="f"+f_val;
f_test_str=f_val;




figure
G = digraph(edg(:,1),edg(:,2),f_val);
    p0=plot(G,'Layout', 'force','Marker', 'o', 'EdgeAlpha', 1,'EdgeLabel',f_test_str,...
        'EdgeFontSize',10,'MarkerSize', 10,'LineWidth',2,'NodeFontSize',15,...
        'NodeColor','black','NodeLabelColor','k','ArrowSize',15,'EdgeCData',f_val,...
        'Marker', 'o','UseGravity',true);


% Y=[4,0,0,2,2,0];
% X=[2,0,4,3,1,2];
% G = digraph(edg(:,1),edg(:,2),f_val);
%     p0=plot(G,'XData',X,'YData',Y,'Marker', 'o', 'EdgeAlpha', 1,'EdgeLabel',f_test_str,...
%         'EdgeFontSize',15,'MarkerSize', 12,'LineWidth',2,'NodeFontSize',12,...
%         'NodeColor','black','NodeLabelColor',[0 0 0],'ArrowSize',15,'EdgeCData',f_val,...
%         'Marker', 'O');
%     w=f_val-min(f_val);
%     p0.LineWidth = 4*w/max(w)+1;%G.Edges.LWidths;
% %     p0.ArrowSize=3*p0.LineWidth
%     title(ttag,'FontSize', 15,'Interpreter','latex','position', [0.5, 3.5])
% set(gca,'FontSize', 20)
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% mymap = [.1 .1 0
%     1 0 0
%     0 0 1
%     0 0 1];
% colormap(mymap)
% box off
% axis off

end
% 'auto', 'circle', 'force', 'layered', 'subspace', 'subspace3', 'force3'