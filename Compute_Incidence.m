function [B1,B2]=Compute_Incidence(edg,param)

G = digraph(edg(:,1),edg(:,2));

B1=full(G.incidence); % Node to Edge Incidence Matrix

%Find the 2-simplices (triangles)
NumEdge=length(edg(:,1));
AllNodes=(1:max(edg(:)))';
N=length(AllNodes);

two_simplexes=[];
for i=1:NumEdge
    edg_sel=edg(i,:);
    edg_pool_i=setdiff (edg, edg_sel, 'rows');
    edg_pool_i=[edg_pool_i;flip(edg_pool_i,2)];
    node_1=edg_sel(1);
    node_2=edg_sel(2);
    possible_edges_i=[node_1*ones(N,1),AllNodes;node_2*ones(N,1),AllNodes];
    edge_detection=ismember(edg_pool_i,possible_edges_i,'rows');
    edge_detected_idx=find(edge_detection==1);
    edge_detected=edg_pool_i(edge_detected_idx,:);
    for j=1:size(edge_detected,1)
        MyEdge=edge_detected(j,:);
        OtherNode=setdiff(edg_sel,MyEdge(1,1));
        OtherEdge=[OtherNode,MyEdge(1,2)];
        triangle_flag=nnz(ismember(edg_pool_i,OtherEdge,'rows'));
        if triangle_flag==1
            triangle=union(MyEdge,edg_sel);
            two_simplexes=[two_simplexes;triangle];
        end
    end
end
two_simplexes=unique(two_simplexes,'row');


% B2 Calculation
e_flip=flip(edg,2);
B2=[];
for k=1:size(two_simplexes,1)
tk=two_simplexes(k,:);
tk=[tk,tk(1)];
tk_e=[tk(1:2);tk(2:3);tk(3:4)];
incid_plus_zero_k= sum(abs(edg-repmat(tk_e(1,:),size(edg,1),1)),2).*...
              sum(abs(edg-repmat(tk_e(2,:),size(edg,1),1)),2).*...
              sum(abs(edg-repmat(tk_e(3,:),size(edg,1),1)),2);
incid_plus_k=ismember(incid_plus_zero_k,0);
incid_minus_zero_k= sum(abs(e_flip-repmat(tk_e(1,:),size(edg,1),1)),2).*...
               sum(abs(e_flip-repmat(tk_e(2,:),size(edg,1),1)),2).*...
               sum(abs(e_flip-repmat(tk_e(3,:),size(edg,1),1)),2);
incid_minus_k=ismember(incid_minus_zero_k,0);
incid_k=incid_plus_k-incid_minus_k;
B2=[B2, incid_k];
end



%Print Incidence Matrices in Table form.
if(param.IncidDisplay==1)
    B1_r="n"+ (1:N)';
    B1_c="e"+ (1:NumEdge)';
    B1Table = array2table(B1,'rowNames',B1_r,'VariableNames',B1_c)

    T=size(B2,2);
    B2_c1=strcat('(',num2str(two_simplexes(:,1)),',',num2str(two_simplexes(:,2)),',',num2str(two_simplexes(:,3)),')');
    B2_c=cellstr(B2_c1);
    B2_r="e"+ (1:NumEdge)';
    % B2_c="t"+ (1:T)';
    B2Table = array2table(B2,'rowNames',B2_r,'VariableNames',B2_c)
end

end