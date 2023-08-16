 %%%%%TEST-2:Cherry Hills%%%%%%%%%%%
div_free=0;
IncidDisplay=param.IncidDisplay;
DisplayGraph=param.DisplayGraph;
T=param.T;
P=param.P;
HodgeNormlzn=param.HodgeNormlzn;

directed_links=[1,2;2,3;2,4;4,5;5,3;...
                    3,6;6,7;7,8;7,10;8,9;10,11;11,12;12,13;...
                    13,15;15,20;13,14;19,20;14,19;18,19;17,18;14,16;...
                    21,32;32,33;33,34;32,34;...
                    21,22;20,21;22,23;23,24;23,25;25,26;26,27;...
                    27,30;34,35;35,36;...
                    16,17;27,28;28,29;28,30;30,31];
edg_index_unsorted=1:40;                
    
[edg,sort_idx]=sortrows(directed_links); % Lexigocraphic sorting
edg_index_sorted=edg_index_unsorted(sort_idx);

load("CherryHillsData")
f_time_series=CherryHillsData(:,1:T);
            
if DisplayGraph==1           
    DispGraph(edg,edg_index_sorted,"Cherry Hills")
end

[B1,B2]=Compute_Incidence(edg,param);

End_nodes=[1,9,24,29,31,36];


B1(End_nodes,:)=0; % Removing end nodes from flow conservation
L1_lwr=B1'*B1; %Lower Hodge Laplacian (Order-1)
L1_upr=B2*B2'; %Upper Hodge Laplacian (Order-1)


B1_jun=B1;
B1_jun(End_nodes,:)=0;
L1_lwr_jn=B1_jun'*B1_jun;

Hodge.L1_lwr=L1_lwr;
Hodge.L1_upr=L1_upr;
Hodge.L1_lwr_jn=L1_lwr_jn;


if HodgeNormlzn==1
    [Ul,Sl,Vl] = svd(L1_lwr);
    Sl_norm=Sl/(max(Sl(:)));
    L1_lwr=Ul*Sl_norm*Vl';

    [Uu,Su,Vu] = svd(L1_upr);
    Su_norm=Su/(max(Su(:)));
    L1_upr=Uu*Su_norm*Vu';
end















