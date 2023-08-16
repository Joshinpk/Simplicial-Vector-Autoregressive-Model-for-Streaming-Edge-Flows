
directed_links=[1,2;2,6;3,1;3,4;4,5;5,9;6,5;6,8;7,8;8,9;8,16;9,10;10,16;10,17;...
                11,4;11,10;12,3;12,11;13,12;14,11;14,23;15,10;15,14;16,18;17,16;...
                18,7;19,15;19,17;20,18;20,19;21,20;21,22;22,15;22,20;23,22;23,24;24,13;24,21];
f_tag=1:38;

[edg,sort_idx]=sortrows(directed_links); %Lexigocraphic sorting
f_tag_s=f_tag(sort_idx);
if DisplayGraph==1    
    DispGraph(edg,f_tag_s,"Sioux Falls transportation network")
end
NumEdge=size(edg,1);

f=0.01*randn(NumEdge,P); %Intialization for first P samples

%B1: node to edge incidence matrix, B2: edge to triangle incidence matrix. 
[B1,B2]=Compute_Incidence(edg,param);

L1_lwr=B1'*B1; %Lower Hodge Laplacian (Order-1)
L1_upr=B2*B2'; %Upper Hodge Laplacian (Order-1)

sigma=1.0;%noise power
for p=1:P
    A(:,:,p)=VAR_Coeff_Gen(ModelNum,NumEdge,B1,B2);
end
for t=P+1:T
    for p=1:P
        if mod(t,100)==0
            A(:,:,p)=VAR_Coeff_Gen(ModelNum,NumEdge,B1,B2);
        end
        f(:,t)=A(:,:,p)*f(:,t-p)+sigma*randn(NumEdge,1);% check
    end
end

f_time_series=f;

Hodge.L1_lwr=L1_lwr;
Hodge.L1_upr=L1_upr;
L1_lwr_jn=L1_lwr;
Hodge.L1_lwr_jn=L1_lwr_jn;

function A=VAR_Coeff_Gen(ModelNum,NumEdge,B1,B2)
    if ModelNum==1
       C=inv(rand(1)*eye(NumEdge)+rand(1)*B1'*B1+rand(1)*B2*B2'+rand(1)*(B1'*B1)^2+rand(1)*(B2*B2')^2);
    elseif ModelNum==2
         C=rand(NumEdge,NumEdge);
    else
       error("Choosel ModelNum as 1 or 2") 
    end
        
    B=C+C';
    [U,S]=eig(B);
    S=real((1*S/max(S(:))));
    dia=diag(S);
    dia(end-15:end)=0.95;
    S=diag(dia);
    A=U*S*U';
    A=A;
end
