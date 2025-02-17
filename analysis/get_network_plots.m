
graph_option='Strength'; % use 'Strength' for mean of reliable Js between groups; 'Density' for density of reliable Js between groups.
thr_rel=1.5; % threshold defining reliable couplings: J is reliable if |J|>thr_rel*errJ

dataset = fileread("param_dataset");
pne = fileread("param_pne");
dataset = dataset(1:end-1);
pne = pne(1:end-1);

use_cellist = ['../../datasets/',dataset,'/cellist.mat'];
use_dataset = ['../../datasets/',dataset,'/cellfiles'];

session_length = 8;
if strcmp(dataset,'ER')
    session_length = 9;
end

pos_neg_encoding=str2num(pne); % 1: separate positive and negative encoding cells; otherwise: pool positive and negative encoding cells.

node_names=1; % 1: add node names to the plot; otherwise: doesn't add node names
myMap=[0,0,1;1,0,0]; % blue, red: first color for negative (eSign<0) edges; second color for positive (eSign>0) edges
LineWidthRange=10;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(pos_neg_encoding)
    % 1=untuned; 2=offer value A positive; 3=offer value A negative; 4=offer value B positive; 5=offer value B negative; 6=chosen value positive; 
    % 7=chosen value negative; 8=chosen juice A; 9=chosen juice B
    nb_cell_types=9;
    pos_neg_string='_pn';
else
    % 1=untuned; 2=offer value A; 3=offer value B; 4=chosen value; 5=chosen juice A; 6= chosen juice B
    nb_cell_types=6;
    pos_neg_string='';
end
load(['../../datasets/',dataset,'/matlab_cache/par_groups',pos_neg_string],'J_dJ_groups_alls')
if(strcmp(graph_option,'Strength'))
    Jgraph=cell(1,nb_cell_types^2);
    for cc=1:nb_cell_types^2
        if(~isempty(J_dJ_groups_alls{1,cc}))
            ind=find(abs(J_dJ_groups_alls{1,cc}(:,1))>thr_rel*J_dJ_groups_alls{1,cc}(:,2));
            Jgraph{1,cc}=J_dJ_groups_alls{1,cc}(ind,1);
        end
    end
elseif(strcmp(graph_option,'Density'))
    Jgraph=cell(1,nb_cell_types^2);
    for cc=1:nb_cell_types^2
        if(~isempty(J_dJ_groups_alls{1,cc}))
            ind=double(abs(J_dJ_groups_alls{1,cc}(:,1))>thr_rel*J_dJ_groups_alls{1,cc}(:,2));
            Jgraph{1,cc}=sum(ind)/length(J_dJ_groups_alls{1,cc}(:,1));
        end
    end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean_Jgraph=zeros(1,nb_cell_types^2);
s=[];
t=[];
Weight=[];
eSign=[];
s_nocells=[];
t_nocells=[];
for ct1=1:nb_cell_types
    for ct2=ct1:nb_cell_types
        if(ct1==ct2)
            cc=ct1;
        else           
            cc=nb_cell_types+(nb_cell_types-1)*(ct1-1)+ct2-1;
        end
        mean_Jgraph(1,cc)=mean(Jgraph{1,cc});
        if(isnan(mean_Jgraph(1,cc)))
            s_nocells=[s_nocells;ct1];
            t_nocells=[t_nocells;ct2];
        elseif((mean_Jgraph(1,cc)~=0))
            s=[s;ct1];
            t=[t;ct2];
            Weight=[Weight;abs(mean_Jgraph(1,cc))];
            eSign=[eSign;sign(mean_Jgraph(1,cc))];
        end    
    end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Node properties:
if(pos_neg_encoding)    
    ds_pos=1/2;
    ds_neg=1/4;
    OVApxy=[0,1]-[0,ds_neg];
    OVAnxy=[0,1]+[0,ds_pos];
    OVBpxy=[3,1]-[0,ds_neg];
    OVBnxy=[3,1]+[0,ds_pos];
    CJAxy=[1,0];
    CJBxy=[2,0];
    CVpxy=[1.5,2]-[ds_neg,0];
    CVnxy=[1.5,2]+[ds_pos,0];
    NSxy=[1.5,1];
    % Node names:
    Name={'NS';'OVAp';'OVAn';'OVBp';'OVBn';'CVp';'CVn';'CJA';'CJB'}; % same order as in J_cellTypes_GLM_8.m
    % Node coordinates:
    nodexy_all=[NSxy;OVApxy;OVAnxy;OVBpxy;OVBnxy;CVpxy;CVnxy;CJAxy;CJBxy]; 
else   
    OVAxy=[0,1];
    OVBxy=[3,1];
    CJAxy=[1,0];
    CJBxy=[2,0];
    CVxy=[1.5,2];
    NSxy=[1.5,1];
    % Node names:
    Name={'NS';'OVA';'OVB';'CV';'CJA';'CJB'};  % same order as in J_cellTypes_GLM_8.m
    % Node coordinates:
    nodexy_all=[NSxy;OVAxy;OVBxy;CVxy;CJAxy;CJBxy]; 
end
nodex=nodexy_all(:,1);
nodey=nodexy_all(:,2);

% Edge properties and Graph G:
EndNodes=[s,t]; % EndNodes(i,:) = [source,target] nodes of edge i; Weight(i,1) = strength of edge i; eSign(i,1) = sign of edge i
EdgeTable=table(EndNodes,Weight,eSign);
NodeTable=table(Name);  
G= graph(EdgeTable,NodeTable);
LWidths = LineWidthRange*G.Edges.Weight/max(G.Edges.Weight);
EdgeColors = G.Edges.eSign;

% Edge properties and Graph G_nocells (indicating connections that are not sampled due to absence of cells in the groups):
if(~isempty(s_nocells))
    EndNodes=[s_nocells,t_nocells];
    Weight=ones(length(s_nocells),1); 
    EdgeTable=table(EndNodes,Weight);
    G_nocells= graph(EdgeTable,NodeTable);
    LWidths_nocells = LineWidthRange/5;
    EdgeColors_nocells = [0.2,0.2,0.2];
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure
hold on
colormap(myMap);
% plot all edges including self-loops
p=plot(G,'XData',nodex,'YData',nodey,'ShowArrows',false,'LineWidth',LWidths,'EdgeCData',EdgeColors,'LineStyle','-','EdgeAlpha',1,'Marker','o','MarkerSize',7,'NodeColor','k','NodeLabelColor','k');
if(~isempty(s_nocells))
    p_nocells=plot(G_nocells,'XData',nodex,'YData',nodey,'ShowArrows',false,'LineWidth',LWidths_nocells,'EdgeColor',EdgeColors_nocells,'LineStyle','--','EdgeAlpha',0.3,'Marker','o','MarkerSize',7,'NodeColor','k','NodeLabelColor','k');
end
if(node_names)
    p.NodeFontSize=15;
    p.NodeFontWeight='bold';
    p_nocells.NodeFontSize=15;
    p_nocells.NodeFontWeight='bold';
else
    p.NodeFontSize=0.01;
    p_nocells.NodeFontSize=0.01;
end
savefile = ['../../figures/new_network_plots/network_plot_1.png'];
saveas(gca,savefile);
hold off


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Wang's model

option_digraph=1;
option_arrows=0;

% Node properties for the Wang's model
OVAxy=[0,1];
OVBxy=[3,1];
CJAxy=[1,0];
CJBxy=[2,0];
CVxy=[1.5,2];
NSxy=[1.5,1];
% Node names:
Name={'NS';'OVA';'OVB';'CV';'CJA';'CJB'};  
% Node coordinates:
nodexy_all=[NSxy;OVAxy;OVBxy;CVxy;CJAxy;CJBxy]; 
nodex=nodexy_all(:,1);
nodey=nodexy_all(:,2);

% Edge properties and Graph G_Wang
% node code:
NS=1;
OVA=2;
OVB=3;
CV=4;
CJA=5;
CJB=6;
if(option_digraph==1)&&(option_arrows)
    Connected_CellsSign=[OVA,CJA,1;OVB,CJB,1;CJA,CJA,1;CJB,CJB,1;CJA,CJB,1;CJA,NS,1;CJB,NS,1;NS,NS,1;CJA,CV,1;CJB,CV,1;NS,CV,1;CV,CV,-1;CV,CJA,-1;CV,CJB,-1;CV,NS,-1];
elseif(option_digraph)
    Connected_CellsSign=[OVA,CJA,1;OVB,CJB,1;CJA,CJA,1;CJB,CJB,1;CJA,CJB,1;CJA,NS,1;CJB,NS,1;NS,NS,1;CJA,CV,1;CJB,CV,-1;NS,CV,1;CV,CV,-1;CV,CJA,-1;CV,CJB,1;CV,NS,-1]; % to better render the symmetry of the graph
else
    Connected_CellsSign=[OVA,CJA,1;OVB,CJB,1;CJA,CJA,1;CJB,CJB,1;CJA,CJB,1;CJA,NS,1;CJB,NS,1;NS,NS,1;CJA,CV,1;CJB,CV,1;NS,CV,1;CV,CV,-1];
end

s_Wang=Connected_CellsSign(:,1);
t_Wang=Connected_CellsSign(:,2);
eSign=Connected_CellsSign(:,3);
EndNodes=[s_Wang,t_Wang];
Weight=ones(length(s_Wang),1); 
EdgeTable=table(EndNodes,Weight,eSign);
NodeTable=table(Name);
if(option_digraph==1)
    G_Wang=digraph(EdgeTable,NodeTable);
    G_Wang_noselfloops=digraph(EdgeTable,NodeTable,'omitselfloops');
    EdgeColors_noselfloops = G_Wang_noselfloops.Edges.eSign;
else
    G_Wang=graph(EdgeTable,NodeTable);
end
EdgeColors = G_Wang.Edges.eSign;
LWidths = LineWidthRange/2;

figure
hold on
colormap(myMap);
p=plot(G_Wang,'XData',nodex,'YData',nodey,'ShowArrows',false,'LineWidth',LWidths,'EdgeCData',EdgeColors,'LineStyle','-','EdgeAlpha',1,'Marker','o','MarkerSize',7,'NodeColor','k','NodeLabelColor','k','NodeFontSize',15,'NodeFontWeight','bold');
if((option_digraph)&&(option_arrows))
    p_arrows=plot(G_Wang_noselfloops,'XData',nodex,'YData',nodey,'ShowArrows',true,'ArrowSize',35,'LineWidth',LWidths,'EdgeCData',EdgeColors_noselfloops,'LineStyle','-','EdgeAlpha',1,'Marker','o','MarkerSize',7,'NodeColor','k','NodeLabelColor','k','NodeFontSize',15,'NodeFontWeight','bold');
end
if(node_names)
    p.NodeFontSize=15;
    p.NodeFontWeight='bold';
    if((option_digraph)&&(option_arrows))
        p_arrows.NodeFontSize=15;
        p_arrows.NodeFontWeight='bold';
    end
else
    p.NodeFontSize=0.01;
    if((option_digraph)&&(option_arrows))
        p_arrows.NodeFontSize=0.01;
    end    
end
savefile = ['../../figures/new_network_plots/network_plot_2.png'];
saveas(gca,savefile);
hold off








