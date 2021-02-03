function [modularity_table] = create_mod_figure(SC_mats, Ci_max, conn_hubs, prov_hubs, outdir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB script to create figure of module configuration and connector hub
% connectivity and get the modularity table with the regions belonging to 
% each module
% 
% INPUTS:
%   - SC_mats:          3D weighted structural connectivity matrices      
% 
%   - Ci_max:           modularity partiotioning. Vector of module-codes
%                       assigned to each area of the parcellation.
%                         
%   - conn_hubs:        vector of connector hubs.
%
%   - prov_hubs:        vector of provincial hubs.
%
%   - outdir:           output directory to which the figure will be saved.
%
% OUTPUTS:
%   - modularity_table: table with identification of the regions belonging
%                       to each module
%
%
% USAGE: 
%       ex: modularity_analysis(subjects_PBD.SC_3D,modules_HC.Ci,...
%               modules_PBD.cHubs,'modules')
%
%
% Ana Coelho (01/2021)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ci_max = m3_mod_results.results.Ci;
conn_hubs = m3_mod_results.results.conn_hubs;
prov_hubs = m3_mod_results.results.prov_hubs;
SC_mats = m3_mats.m3_mats_cons_thr;

n_aal  = size(SC_mats,2);
n_pop1 = 1:size(SC_mats,3);

SC_avg = zeros(size(SC_mats,2));

for n=1:n_aal
    for p=1:n_aal
        SC_avg(n,p) = mean(squeeze(SC_mats(n,p,n_pop1)));
    end
end

SC_norm = SC_avg/max(abs(SC_avg(:)));

% Loading Labels and COGs
sourcedir ='/Users/anacoelho/Documents/MATLAB/';
load([sourcedir 'utils/parcelations/silhouette_thr300_cog.txt']);
load([sourcedir 'utils/parcelations/silhouette_thr300_lat_labels.mat']);
label90 = silhouette_lat_labels;
MNI_coord = silhouette_thr300_cog;
clear silhouette_thr300_cog
n_aal  = size(SC_norm,1);

addpath(genpath([sourcedir 'utils']));

hubsFig = figure('units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto','color','white');
sub(1)  = subplot(1,2,1); hold on
sub(2)  = subplot(1,2,2); hold on
colormap HSV
cm_nodes = colormap;
% Initialize modules structure
modules.mod(1).counter=0; modules.mod(1).areas=[];
modules.mod(2).counter=0; modules.mod(2).areas=[];
modules.mod(3).counter=0; modules.mod(3).areas=[];
modules.mod(4).counter=0; modules.mod(4).areas=[];
modules.mod(5).counter=0; modules.mod(5).areas=[];
modules.mod(6).counter=0; modules.mod(6).areas=[];
modules.mod(7).counter=0; modules.mod(7).areas=[];
modules.mod(8).counter=0; modules.mod(8).areas=[];
modules.mod(9).counter=0; modules.mod(9).areas=[];
modules.mod(10).counter=0; modules.mod(10).areas=[];

for n=1:n_aal
    switch Ci_max(n)
        case 1
            %             colour=[1 1 0];     % yellow
            colour=cm_nodes(49,:);
            modules.mod(1).counter=modules.mod(1).counter+1;
            modules.mod(1).areas=[modules.mod(1).areas n];
        case 2
            %             colour=[1 0 1];     % magenta
            colour=cm_nodes(7,:);
            modules.mod(2).counter=modules.mod(2).counter+1;
            modules.mod(2).areas=[modules.mod(2).areas n];
        case 3
            %             colour=[0 1 1];     % cyan
            colour=cm_nodes(37,:);
            modules.mod(3).counter=modules.mod(3).counter+1;
            modules.mod(3).areas=[modules.mod(3).areas n];
        case 4
            %             colour=[1 0 0];     % red
            colour=cm_nodes(19,:);
            modules.mod(4).counter=modules.mod(4).counter+1;
            modules.mod(4).areas=[modules.mod(4).areas n];
        case 5
            %             colour=[0 1 0];     % green
            colour=cm_nodes(43,:);
            modules.mod(5).counter=modules.mod(5).counter+1;
            modules.mod(5).areas=[modules.mod(5).areas n];
        case 6
            %           colour=[1 1 1];     % white
            colour=cm_nodes(55,:);
            modules.mod(6).counter=modules.mod(6).counter+1;
            modules.mod(6).areas=[modules.mod(6).areas n];
        case 7
            %             colour=[0 0 1];     % blue
            colour=cm_nodes(13,:);
            modules.mod(7).counter=modules.mod(7).counter+1;
            modules.mod(7).areas=[modules.mod(7).areas n];
        case 8
            %             colour=[0 0 0];     % black
            colour=cm_nodes(1,:);
            modules.mod(8).counter=modules.mod(8).counter+1;
            modules.mod(8).areas=[modules.mod(8).areas n];
        case 9
            %             colour=[0 0 0];     % black
            colour=[1, 1, 0.16];
            modules.mod(9).counter=modules.mod(9).counter+1;
            modules.mod(9).areas=[modules.mod(9).areas n];
        case 10
            colour=[1, 0.898, 0.84];
            modules.mod(10).counter=modules.mod(10).counter+1;
            modules.mod(10).areas=[modules.mod(10).areas n];
       case 11
            colour=[0.61, 0.59, 1];
            modules.mod(11).counter=modules.mod(11).counter+1;
            modules.mod(11).areas=[modules.mod(11).areas n];
       case 12
            colour=[0.60, 1, 0.878];
            modules.mod(12).counter=modules.mod(12).counter+1;
            modules.mod(12).areas=[modules.mod(12).areas n];
       case 13
            colour=[0.56, 0.56, 0.56];
            modules.mod(13).counter=modules.mod(13).counter+1;
            modules.mod(13).areas=[modules.mod(13).areas n];
       case 14
            colour=[0.2, 0.56, 0.37];
            modules.mod(14).counter=modules.mod(14).counter+1;
            modules.mod(14).areas=[modules.mod(14).areas n];
       case 15
            colour=[0.5, 0.36, 0.25];
            modules.mod(15).counter=modules.mod(15).counter+1;
            modules.mod(15).areas=[modules.mod(15).areas n];
        otherwise
            error('SWITCH ONLY PREPARED FOR MAX OF 6 MODULES');
    end
    for p = n+1:n_aal
        % EFFICIENCY
        % Inter-modularity - to/from connector hubs only
        if (ismember(n, conn_hubs) || ismember(p, conn_hubs)) && Ci_max(n)~=Ci_max(p) && SC_norm(n,p)>0 % intermodularity
            [X,Y,Z]=cylinder1(MNI_coord(n,:),MNI_coord(p,:),3*SC_norm(n,p),100);
            surf(X,Y,Z,'FaceColor',[0 0 0],'EdgeColor','none','parent',sub(2));
        end
        % Intra-modularity
        if Ci_max(p)==Ci_max(n)
            [X,Y,Z]=cylinder1(MNI_coord(n,:),MNI_coord(p,:),0.01,100);
            surf(X,Y,Z,'FaceColor',colour,'EdgeColor',colour,'parent',sub(1));
        end
    end
    % PLOT NODES ----------------------------------
    % NODES - MODULARITY
    sphere2(MNI_coord(n,1),MNI_coord(n,2),MNI_coord(n,3),1,colour,0,sub(2));
    %     text(MNI_coord(n,1),MNI_coord(n,2)+3,MNI_coord(n,3)+3,num2str(n),'parent',sub(1),'FontWeight','bold','FontSize',7);
    sphere2(MNI_coord(n,1),MNI_coord(n,2),MNI_coord(n,3),1,colour,0,sub(1));
    % Create 3D discs around nodes to identify both connector and provincial HUBS
    if ismember(n, conn_hubs)
        %create3Ddisc(MNI_coord(n,:), 0.9, 2, [0.5 0.5 0.5], [0.5 0.5 0.5], sub(1));
        %create3Ddisc(MNI_coord(n,:), 0.9, 2, [0 0 0], [0 0 0], sub(2));
        create3Ddisc(MNI_coord(n,:), 2, 4, [0.5 0.5 0.5], sub(2));
        create3Ddisc(MNI_coord(n,:), 2, 4, [0 0 0], sub(1));
    elseif ismember(n, prov_hubs)
        %create3Ddisc(MNI_coord(n,:), 0.3, 2.3, [0.5 0.5 0.5], [0.5 0.5 0.5], sub(1));
        %create3Ddisc(MNI_coord(n,:), 0.3, 2.3, [0 0 0], [0 0 0], sub(2));
        create3Ddisc(MNI_coord(n,:), 2, 2, [0.5 0.5 0.5], sub(2));
        create3Ddisc(MNI_coord(n,:), 2, 2, [0 0 0], sub(1));
    end
    clear colour;
end
axis(sub, 'off')
axis(sub, 'equal')

modularity_xlslabels{1} = {'Module ID','Number of nodes','Node members','Node members'};
modularity_xlslabels{2} = {'','','number','label'};

modularity_table        = [modularity_xlslabels{1}; modularity_xlslabels{2}];
modules.number     = max(Ci_max);

for g=1:modules.number 
    for h=1:length(modules.mod(g).areas)
        if h==1
            modularity_row   = [num2cell(g), modules.mod(g).counter, num2cell(modules.mod(g).areas(h)), label90(modules.mod(g).areas(h),:)];
            modularity_table = [modularity_table; modularity_row];
        else
            modularity_row   = ['.', '.', num2cell(modules.mod(g).areas(h)), label90(modules.mod(g).areas(h),:)];
            modularity_table = [modularity_table; modularity_row];        
        end
    end
    
end 

% save figure
set(0, 'CurrentFigure', hubsFig);
print(strcat(out_dir,'/modularity.png'),'-dpng','-r300');  

end