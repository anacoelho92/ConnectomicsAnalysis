function [results] = modularity_and_hubs(SC_mats, niterLouvain, out_dir)

% Function to find optimal community structure, by finding the one with
% higher number of occurences
% Since the same modularity structure can appear with modules in a
% different order, it's necessary to perform a match of equal structures
% and then find unique community structures and count how many times each
% one occurs
% The structure with higher number of counts with the one chosen as optimal
%
% Then it computes provincial and connector hubs and creates images of
% modularity structure and connector hubs connectivity profile
%
% Inputs:   SC_mats: 3D weighted connectivity matrices
%
%           niterLouvain: number of iterations to run community Louvain
%           algorithm
%
%           out_dir: output directory to save figures and results
%
% Outputs:  results: structure with the following outpus:
%
%           -   Ci: Affiliation vector of optimal community structure
%
%           -   patterns: IDs, number of occurences and mean(Q) of each unique structure
%
%           -   Q_mean: mean value (for all occurences of the optimal structure)
%           of the modularity statistic
%
%           -   conn_hubs: list of connector hubs identified (ids)
%           -   connectorHubs_lst: list of connector hubs identified (labels)
%
%           -   prov_hubs: list of provincial hubs identified (ids)
%           -   provincialHubs_lst: list of provincial hubs identified (labels)
%
%
% Henrique Fernandes & Ana Coelho 2018
% (last update 04/2020)
%
% email: henrique.fernandes@clin.au.dk
%


% convert subjects connectivity matrices to normalized group matrix
n_aal  = size(SC_mats,2);
n_pop1 = 1:size(SC_mats,3);

SC_avg = zeros(size(SC_mats,2));

for n=1:n_aal
    for p=1:n_aal
        SC_avg(n,p) = mean(squeeze(SC_mats(n,p,n_pop1)));
    end
end

SC_norm = SC_avg/max(abs(SC_avg(:)));


% Louvain algorithm to compute modularity coefficient

% struct to save results
results = struct();

% 3D matrix to save affiliation vector at each iteration (maximum number of modules 20)
Ci_bin_all=zeros(20,size(SC_norm,1),niterLouvain);

% array to save modularity statistic at each iteration
Q_all = zeros(niterLouvain,1);

for i=1:niterLouvain
    [Ci_arr, Q_arr] = community_louvain(SC_norm);
%     fprintf('Iteration %d \n',i);
    % binary matrix of community structure (each row is a module, each
    % column is a brain region)
    Ci_bin = zeros(20,size(Ci_arr,1));
    
    for j=1:size(Ci_arr,1)
        Ci_bin(Ci_arr(j),j)=1;
    end
    
    % save community structure of this iteration to the 3D matrix
    Ci_bin_all(:,:,i) = Ci_bin;
    % save modularity coefficient of this iteration to array
    Q_all(i,1) = Q_arr;
    
end


% get number of modules of each partition
n_mod = zeros(niterLouvain,1);
for i=1:size(Ci_bin_all,3)
    n_mod(i,1)=size(Ci_bin_all,1)-length(find(all(Ci_bin_all(:,:,i)==0,2)));
end

% find unique community structures and count occurences
Ci_unique = zeros(niterLouvain,1);
counts_unique = zeros(niterLouvain,1);
Qmean_unique = zeros(niterLouvain,1);


ids=1:size(Ci_bin_all,3);

while ~isempty(ids)
    i=ids(1);
    comm = [];
    % compare current structure with the others
    for j=1:size(Ci_bin_all,3)
        if j~=i
            
            % remove rows with zeros (non-existent modules)
            Ci1 = Ci_bin_all(:,:,i);
            rows_zeros = find(all(Ci1==0,2));
            for k=length(rows_zeros):-1:1
                Ci1(rows_zeros(k),:)=[];
            end
            
            % remove rows with zeros (non-existent modules)
            Ci2 = Ci_bin_all(:,:,j);
            rows_zeros = find(all(Ci2==0,2));
            for k=length(rows_zeros):-1:1
                Ci2(rows_zeros(k),:)=[];
            end
            
            if size(Ci2,1) == size(Ci1,1)
                equal = check_structures_equal(Ci1,Ci2,n_mod(i,1));
            else
                equal = 0;
            end
            
            % if structures are equal save the index
            if equal
                comm = [comm, j];
            end
        end
    end
    
    % remove indices that were found to be equal
    ids_copy = ids;
    comm = cat(2,i, comm);
    ids = setdiff(ids_copy,comm);
    
    % add counts of equal structure and compute mean modularity coefficient
    Ci_unique(i,1) = i;
    counts_unique(i,1) = length(comm);
    for k=1:length(comm)
        Qmean_unique(i,1) = Qmean_unique(i,1) + Q_all(comm(k),1);
    end
    Qmean_unique(i,1) = Qmean_unique(i,1) / counts_unique(i,1);
    disp(numel(ids));
end

% optimal structure is the one with higher number of occurences
Ci_bin_max = Ci_bin_all(:,:,find(counts_unique==max(counts_unique)));
rows_zeros = find(all(Ci_bin_max==0,2));
for k=length(rows_zeros):-1:1
    Ci_bin_max(rows_zeros(k),:)=[];
end

% convert to vector
for i=1:size(Ci_bin_max,1)
    Ci_bin_max(i,:) = Ci_bin_max(i,:)*i;
end

Ci_max = transpose(sum(Ci_bin_max,1));
results.Ci = Ci_max;

patterns = table(find(counts_unique),counts_unique(find(counts_unique)),Qmean_unique(find(counts_unique)),'VariableNames',{'ID','Counts','Qmean'});
results.patterns = patterns;

%counts_max = max(counts_unique);
%results.counts = counts_max;

Q_max = Qmean_unique(find(counts_unique==max(counts_unique)));
results.Q_mean = Q_max;


modules.number     = max(Ci_max);
modules.modularity = Q_max;

%%%%%% PROVINCIAL AND CONNECTOR HUBS %%%%%%

% Loading Labels and COGs
sourcedir ='/Users/anacoelho/Documents/MATLAB/';
load([sourcedir 'utils/parcelations/silhouette_thr300_cog.txt']);
load([sourcedir 'utils/parcelations/silhouette_thr300_lat_labels.mat']);
label90 = silhouette_lat_labels;
MNI_coord = silhouette_thr300_cog;
clear silhouette_thr300_cog
n_aal  = size(SC_norm,1);

addpath(genpath([sourcedir 'utils']));

% participation coefficient threshold
p_thr=0.3;

% Measure diversity of intermodular connections of individual nodes
mod_deg_zscore_vec = module_degree_zscore(SC_norm, Ci_max);

mod_zcore_avg   = mean(mod_deg_zscore_vec);
mod_zcore_stdev = std2(mod_deg_zscore_vec);
mod_zscore_thr  = mod_zcore_avg+mod_zcore_stdev;

% vector with the indexes of the nodes that show relevant intramodular connectivity - hubs (connector hubs + provincial hubs)
intra_hubs = find(mod_deg_zscore_vec>mod_zscore_thr);

% Measure of diversity of intermodular conditions of individual nodes
p_coef_vec = participation_coef(SC_norm, Ci_max);

% vector with the indexes of the nodes that show relevant intermodular connectivity - potential connector hubs
inter_hubs = find(p_coef_vec>p_thr);

% Defining the connector and provincial hubs - AAL areas number
conn_hubs  = intersect(intra_hubs, inter_hubs);

prov_hubs  = (~ismember(intra_hubs, inter_hubs).*intra_hubs);
prov_hubs(prov_hubs==0) = [];

% List of connector Hubs:
connectorHub_lst  = label90(conn_hubs,:);
provincialHub_lst = label90(prov_hubs,:);

results.conn_hubs = conn_hubs;
results.connectorHubs_lst = connectorHub_lst;

results.prov_hubs = prov_hubs;
results.provincialHubs_lst = provincialHub_lst;

results.mod_deg_zscore = mod_deg_zscore_vec;
results.p_coef = p_coef_vec;

% save results to file
save(strcat(out_dir,'/results.mat'),'results');

%%%%%% FIGURE %%%%%%

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
modules.mod(11).counter=0; modules.mod(11).areas=[];
modules.mod(12).counter=0; modules.mod(12).areas=[];
modules.mod(13).counter=0; modules.mod(13).areas=[];
modules.mod(14).counter=0; modules.mod(14).areas=[];
modules.mod(15).counter=0; modules.mod(15).areas=[];


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

% save figure
set(0, 'CurrentFigure', hubsFig);
print(strcat(out_dir,'/modularity.png'),'-dpng','-r300');

function [equal] = check_structures_equal(Ci1, Ci2, nmod)
c = intersect(Ci1,Ci2,'rows');
diff = nmod - size(c,1);
if diff == 0
    equal = 1;
else
    equal = 0;
end