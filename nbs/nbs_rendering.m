function [component] = nbs_rendering(connmat_pop1, connmat_pop2, tstats_thr, pval_thr,atlas, test, covariates)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB script for running and rending outputs of NBS, for the alternative
% hypothesis that the means of both groups are not equal. 
%
% Please note that thm
% All significant connected components will be plotted and colour-coded
% differently. In case only a single connected component is found, the
% edges will be colour-coded according to their directionality, i.e.
% increase (SC1<SC2) or decrease (SC1>SC2) in connectivity.
% 
% Inputs:
%   - connmat_pop1/2:   3D weighted connectivity matrices for each group.
%      (mandatory)      If FC matrices, negative values will be nulled 
%                       (optional: negative values converted to positive) 
% 
%   - tstats_thr:       threshold for the t-stats used by nbs_bct.m to
%                       identify connected components.
%
%   - pval_thr:         determines the threshold for determining the
%                       significance of a connected component.
%
%   - atlas:            atlas to be used: 1 for shen atlas, 2 for aal, 
%                       3 for myatlas, 4 for shen with no cerebellum
%                       5 for silhouette parcellation
%                      
%   - test:             type of test to be performed: '2sample' for 2sample
%                       t-test, 'paired' for paired sample t-test
%
% Output:
%   - component:        Structure with information on the connections 
%                       comprising each of the identified significant
%                       components. Additionally, adjacent, weighted,
%                       weighted-normalised, and binary matrices are
%                       stored for each of the resulting components.
%
%
% Henrique Fernandes & Ana Coelho 2017
% (last update 09/2020)
%
% email: henrique.fernandes@psych.ox.ac.uk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: include tailed option
%   I - CHECKING INPUTS
if nargin<3
    error('Insufficient number of inputs. Check the required input structure before running the script again.')
elseif nargin<4
    if ~isequal(size(connmat_pop1), size(connmat_pop2)); error('Different size matrices. Please check that the connectivity matrices from pop1 and pop2 are of equal size.'); end
    if isempty(pval_thr); pval_thr = 0.05; end
end

%   II - LOADING: cortex & labels & COGs
addpath(genpath('./utils'));
load('./utils/cortex/Cortex_patch.mat.');

switch atlas %load labels and COGs of respective atlas
    case 1
        load('./utils/parcelations/SH_labels_lat_H.mat');
        load('./utils/parcelations/SH_cog.txt');
        MNI_coord = SH_cog; clear SH_cog
    case 2
        load('./utils/parcelations/AAL_labels_lat_H_90.mat');
        load('./utils/parcelations/aal-90_cog.txt');
        MNI_coord = aal_90_cog; clear aal_90_cog
    case 3
        load('./utils/parcelations/myatlas_labels_lat_H.mat');
        load('./utils/parcelations/myatlas_cog.txt');
        MNI_coord = myatlas_cog; clear myatlas_cog
    case 4
        load('./utils/parcelations/SH_nocerebellum_labels_lat_H.mat');
        load('./utils/parcelations/SH_nocerebellum_cog.txt');
        MNI_coord = SH_nocerebellum_cog; clear SH_nocerebellum_cog 
    case 5
        load('./utils/parcelations/silhouette_thr300_lat_labels.mat');
        load('./utils/parcelations/silhouette_thr300_cog.txt');
        MNI_coord = silhouette_thr300_cog; clear silhouette_thr300_cog
end


%   III - FIGURE PROPERTIES: glass brain rendering
figure('units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto','color','white'); hold on;
patch('vertices',cortexvert,'faces',cortexface,'EdgeColor','none','FaceColor',[0.4 0.4 .4],'FaceAlpha',0.08); %plots a mesh surface for the cortex
material dull; 

%   IV - EDGE PROPERTIES: colour and scaling factor
colormap HSV
cm_cmp     = colormap;
colour_arr = [49 37 7 19 43 55 13 1 62];
s_factor   = 3;
monocrom   = 0;

%   V - CONNECTOME TYPE VERIFICATION: checking if FC matrix with negative values
if sum(connmat_pop1(:)<0)>0 || sum(connmat_pop2(:)<0)>0
    fprintf('Negative values have been found in one or either of the input matrices. Inputs are assumed to be FC matrices\n');
    % null negs
    fprintf('Negative values will be nulled.');
    connmat_pop1(connmat_pop1<0) = 0;
    connmat_pop2(connmat_pop2<0) = 0;
%     % absolute
%     fprintf('Negative values will be transformed to positive.\n');
%     connmat_pop1 = abs(connmat_pop1);
%     connmat_pop2 = abs(connmat_pop2);
end

%   VI - RUNNING NBS - alternative hypothesis is: means are not equal.
if ~isempty(covariates)
    [nbs_pvals,nbs_adj] = nbs_bct(connmat_pop1,connmat_pop2, covariates, tstats_thr, test);
else
    [nbs_pvals,nbs_adj] = nbs_bct(connmat_pop1,connmat_pop2, tstats_thr, test);
end
    % load('test/nbs_adj.mat');load('test/nbs_pvals')

%   VII - CHECKING COMPONENTS
if sum(nbs_pvals<pval_thr) == 0
    fprintf('No significant components found. Try using different t-stats threshold.')
    component={}; % empty component to retrieve when there aren't significant results
    return;
elseif sum(nbs_pvals<pval_thr) > 1
    monocrom = 1;
end
%   matrix of (difference) weights
connmat_diff     = mean(connmat_pop2,3)-mean(connmat_pop1,3);
%   vector of codes for existing significant components
cmp_sign_idx_vec = find(nbs_pvals<pval_thr);


%   VIII - PRINTING COMPONENTS AND TABLES
for cmp=1:length(cmp_sign_idx_vec)
    
    nbs_adj_tmp = full(nbs_adj);
    nbs_adj_tmp(nbs_adj_tmp~=cmp_sign_idx_vec(cmp)) = 0;
    nbs_adj_tmp(nbs_adj_tmp==cmp_sign_idx_vec(cmp)) = 1; %binarize
    mask = tril(ones(size(nbs_adj_tmp,1)),-1);

    component{cmp}.adj      = nbs_adj_tmp;
    component{cmp}.bin      = component{cmp}.adj.*mask;
    component{cmp}.wei      = component{cmp}.bin.*connmat_diff;
    component{cmp}.wei_norm = component{cmp}.wei/max(abs(connmat_diff(:)));
    component{cmp}.pval = nbs_pvals(cmp_sign_idx_vec(cmp));
    
    [component{cmp}.area1.number,component{cmp}.area2.number,component{cmp}.diff_vec] = find(component{cmp}.wei);
    
    
%    component{cmp}.area1.label = aal_lat_labels(component{cmp}.area1.number);
%    component{cmp}.area2.label = aal_lat_labels(component{cmp}.area2.number); 
    
    switch atlas %get labels of respective atlas
        case 1
            component{cmp}.area1.label = SH_lat_labels(component{cmp}.area1.number);
            component{cmp}.area2.label = SH_lat_labels(component{cmp}.area2.number);
            
        case 2
            component{cmp}.area1.label = aal_lat_labels(component{cmp}.area1.number);
            component{cmp}.area2.label = aal_lat_labels(component{cmp}.area2.number);
            
        case 3
            component{cmp}.area1.label = myatlas_lat_labels(component{cmp}.area1.number);
            component{cmp}.area2.label = myatlas_lat_labels(component{cmp}.area2.number);
            
        case 4
            component{cmp}.area1.label = SH_nocerebellum_lat_labels(component{cmp}.area1.number);
            component{cmp}.area2.label = SH_nocerebellum_lat_labels(component{cmp}.area2.number);
    end

    component{cmp}.conn_table  = [ arrayfun(@num2str,component{cmp}.area1.number, 'UniformOutput', false) component{cmp}.area1.label arrayfun(@num2str,component{cmp}.area2.number, 'UniformOutput', false) component{cmp}.area2.label];
 
    component{cmp}.colour=cm_cmp(colour_arr(cmp_sign_idx_vec(cmp)),:); 
    
    for l=1:length(component{cmp}.area1.number)
        if monocrom
            edge_colour = component{cmp}.colour;
        else
            if component{cmp}.wei_norm(component{cmp}.area1.number(l),component{cmp}.area2.number(l))<0
                edge_colour = [0 0.34 1];  % blue; decrease in connectivity; SC1>SC2
            else
                edge_colour = [1 0.18 0];  % red; increase in connectivity; SC1<SC2
            end
        end
        edge_width  = component{cmp}.wei_norm(component{cmp}.area1.number(l),component{cmp}.area2.number(l))*s_factor;
        [X1,Y1,Z1]  = cylinder1(MNI_coord(component{cmp}.area1.number(l),:),MNI_coord(component{cmp}.area2.number(l),:),edge_width,100);
        wei_conn_lr = surf(X1,Y1,Z1,'FaceColor',edge_colour,'EdgeColor','none');
        sphere2(MNI_coord(component{cmp}.area1.number(l),1),MNI_coord(component{cmp}.area1.number(l),2),MNI_coord(component{cmp}.area1.number(l),3),edge_width,edge_colour,1); %text(MNI_coord(i,1)-2,MNI_coord(i,2)-3,MNI_coord(i,3),int2str(i),'FontWeight','bold','FontSize',12);
        sphere2(MNI_coord(component{cmp}.area2.number(l),1),MNI_coord(component{cmp}.area2.number(l),2),MNI_coord(component{cmp}.area2.number(l),3),edge_width,edge_colour,1); %text(MNI_coord(j,1)-2,MNI_coord(j,2)-6,MNI_coord(j,3),int2str(j),'FontWeight','bold','FontSize',12);
    end
    
    fprintf('\n\nNetwork of %s Connected components (pval: %s) comprising the following connections:\n\n', int2str(l), num2str(nbs_pvals(cmp_sign_idx_vec(cmp))));
    component{cmp}.conn_table
end

axis equal
axis off
camlight;


