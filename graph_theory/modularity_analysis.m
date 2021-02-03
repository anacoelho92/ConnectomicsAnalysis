
function [modules] = modularity_analysis(SC, Ci, cHubs, pop_label,outdir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB script for running intra- and inter-modularity connectivity 
% analysis given SC matrices, and a community partition structure. Optional
% connector Hub to modules connectivity analysis.
% 
% INPUTS:
%   - SC:               3D weighted structural connectivity matrices for 
%                       each group (subject's individual SC matrices), or 
%                       2D group SC matrix.       
% 
%   - Ci:               modularity partiotioning. Vector of module-codes
%                       assigned to each area of the parcellation.
%                         
%   - cHubs:            vector of connector hubs.
%     (optional)
%
%   - pop_label:        Group ID label. Used to output variable and file
%                       name definition.
%
%   - outdir:           output directory to which the resulting analysis
%     (optional)        output file (structure-type variable) will be
%                       saved. If empty, current directory will be assumed
%                       as output directory.
%
% OUTPUTS:
%   - modules:        structure containing the fields 'Ci', 'cHubs', and a
%                     summary of the modularity connectiivty analyses:
%                     - 'modules': intra/inter-modular (subject- and 
%                                  group-level)
%                     - 'cHubs2modules': connectorHub to modules
%
%
% USAGE: 
%       ex: modularity_analysis(subjects_PBD.SC_3D,modules_HC.Ci,...
%               modules_PBD.cHubs,'modules')
%
%
% Henrique Fernandes & Ana Coelho 2018
% (last update 01/2021)
%
% email: henrique.fernandes@clin.au.dk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(outdir)
    oudir = pwd;
end

%% Computing subject/group Modular Network Analysis

for mods=1:size(unique(Ci))
    
    % Defining modules as binary vectors (i.e. areas belonging to each module) 
    modules{mods}.nodes   = find(Ci==mods);
    modules{mods}.nodevec = ismember(Ci,mods);

    % If 3D matrix (i.e. individual subject's matrix)
    if size(SC,3)>1
        for subj=1:size(SC,3)
            % WITHIN-MODULE GLOBAL (MEAN) METRICS
            % Network
            modules{mods}.subjects.SCmat{subj}               = SC(modules{mods}.nodevec,modules{mods}.nodevec',subj);
            % Strength
            modules{mods}.subjects.strength_total(subj)      = sum(sum(SC(modules{mods}.nodevec,modules{mods}.nodevec',subj),2));
            modules{mods}.subjects.strength_total_mean(subj) = mean(sum(SC(modules{mods}.nodevec,modules{mods}.nodevec',subj),2));
            % Mean Connectivity
            modules{mods}.subjects.strength_mean(subj)       = mean(mean(SC(modules{mods}.nodevec,modules{mods}.nodevec',subj),2));
            % Degree
            modules{mods}.subjects.degree(subj)              = sum(degrees_und(SC(modules{mods}.nodevec,modules{mods}.nodevec',subj)')');
            modules{mods}.subjects.degree_mean(subj)         = mean(degrees_und(SC(modules{mods}.nodevec,modules{mods}.nodevec',subj)')');
            % Within-module Degree z-score (within module version of degree centrality)
            moddeg_tmp = module_degree_zscore(SC(:,:,subj),Ci);
            modules{mods}.subjects.module_deg_zcore(subj)    = mean(moddeg_tmp(modules{mods}.nodevec)); 
            % Participation coefficient - measure of diversity of intermodular connections of individual nodes 
            pcoef_tmp  = participation_coef(SC(:,:,subj),Ci);
            modules{mods}.subjects.pcoef(subj)               = mean(pcoef_tmp(modules{mods}.nodevec));
            % Betweenness centrality  - fraction of all shortest paths in the network that contain a given node
            bwc_tmp    = betweenness_wei(SC(:,:,subj));
            modules{mods}.subjects.bwc(subj)                 = mean(bwc_tmp(modules{mods}.nodevec)); 
            clear moddeg_tmp pcoef_tmp moddeg_tmp;
        end
        SCavg = mean(SC,3);
    else
        SCavg = SC;
    end
        SCavg_norm = SCavg/max(abs(SCavg(:)));
 
    % WITHIN-MODULE GLOBAL (MEAN) METRICS
    % Network
    modules{mods}.group.SCmat               = SCavg_norm(modules{mods}.nodevec,modules{mods}.nodevec');
    % Node count
    modules{mods}.group.N_nodes             = sum(modules{mods}.nodevec);
    % Strength
    modules{mods}.group.strength_total      = sum(sum(SCavg_norm(modules{mods}.nodevec,modules{mods}.nodevec'),2));
    modules{mods}.group.strength_total_mean = mean(sum(SCavg_norm(modules{mods}.nodevec,modules{mods}.nodevec'),2));
    % Mean Connectivity
    modules{mods}.group.strength_mean       = mean(mean(SCavg_norm(modules{mods}.nodevec,modules{mods}.nodevec'),2));
    % Degree
    modules{mods}.group.degree              = sum(degrees_und(SCavg_norm(modules{mods}.nodevec,modules{mods}.nodevec')')');
    modules{mods}.group.degree_mean         = mean(degrees_und(SCavg_norm(modules{mods}.nodevec,modules{mods}.nodevec')')');
    % Within-module Degree z-score (within module version of degree centrality)
    moddeg_tmp = module_degree_zscore(SCavg_norm(:,:),Ci);
    modules{mods}.group.module_deg_zcore    = mean(moddeg_tmp(modules{mods}.nodevec)); 
    % Participation coefficient - measure of diversity of intermodular connections of individual nodes 
    pcoef_tmp  = participation_coef(SCavg_norm(:,:),Ci);
    modules{mods}.group.pcoef               = mean(pcoef_tmp(modules{mods}.nodevec));
    % Betweenness centrality  - fraction of all shortest paths in the network that contain a given node
    bwc_tmp    = betweenness_wei(SCavg_norm(:,:));
    modules{mods}.group.bwc                 = mean(bwc_tmp(modules{mods}.nodevec)); 
end


% INTER- AND INTRA-MODULAR network properties (strength and degree)
for m=1:size(unique(Ci))
    for n=1:size(unique(Ci))
        if m~=n   % INTER-MODULE
            modularity.modConn.intermod_strength(m,n) = sum(sum(SCavg_norm(modules{m}.nodes,modules{n}.nodes)));
            modularity.modConn.intermod_degree(m,n)   = sum(degrees_und(SCavg_norm(modules{m}.nodes,modules{n}.nodes)));
            
            % Connectity between the connector Hubs of a module and all other modules 
            if ~isempty(cHubs)
                mod_cHubs = intersect(cHubs,modules{m}.nodes);
                modularity.modConn.intermod_cHub_strength(m,n) = sum(sum(SCavg_norm(mod_cHubs,modules{n}.nodes)));
                modularity.modConn.intermod_cHub_degree (m,n)  = sum(degrees_und(SCavg_norm(mod_cHubs,modules{n}.nodes)));
            end
        else       % INTRA-MODULE
            modularity.modConn.intramod_strength(m,n) = sum(sum(SCavg_norm(modules{m}.nodes,modules{n}.nodes)));
            modularity.modConn.intramod_degree(m,n)   = sum(degrees_und(SCavg_norm(modules{m}.nodes,modules{n}.nodes)));
        end
        modularity.modConn.strength(m,n) = sum(sum(SCavg_norm(modules{m}.nodes,modules{n}.nodes)));
        modularity.modConn.degree(m,n)   = sum(degrees_und(SCavg_norm(modules{m}.nodes,modules{n}.nodes)));
        
        clear mod_cHubs;
    end
    
end

%% Plotting connectivity: INTRA/INTER-MODULAR + CHUBS-MODULES 

% INTRAMODULAR CONNECTIVITY
% Degree
%figure; imagesc((modularity.modConn.intramod_degree))
%figure; imagesc(sum(modularity.modConn.intramod_degree))
%intramod_deg_fig = figure; heatmap(modularity.modConn.intramod_degree,'Colormap',jet,'XLabel','modules','YLabel','modules')
intramod_deg_fig = create_heatmap(modularity.modConn.intramod_degree,size(unique(Ci)),size(unique(Ci)),'modules','modules','%d');

% Strength
%figure; imagesc((modularity.modConn.intramod_strength))
%figure; imagesc(sum(modularity.modConn.intramod_strength))
%intramod_str_fig = figure; heatmap(modularity.modConn.intramod_strength,'Colormap',jet,'XLabel','modules','YLabel','modules')
intramod_str_fig = create_heatmap(modularity.modConn.intramod_strength,size(unique(Ci)),size(unique(Ci)),'modules','modules','%0.2f');

% INTERMODULAR CONNECTIVITY
% Degree
%figure; imagesc((modularity.modConn.intermod_degree))
%figure; imagesc(sum(modularity.modConn.intermod_degree))
%intermod_deg_fig = figure; heatmap(modularity.modConn.intermod_degree,'Colormap',jet,'XLabel','modules','YLabel','modules')
intermod_deg_fig = create_heatmap(modularity.modConn.intermod_degree,size(unique(Ci)),size(unique(Ci)),'modules','modules','%d');

%intermod_deg_mean_fig = figure('Renderer', 'painters', 'Position', [10 10 100 450]); heatmap(sum(modularity.modConn.intermod_degree,2),'Colormap',jet,'XLabel','mean','XDisplayLabels',{''},'YDisplayLabels',repmat({''},1,length(unique(Ci))))
intermod_deg_mean_fig = create_heatmap(sum(modularity.modConn.intermod_degree,2),1,size(unique(Ci)),'mean','','%d');
set(0, 'currentfigure', intermod_deg_mean_fig);
set(gcf,'position',[10,10,100,400]);
set(gca,'XTick',[],'xticklabel',[],'YTick',[],'yticklabel',[],'xaxisLocation','top');
%ax = gca;
%axp = struct(ax);       
%axp.Axes.XAxisLocation = 'top';

% Strength
%figure; imagesc((modularity.modConn.intermod_strength))
%figure; imagesc(sum(modularity.modConn.intermod_strength))
%intermod_str_fig = figure; heatmap(modularity.modConn.intermod_strength,'Colormap',jet,'XLabel','modules','YLabel','modules')
intermod_str_fig = create_heatmap(modularity.modConn.intermod_strength,size(unique(Ci)),size(unique(Ci)),'modules','modules','%0.2f');

%intermod_str_mean_fig = figure('Renderer', 'painters', 'Position', [10 10 100 450]); heatmap(sum(modularity.modConn.intermod_strength,2),'Colormap',jet,'XLabel','mean','XDisplayLabels',{''},'YDisplayLabels',repmat({''},1,length(unique(Ci))))
intermod_str_mean_fig = create_heatmap(sum(modularity.modConn.intermod_strength,2),1,size(unique(Ci)),'mean','','%0.2f');

set(0, 'currentfigure', intermod_str_mean_fig);
set(gcf,'position',[10,10,100,400]);
set(gca,'XTick',[],'xticklabel',[],'YTick',[],'yticklabel',[],'xaxisLocation','top');

%ax = gca;
%axp = struct(ax);       
%axp.Axes.XAxisLocation = 'top';

% CONNECTOR-HUB DRIVEN INTERMODULARITY
% Degree
%figure; imagesc(modularity.modConn.intermod_cHub_degree)
%figure; imagesc(sum(modularity.modConn.intermod_cHub_degree,2))
%connintermod_deg_fig = figure; heatmap(modularity.modConn.intermod_cHub_degree,'Colormap',jet,'XLabel','modules','YLabel','modules')
connintermod_deg_fig = create_heatmap(modularity.modConn.intermod_cHub_degree,size(unique(Ci)),size(unique(Ci)),'modules','modules','%d');

%connintermod_deg_mean_fig = figure('Renderer', 'painters', 'Position', [10 10 100 450]); heatmap(sum(modularity.modConn.intermod_cHub_degree,2),'Colormap',jet,'XLabel','mean','XDisplayLabels',{''},'YDisplayLabels',repmat({''},1,length(unique(Ci))))
connintermod_deg_mean_fig = create_heatmap(sum(modularity.modConn.intermod_cHub_degree,2),1,size(unique(Ci)),'mean','','%d');

set(0, 'currentfigure', connintermod_deg_mean_fig);
set(gcf,'position',[10,10,100,400]);
set(gca,'XTick',[],'xticklabel',[],'YTick',[],'yticklabel',[],'xaxisLocation','top');

% ax = gca;
% axp = struct(ax);       
% axp.Axes.XAxisLocation = 'top';

% Strength
%figure; imagesc(modularity.modConn.intermod_cHub_strength)
%figure; imagesc(sum(modularity.modConn.intermod_cHub_strength,2))
%connintermod_str_fig = figure; heatmap(modularity.modConn.intermod_cHub_strength,'Colormap',jet,'XLabel','modules','YLabel','modules')
%connintermod_str_fig = create_heatmap(modularity.modConn.intermod_cHub_strength,size(unique(Ci)),size(unique(Ci)),'modules','modules','%0.2f');

%connintermod_str_mean_fig = figure('Renderer', 'painters', 'Position', [10 10 100 450]); heatmap(sum(modularity.modConn.intermod_cHub_strength,2),'Colormap',jet,'XLabel','mean','XDisplayLabels',{''},'YDisplayLabels',repmat({''},1,length(unique(Ci))))
%connintermod_str_mean_fig = create_heatmap(sum(modularity.modConn.intermod_cHub_strength,2),1,size(unique(Ci)),'mean','','%0.2f');

%set(0, 'currentfigure', connintermod_str_mean_fig);
%set(gcf,'position',[10,10,100,400]);
%set(gca,'XTick',[],'xticklabel',[],'YTick',[],'yticklabel',[],'xaxisLocation','top');

%ax = gca;
%axp = struct(ax);       
%axp.Axes.XAxisLocation = 'top';


set(0, 'currentfigure', intramod_deg_fig); print(intramod_deg_fig,strcat(outdir,'intramod_deg.png'),'-dpng','-r300');
%set(0, 'currentfigure', intramod_str_fig); print(intramod_str_fig,strcat(outdir,'intramod_str.png'),'-dpng','-r300');
set(0, 'currentfigure', intermod_deg_fig); print(intermod_deg_fig,strcat(outdir,'intermod_deg.png'),'-dpng','-r300');
set(0, 'currentfigure', intermod_deg_mean_fig);
oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
     'PaperPosition',[newpos(1),newpos(2),1.5*newpos(3),newpos(4)])
print(intermod_deg_mean_fig,strcat(outdir,'intermod_deg_mean.png'),'-dpng','-r300');
drawnow
% set(gcf,'Units',oldscreenunits,...
%      'PaperUnits',oldpaperunits,...
%      'PaperPosition',oldpaperpos)
% %set(0, 'currentfigure', intermod_str_fig); print(intermod_str_fig,strcat(outdir,'intermod_str.png'),'-dpng','-r300');
% set(0, 'currentfigure', intermod_str_mean_fig);
% oldscreenunits = get(gcf,'Units');
% oldpaperunits = get(gcf,'PaperUnits');
% oldpaperpos = get(gcf,'PaperPosition');
% set(gcf,'Units','pixels');
% scrpos = get(gcf,'Position');
% newpos = scrpos/100;
% set(gcf,'PaperUnits','inches',...
%      'PaperPosition',[newpos(1),newpos(2),1.5*newpos(3),newpos(4)])
% print(intermod_str_mean_fig,strcat(outdir,'intermod_str_mean.png'),'-dpng','-r300');
% drawnow
% set(gcf,'Units',oldscreenunits,...
%      'PaperUnits',oldpaperunits,...
%      'PaperPosition',oldpaperpos)
set(0, 'currentfigure', connintermod_deg_fig); print(connintermod_deg_fig,strcat(outdir,'connintermod_deg.png'),'-dpng','-r300');
set(0, 'currentfigure', connintermod_deg_mean_fig);
oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
     'PaperPosition',[newpos(1),newpos(2),1.5*newpos(3),newpos(4)])
print(connintermod_deg_mean_fig,strcat(outdir,'connintermod_deg_mean.png'),'-dpng','-r300');
drawnow
% set(gcf,'Units',oldscreenunits,...
%      'PaperUnits',oldpaperunits,...
%      'PaperPosition',oldpaperpos)
%  
% %set(0, 'currentfigure', connintermod_str_fig); print(connintermod_str_fig,strcat(outdir,'connintermod_str.png'),'-dpng','-r300');
% set(0, 'currentfigure', connintermod_str_mean_fig); 
% %set(0, 'currentfigure', connintermod_deg_mean_fig);
% oldscreenunits = get(gcf,'Units');
% oldpaperunits = get(gcf,'PaperUnits');
% oldpaperpos = get(gcf,'PaperPosition');
% set(gcf,'Units','pixels');
% scrpos = get(gcf,'Position');
% newpos = scrpos/100;
% set(gcf,'PaperUnits','inches',...
%      'PaperPosition',[newpos(1),newpos(2),1.5*newpos(3),newpos(4)])
% print(connintermod_str_mean_fig,strcat(outdir,'connintermod_str_mean.png'),'-dpng','-r300');
% drawnow
% set(gcf,'Units',oldscreenunits,...
%      'PaperUnits',oldpaperunits,...
%      'PaperPosition',oldpaperpos)
%% Creating a variable name based on 'pop_label'

eval(sprintf('modules_%s.Ci = Ci', pop_label))
eval(sprintf('modules_%s.cHubs = cHubs', pop_label))
eval(sprintf('modules_%s.modules = modules', pop_label))
eval(sprintf('modules_%s.cHubs2modules = modularity.modConn', pop_label))

%% Saving structure to outdir

save([outdir '/module_fingerprint_' pop_label],['modules_' pop_label])


end