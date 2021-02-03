function [results, mean_std_stats_pop1, mean_std_stats_pop2] = stats_graph_metrics(metrics_pop1, metrics_pop2, type)
%
% Function to perform statistical analysis of graph metrics
%
% Metrics analysed:
%   - global efficiency
%   - local efficiency
%   - connection density
%   - degree
%   - mean connectivity
%
% Inputs:
%           metrics_pop1, metrics_pop2: structure with graph metrics for
%           each population (output of graph_theory_metrics.m)
%
%           type: type of test ('ttest2' for two sample t-test,
%           'paired-ttest' for paired t-test)
%
% Outputs:
%           results: struct with test statistic and p-values for the
%           different tests
%
%           mean_std_stats_pop1/pop2: struct with mean and std for each
%           metric
%
% Ana Coelho 10-03-2020
%


switch (type)
    case 'ttest2'
        % Global Efficiency
        [h,p,ci,stats]=ttest2(metrics_pop1.Eglob,metrics_pop2.Eglob);
        results.Eglob=[stats.tstat,p];

        % Local Efficiency
        [h,p,ci,stats]=ttest2(metrics_pop1.Elocal,metrics_pop2.Elocal);
        results.Elocal=[stats.tstat,p];

        % Connection density
        [h,p,ci,stats]=ttest2(metrics_pop1.density,metrics_pop2.density);
        results.density=[stats.tstat,p];

        % Degree
        [h,p,ci,stats]=ttest2(mean(metrics_pop1.degree,2),mean(metrics_pop2.degree,2));
        results.degree=[stats.tstat,p];

        % Mean Connectivity
        [h,p,ci,stats]=ttest2(metrics_pop1.mean_connectivity,metrics_pop2.mean_connectivity);
        results.mean_connectivity=[stats.tstat,p];
        
    case 'paired-ttest'
        % Global Efficiency
        [h,p,ci,stats]=ttest(metrics_pop1.Eglob,metrics_pop2.Eglob);
        results.Eglob=[stats.tstat,p];

        % Local Efficiency
        [h,p,ci,stats]=ttest(metrics_pop1.Elocal,metrics_pop2.Elocal);
        results.Elocal=[stats.tstat,p];

        % Connection density
        [h,p,ci,stats]=ttest(metrics_pop1.density,metrics_pop2.density);
        results.density=[stats.tstat,p];

        % Degree
        [h,p,ci,stats]=ttest(mean(metrics_pop1.degree,2),mean(metrics_pop2.degree,2));
        results.degree=[stats.tstat,p];

        % Mean Connectivity
        [h,p,ci,stats]=ttest(metrics_pop1.mean_connectivity,metrics_pop2.mean_connectivity);
        results.mean_connectivity=[stats.tstat,p];
end

% compute mean and std of each metric
mean_std_stats_pop1=struct();
mean_std_stats_pop2=struct();

% Mean Connectivity
m_conn = mean(metrics_pop2.mean_connectivity);
std_conn = std(metrics_pop2.mean_connectivity);
mean_std_stats_pop2.mean_connectivity=[m_conn,std_conn];

m_conn = mean(metrics_pop1.mean_connectivity);
std_conn = std(metrics_pop1.mean_connectivity);
mean_std_stats_pop1.mean_connectivity=[m_conn,std_conn];

% Degree
m_deg = mean(metrics_pop2.sum_degree);
std_deg = std(metrics_pop2.sum_degree);
mean_std_stats_pop2.degree=[m_deg,std_deg];

m_deg = mean(metrics_pop1.sum_degree);
std_deg = std(metrics_pop1.sum_degree);
mean_std_stats_pop1.degree=[m_deg,std_deg];

% Connection Density
m_dens = mean(metrics_pop2.density);
std_dens = std(metrics_pop2.density);
mean_std_stats_pop2.density=[m_dens,std_dens];

m_dens = mean(metrics_pop1.density);
std_dens = std(metrics_pop1.density);
mean_std_stats_pop1.density=[m_dens,std_dens];

% Global Efficiency
m_eglob = mean(metrics_pop2.Eglob);
std_eglob = std(metrics_pop2.Eglob);
mean_std_stats_pop2.Eglob=[m_eglob,std_eglob];

m_eglob = mean(metrics_pop1.Eglob);
std_eglob = std(metrics_pop1.Eglob);
mean_std_stats_pop1.Eglob=[m_eglob,std_eglob];

% Local Efficiency
m_eloc = mean(metrics_pop2.Elocal);
std_eloc = std(metrics_pop2.Elocal);
mean_std_stats_pop2.Elocal=[m_eloc,std_eloc];

m_eloc = mean(metrics_pop1.Elocal);
std_eloc = std(metrics_pop1.Elocal);
mean_std_stats_pop1.Elocal=[m_eloc,std_eloc];