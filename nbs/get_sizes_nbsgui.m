function [n_links, n_nodes] = get_sizes_nbsgui(results, tstats, out_file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB script to extract number of nodes, links and components of nbs 
% results obtained using the GUI 
%
% input:
%         results: cell with results from different threshold values
%         tstats: array of tstat values used to get results
%         out_file: basename for output files to save plots
%
% outputs:
%         n_links: number of links of the largest component for each tstat
%         n_nodes: number of nodes of the largest component for each tstat
%
% Ana Coelho 10-01-2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n_links = zeros(length(results),1);
n_nodes = zeros(length(results),1);
n_cmps = zeros(length(results),1);

for k=1:length(results)
    
    % get the largest component
    max_size = 0;
    max_cmp = 0;
    
    if results{1,k}.NBS.n>0
        for i=1:results{1,k}.NBS.n
            adj=results{1,k}.NBS.con_mat{1,i};
            adj_final=full(adj);
            sz_adj = length(find(adj_final ~= 0)); % number of links in the component
            
            % compare number of links with maximum size
            if sz_adj>max_size
                max_cmp = i;
                max_size=sz_adj;
            end
            
        end
        
        % get adjacent matrix of the largest component
        adj = results{1,k}.NBS.con_mat{1,max_cmp};
        adj_final = full(adj); % lower triangle of matrix
        
        [a,b]=find(adj_final~=0); % get connections different from zero
        nodes = cat(1,a,b); % concatenate nodes of all connections
        n_nodes(k,1) = length(unique(nodes)); % number of nodes
        n_links(k,1) = length(a); % number of links
        n_cmps(k,1) = results{1,k}.NBS.n; % number of significant components
    end
end

% plot for number of nodes
plot(tstats,n_nodes,'-ob','MarkerFaceColor','blue')
grid on
xlabel('F-threshold')
ylabel('Number of nodes')
ax=gca;
ax.FontSize=13;
print(gcf,strcat(out_file,'_nnodes.png'),'-dpng','-r300');

% plot for number of links
plot(tstats,n_links,'-ob','MarkerFaceColor','blue')
grid on
xlabel('F-threshold')
ylabel('Number of connections')
ax=gca;
ax.FontSize=13;
print(gcf,strcat(out_file,'_nlinks.png'),'-dpng','-r300');

% plot for number of components
plot(tstats,n_cmps,'-ob','MarkerFaceColor','blue')
grid on
xlabel('F-threshold')
ylabel('Number of components')
ax=gca;
ax.FontSize=13;
print(gcf,strcat(out_file,'_ncmps.png'),'-dpng','-r300');

end