function [Eff_norm_asc,hbar_labels, Eff_error_asc, global_hubs] = get_global_hubs(Eff, labels)

% Function to get global hubs and sorted nodal efficiency normalized values
%
%
% Inputs:   Eff: matrix with local efficiency values (n_nodes x n_subjects)
%
%           labels: labels of the parcellation used
%
%
% Outputs:  Eff_norm_asc: normalized local efficiency values in ascending
%           order
%
%           hbar_labels : node labels sorted in ascending order of
%           normalized local efficiency values
%
%           Eff_error_asc: standard deviation of normalized local 
%           efficiency values in ascending order
%
%           global_hubs: IDs of nodes identified as global hubs
%
%
% Ana Coelho 11-01-2021
%
%

n_nodes = size(Eff,1);
Eff_lo.mean=zeros(n_nodes,1);
Eff_lo.stdev=zeros(n_nodes,1);

for n=1:n_nodes
    Eff_lo.mean(n)=mean(squeeze(Eff(n,:)));
    Eff_lo.stdev(n)=std2(squeeze(Eff(n,:)));
end

[Eff_nodal,Eff_norm]=nodal_efficiency_lo_HF(Eff_lo);
global_hubs = find(Eff_nodal);

[Eff_norm_asc, IDx] = sort(Eff_norm.mean,'ascend');
hbar_labels         = labels(IDx);
Eff_error_asc       = Eff_norm.stdev(IDx);

end