function [metrics] = graph_theory_metrics(in_file, nnodes)

%
% Function to compute graph theory metrics
%
% Metrics calculated:
%   - global, nodal and local efficiency
%   - normalised nodal efficiency
%   - connection density
%   - degree (for each node and sum of all nodes)
%   - mean connectivity
%
% Inputs:
%           in_file: file with connectivity matrices (3D array: subjects'
%           connectivity matrices concatenated in the 3rd dimension)
%
%           nnodes: number of nodes of the parcellation used
%
% Outputs:
%           metrics: struct with graph theory metrics
%
%
% Ana Coelho 06-02-2020
%

% load file with connectivity matrices
mats_file = load(in_file);
name = fieldnames(mats_file);
mats = mats_file.(name{1,1});

% Graph theory metrics

% global efficiency
Eglob=zeros(size(mats,3),1);
for i=1:size(mats,3)
    Eglob(i,1)=efficiency_wei(mats(:,:,i));
end
metrics.Eglob = Eglob;

% nodal efficiency
Enodal=zeros(size(mats,3),nnodes);
for i=1:size(mats,3)
    Enodal(i,:)=efficiency_wei(mats(:,:,i),1);
end
metrics.Enodal = Enodal;

% local efficiency
Elocal = mean(Enodal,2);
metrics.Elocal = Elocal;

% normalised nodal efficiency
Enodal_norm=zeros(size(mats,3),nnodes);
for i=1:size(mats,3)
    Enodal_norm(i,:)=efficiency_wei(mats(:,:,i),1)/metrics.Elocal(i,1);
end
metrics.Enodal_norm = Enodal_norm;

% connection density
density=zeros(size(mats,3),1);
for i=1:size(mats,3)
    density(i,1)=density_und(mats(:,:,i));
end
metrics.density = density;

% degree (each node)
degree=zeros(size(mats,3),nnodes);
for i=1:size(mats,3)
    degree(i,:)=degrees_und(mats(:,:,i));
end
metrics.degree = degree;

% degree (sum of all nodes)
sum_deg = sum(degree,2);
metrics.sum_degree = sum_deg;

% mean connectivity
mean_sc=zeros(size(mats,3),1);
for i=1:size(mats,3)
    m=triu(mats(:,:,i),1);
    mean_sc(i,1)=mean(mean(m));
end
metrics.mean_connectivity = mean_sc;

end