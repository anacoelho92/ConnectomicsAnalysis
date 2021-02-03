function [mats_thr_cons] = consistency_analysis(in_file, out_dir, p, range)

%
% Function to create plots to analyse results of threshold consistency
%
% Inputs:
%           in_file: file with connectivity matrices (3D array: subjects' 
%           connectivity matrices concatenated in the 3rd dimension)
%           out_dir: directory to save output files
%           p: proportion of weights to preserve in threshold consistency
%           range: step of connection strength values for plotting 
%           frequencies of connections lost 
%
% Outputs:
%           connections_lost: plot with strength of connections lost when
%           applying the consistency threshold
%
%           connections_mask: plot with strength of connections from group
%           consistency mask that are not present in all subjects
%
%           counts_lost: plot with number of connections lost in each
%           individual when applying the consistency threshold
%
%           counts_mask: plot with number of connections from group
%           consistency mask that are not present in each individual 
%
%           mats_thr_cons: 3D array with each subject connectivity matrix
%           masked with group consistency mask
%
% Ana Coelho 05-02-2020
%

% load file with connectivity matrices
mats_file = load(in_file);
name = fieldnames(mats_file);
mats = mats_file.(name{1,1});

% compute group mean connectivity matrix
mean_group_mats = mean(mats,3);
% apply threshold consistency
mats_thr_cons = threshold_consistency(mats,p);

% matrix with connections lost when applying threshold consistency
diff = mean_group_mats - mats_thr_cons;

% range of connection strengths
x = 0:range:1; % is this range applicable to any case???

% count frequency of connections in each range of strengths
counts_diff=[];
counts_diff(1,1)=length(find(diff==0));
for i=1:length(x)-1
    counts_diff(i+1,1)=length(find(diff>x(i)& diff<=x(i+1)));
end

% get the index from which all counts are 0 
last_id=0;
for i=1:size(counts_diff,1)
    sum_counts = sum(counts_diff(i:end));
    if sum_counts==0
        last_id=i;
        break;
    end
end

% plot histogram with connections lost and save to file
plot(x(2:last_id),counts_diff(2:last_id)); 
hold on;
xlabel('Connection Strength')
ylabel('Counts')
title('Connections lost when applying consistency threshold')
ax=gca;
ax.FontSize=13;
grid on;
print(gcf,strcat(out_dir,'/connections_lost.png'),'-dpng','-r300');
hold off;

%create mask of group consistency
consistency_mask = mats_thr_cons;
consistency_mask(find(consistency_mask>0))=1;

%create inverse of group consistency mask
ones_mat = ones(size(consistency_mask));
inv_cons_mask = ones_mat - consistency_mask;

% apply inverse consistency mask to each subject's matrix to find
% connections present in each subject but not in the group consistency mask
mats_inv_mask=zeros(size(mats));
for i=1:size(mats,3)
    mat=mats(:,:,i);
    mat_mask = mat.*inv_cons_mask;
    mats_inv_mask(:,:,i)=mat_mask;
end

% count number of connections in each subject's matrix not present in group
% consistency mask
counts_inv_mask=zeros(size(mats,3),1);
for i=1:size(mats_inv_mask,3)
    mat=mats_inv_mask(:,:,i);
    counts_inv_mask(i,1)=length(find(mat~=0));
end

% convert counts to percentage
total_inv_mask = length(find(triu(inv_cons_mask,1)==1))*2; % total number of connections in inverse mask
perc_inv_mask = counts_inv_mask / total_inv_mask * 100;

% scatter plot with number of connections lost in each  individual when applying the consistency threshold
scatter(1:length(perc_inv_mask),perc_inv_mask,'filled');
hold on;
xlim([0 (length(perc_inv_mask)+1)])
xlabel('Subject')
ylabel('Percentage of Connections')
title('Connections lost when applying consistency threshold')
ax=gca;
ax.FontSize=13;
grid on;
print(gcf,strcat(out_dir,'/counts_lost.png'),'-dpng','-r300');
hold off;

% indices of connections present in group consistency
ids_cons_mask=find(consistency_mask==1);

% apply group consistency mask to each subject
mats_cons=zeros(size(mats));
for i=1:size(mats,3)
    mats_cons(:,:,i)=mats(:,:,i).*consistency_mask;
end

% get values of connections of consistency mask in each subject
values_cons_mask = zeros(size(mats,3),size(ids_cons_mask,1));
for i=1:size(mats,3)
    mat=mats_cons(:,:,i);
    values = mat(ids_cons_mask);
    values_cons_mask(i,:)=values;
end

% get elements from consistency mask where at least one subject has a
% 0-connection
col=find(sum(values_cons_mask==0));

% get corresponding value in consistency mask of the connections that are not present in all
% subjects
values_not_mask=zeros(size(col,2),1);
for i=1:size(col,2)
    values_not_mask(i,1)=mats_thr_cons(ids_cons_mask(col(i)));
end

% range of connection strengths
% count frequency of connections in each range of strengths
if length(values_not_mask)<10
    x=values_not_mask;
    counts_not_mask=ones(length(x),1);
else
    x=min(values_not_mask):min(values_not_mask)*10:max(values_not_mask)+min(values_not_mask)*10;
    %x=min(values_not_mask):min(values_not_mask)*0.5:max(values_not_mask)+min(values_not_mask)*0.5;
    counts_not_mask=zeros(length(x)-1,1);
    for i=1:length(x)-1
        counts_not_mask(i,1)=length(find(values_not_mask>x(i) & values_not_mask<=x(i+1)));
        % check if this is correct!!!
    end
    x=x(1:end-1);
end

% plot histogram with connections not present in all subjects and save to file
%plot(x(1:end-1),counts_not_mask); 
plot(x,counts_not_mask); 
hold on;
xlabel('Connection Strength')
ylabel('Counts')
title('Connections from group consistency mask not present in all subjects')
ax=gca;
ax.FontSize=13;
grid on;
print(gcf,strcat(out_dir,'/connections_mask.png'),'-dpng','-r300');
hold off;


% get from each subject number of connections from group consistency mask that are 0
zeros_cons_mask=zeros(size(mats,3),1);
for i=1:size(mats,3)
    mat = mats_cons(:,:,i);
    values = mat(ids_cons_mask);
    zeros_cons_mask(i,1)=length(find(values==0));
end

% convert counts to percentage
total_cons_mask = length(ids_cons_mask); % total number of connections in consistency mask
perc_zeros_cons_mask = zeros_cons_mask./total_cons_mask * 100;

% scatter plot with number of connections from group consistency mask that are not present in each individual 
scatter(1:length(perc_zeros_cons_mask),perc_zeros_cons_mask,'filled');
hold on;
xlim([0 (length(perc_inv_mask)+1)])
xlabel('Subject')
ylabel('Percentage of Connections')
title('Connections from group consistency mask not present in all subjects')
ax=gca;
ax.FontSize=13;
grid on;
print(gcf,strcat(out_dir,'/counts_mask.png'),'-dpng','-r300');
hold off;

% save subjects' connectivity matrices masked with group consistency mask
save(strcat(out_dir,'/mats_thr_consistency.mat'),'mats_cons');

end