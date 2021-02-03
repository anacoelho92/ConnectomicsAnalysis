function [vec_recoded, mean_sim] = relabel_mod_partition(Ci1, Ci2)

%   This scripts allows the relabeling of a code/label vector given two
%   partitions
% 
%
%   INPUTS 
%           Ci1/2: affiliation vector of modularity partition
%
%   OUTPUTS
%           vec_recoded: vector labeled (Ci2 will be relabeled according to
%           Ci1)
%
%
%
%   Henrique Fernandes & Ana Coelho 2018 
%   (last update: 05/2020)
%
%   email: henrique.fernandes@clin.au.dk
%


% convert affiliation vectors to binary (rows=modules x columns=brain
% regions)
Ci1_bin = zeros(max(Ci1),size(Ci1,1));
for j=1:size(Ci1,1)
    Ci1_bin(Ci1(j),j)=1;
end

Ci2_bin = zeros(max(Ci2),size(Ci2,1));
for j=1:size(Ci2,1)
    Ci2_bin(Ci2(j),j)=1;
end


% compute similarity matrix between each pair of modules
sim_matrix=zeros(size(Ci1_bin,1),size(Ci2_bin,1));

for i=1:size(Ci1_bin,1)
    for j=1:size(Ci2_bin,1)
        
        % find regions that don't appear in either module
        sum_vector = Ci1_bin(i,:) + Ci2_bin(j,:);
        ids_zero = find(sum_vector==0);
        ids_nonzero = find(sum_vector>0);
        
        % remove non-existent regions from affiliation vector of each
        % module (only maintain regions that exist at least in one of the
        % modules)
        Ci1_final=Ci1_bin(i,:);
        Ci1_final(ids_zero)=[];
        
        Ci2_final=Ci2_bin(j,:);
        Ci2_final(ids_zero)=[];
        
        % compute difference between the two vectors
        diff = Ci1_final - Ci2_final;
        % regions that appear in the two modules will be 0 in the
        % difference vector
        regions = length(find(diff==0));
        % similarity is the number of shared regions divided by the total
        % number of regions in the two modules
        sim_matrix(i,j) = regions / length(diff);
        
    end
end

% get similarity vector (maximum of each column) and correspondence vector
% of source configuration (row id of maximum of each column)
[sim_vec,targetcode] = max(sim_matrix,[],1);
mean_sim = mean(sim_vec);

if max(Ci2) ~= max(Ci1)
    for i=1:length(targetcode)
        if numel(find(targetcode==i))>1
           ids=find(targetcode==i);
           m = min(sim_vec(ids));
           id_min=find(sim_vec==m);
           targetcode(id_min)=length(targetcode);
        end
    end
    
end
   
    
% vector of target configuration is the ordered columns 
origcode = 1:size(Ci2_bin,1);

    
dummyval    = max(targetcode)+1;

% convert binary matrix to affiliation vector
for i=1:size(Ci2_bin,1)
        Ci2_bin(i,:) = Ci2_bin(i,:)*i;
end
vec_recoded = transpose(sum(Ci2_bin,1));
%vec_recoded = Ci2;

% Checks for potential vector intersections
repvec=find(origcode==targetcode);

if length(targetcode) ~= length(unique(targetcode))
    missing=setdiff(origcode,targetcode);
    [~, ind] = unique(targetcode);
    duplicate_ind = setdiff(1:size(targetcode, 2), ind);
    targetcode(duplicate_ind)=missing(1);
    
end
% if size(Ci2_bin,1) < size(Ci1_bin,1)
%     for i=size(Ci2_bin,1)+1:size(Ci1_bin,1)
%         origcode(i)=i;
%         missing=setdiff(unique(Ci2),targetcode);
%         targetcode(i)=missing(1);
%         
%     end
% end
for i=1:length(repvec)
    targetcode(find(targetcode==repvec(i)))=[];
end

for i=1:length(repvec)
    origcode(find(origcode==repvec(i)))=[];
end

% Iterative relabeling
for a=1:length(origcode)+1
    %if ~ismember(a,repvec)
        if a==1
            from(a) = origcode(a);
            to(a)   = dummyval;
            to(a+1) = from(a);
        else if a==length(origcode)+1
                from(a) = dummyval;
                to(a)   = targetcode(1);
            else    
                from(a) = origcode(find(targetcode==to(a)));
                to(a+1) = from(a); 
            end
        end
        vec_recoded(vec_recoded==from(a))=to(a);
%     else
%         from(a)= origcode(a);
%         to(a+1)=origcode(a+1);
%    end
end

end