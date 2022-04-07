function normalize_mats(wd,subj_list, atlas, out_folder)

% Function to normalize (by waytotal) and make symmetric structural 
% connectivity matrices
%
% Inputs:       
%               - wd: working directory with probtrackx outputs
%
%               - subj_list: list with subjects IDs
%
%               - atlas: name of atlas used (available options: aal, dkt40, shen)
%
%               - out_folder: output folder to save 3D matrix with 
%               normalized and thresholded matrices of all subjects
%
%
% Ana Coelho  19-04-2021
%


clear Q master
clear subjID master

% Open file with subjects ids
flogs = fopen(subj_list);
logno = textscan(flogs,'%s');
fclose(flogs);

% Create array with paths to subjects folders
for index_log=1:size(logno{1},1)
    subj = cell2mat(logno{1}(index_log));
    subj = [wd '/' subj];
    Q{index_log,:} = subj;
end

switch atlas
    case 'aal'
        nrois=90;
    case 'dkt40'
        nrois=76;
    case 'shen'
        nrois=268;
end

%Load network matrix for all subjects, normalize and make symmetric
M = cell(numel(Q),1);

for master = 1:numel(Q)
    %load network matrix
    mat = load([deblank(Q{master,:}),'/fdt_network_matrix']);
    
    % load waytotal file
    waytotal = load([deblank(Q{master,:}),'/waytotal']);
    
    %normalize matrix
    mat_norm = zeros(nrois);
    for i=1:nrois
        for j=1:nrois
            mat_norm(i,j) = mat(i,j)/waytotal(i);
        end
    end
    
    %create symmetric matrix
    sym_mat_norm = zeros(nrois,nrois);
    for i=1:nrois
        for j=1:nrois
            sym_mat_norm(i,j) = (mat_norm(i,j) + mat_norm(j,i))/2;
        end
    end
    
    %save to cell
    M{master} = sym_mat_norm;
    
end

% convert cell to 3D array
m_norm = cat(3,M{1,1},M{2,1});
for i=3:size(M,1)
    m_norm = cat(3,m_norm,M{i,1});
end

if ~exist(out_folder, 'dir')
       mkdir(out_folder)
end
save(strcat(out_folder,'/matrices_norm.mat'),'m_norm');

end
