function threshold_connectivity_mat(in_file, out_dir)

% Function to threshold 1% structural connectivity matrices
%
% Inputs:       
%               - in_file: file with connectivity matrices normalized (
%               3D array: subjects' connectivity matrices concatenated in 
%               the 3rd dimension)
%
%               - out_dir: directory to save thresholded matrices 
%               (3D array)
%
%
% Ana Coelho  19-04-2021
%


% load matrices
M_file = load(in_file);
name = fieldnames(M_file);
M = M_file.(name{1,1});

% normalize each individual matrix to values [0,1]
M_norm={};

for i=1:size(M,3)
    M_norm{i,1}=M(:,:,i)./max(max(M(:,:,i)));
end

% create cell to save thresholded matrices
M_thr=cell(51,1);
thr = 0.01; % since matrices are already normalized to [0,1] we use this threshold that corresponds to 1%
% threshold matrices
for i=1:size(M_norm,1)
    
    % copy matrix to thresholded matrices
    M_thr{i,1}=M_norm{i,1};
    % find matrix elements below the threshold
    ids = find(M_norm{i,1}<thr);
    
    % set to 0 elements below the threshold
    for j=1:size(ids,1)
        M_thr{i,1}(ids(j))=0;
    end
end

% convert cell to 3D array
m_thr_norm = cat(3,M_thr{1,1},M_thr{2,1});
for i=3:size(M_thr,1)
    m_thr_norm = cat(3,m_thr_norm,M_thr{i,1});
end

if ~exist(out_dir, 'dir')
       mkdir(out_dir)
end
save(strcat(out_dir,'/matrices_norm_thr.mat'),'m_thr_norm');

end