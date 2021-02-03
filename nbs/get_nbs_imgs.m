function get_nbs_imgs(cmp, out_dir, atlas)

% Function to create images from NBS results
% The following images will be generated:
%   monocrome images in axial plane of each sub-network
%   sagittal, axial and coronal images of each sub-network with decreases
%   in blue and increases in red
%
% Inputs:
%           cmp: component resulting from NBS analysis
%           out_dir: directory to save output files
%           atlas: atlas to be used (1 for shen atlas, 2 for aal, 3 for myatlas, 4 for shen with no cerebellum, 5 for silhouette parcellation)
%
% Ana Coelho 05-02-2020

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

% Monocrome images
for i=1:size(cmp,2)
    figure('units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto','color','white'); hold on;
    p=patch('vertices',cortexvert,'faces',cortexface,'EdgeColor','none','FaceColor',[0.4 0.4 .4],'FaceAlpha',0.08); %plots a mesh surface for the cortex
    material dull;
    colormap HSV
    cm_cmp     = colormap;
    colour_arr = [49 37 7 19 43 55 13 1 62];
    s_factor   = 3;
    monocrom   = 1;
    cmp{1,i}.colour=cm_cmp(colour_arr(4),:); %green
    
    for l=1:length(cmp{1,i}.area1.number)
        if monocrom
            edge_colour = cmp{1,i}.colour;
        else
            if cmp{1,i}.wei_norm(cmp{1,i}.area1.number(l),cmp{1,i}.area2.number(l))<0
                edge_colour = [0 0.34 1];  % blue; decrease in connectivity; SC1>SC2
            else
                edge_colour = [1 0.18 0];  % red; increase in connectivity; SC1<SC2
            end
        end
        edge_width  = cmp{1,i}.wei_norm(cmp{1,i}.area1.number(l),cmp{1,i}.area2.number(l))*s_factor;
        [X1,Y1,Z1]  = cylinder1(MNI_coord(cmp{1,i}.area1.number(l),:),MNI_coord(cmp{1,i}.area2.number(l),:),edge_width,100);
        wei_conn_lr = surf(X1,Y1,Z1,'FaceColor',edge_colour,'EdgeColor','none');
        sphere2(MNI_coord(cmp{1,i}.area1.number(l),1),MNI_coord(cmp{1,i}.area1.number(l),2),MNI_coord(cmp{1,i}.area1.number(l),3),edge_width,edge_colour,1); %text(MNI_coord(i,1)-2,MNI_coord(i,2)-3,MNI_coord(i,3),int2str(i),'FontWeight','bold','FontSize',12);
        sphere2(MNI_coord(cmp{1,i}.area2.number(l),1),MNI_coord(cmp{1,i}.area2.number(l),2),MNI_coord(cmp{1,i}.area2.number(l),3),edge_width,edge_colour,1); %text(MNI_coord(j,1)-2,MNI_coord(j,2)-6,MNI_coord(j,3),int2str(j),'FontWeight','bold','FontSize',12);
    end
    axis equal
    axis off
    camlight;
    print(gcf,strcat(out_dir,'/net',num2str(i),'monocrome_axial.png'),'-dpng','-r300');
end

% Multicolor images
% Axial
for i=1:size(cmp,2)
    figure('units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto','color','white'); hold on;
    p=patch('vertices',cortexvert,'faces',cortexface,'EdgeColor','none','FaceColor',[0.4 0.4 .4],'FaceAlpha',0.08); %plots a mesh surface for the cortex
    material dull;
    colormap HSV
    cm_cmp     = colormap;
    colour_arr = [49 37 7 19 43 55 13 1 62];
    s_factor   = 3;
    monocrom   = 0;
    
    for l=1:length(cmp{1,i}.area1.number)
        if monocrom
            edge_colour = cmp{1,i}.colour;
        else
            if cmp{1,i}.wei_norm(cmp{1,i}.area1.number(l),cmp{1,i}.area2.number(l))<0
                edge_colour = [0 0.34 1];  % blue; decrease in connectivity; SC1>SC2
            else
                edge_colour = [1 0.18 0];  % red; increase in connectivity; SC1<SC2
            end
        end
        edge_width  = cmp{1,i}.wei_norm(cmp{1,i}.area1.number(l),cmp{1,i}.area2.number(l))*s_factor;
        [X1,Y1,Z1]  = cylinder1(MNI_coord(cmp{1,i}.area1.number(l),:),MNI_coord(cmp{1,i}.area2.number(l),:),edge_width,100);
        wei_conn_lr = surf(X1,Y1,Z1,'FaceColor',edge_colour,'EdgeColor','none');
        sphere2(MNI_coord(cmp{1,i}.area1.number(l),1),MNI_coord(cmp{1,i}.area1.number(l),2),MNI_coord(cmp{1,i}.area1.number(l),3),edge_width,edge_colour,1); %text(MNI_coord(i,1)-2,MNI_coord(i,2)-3,MNI_coord(i,3),int2str(i),'FontWeight','bold','FontSize',12);
        sphere2(MNI_coord(cmp{1,i}.area2.number(l),1),MNI_coord(cmp{1,i}.area2.number(l),2),MNI_coord(cmp{1,i}.area2.number(l),3),edge_width,edge_colour,1); %text(MNI_coord(j,1)-2,MNI_coord(j,2)-6,MNI_coord(j,3),int2str(j),'FontWeight','bold','FontSize',12);
    end
    axis equal
    axis off
    camlight;
    print(gcf,strcat(out_dir,'/net',num2str(i),'_axial.png'),'-dpng','-r300');
    
end

% Coronal
for i=1:size(cmp,2)
    figure('units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto','color','white'); hold on;
    p=patch('vertices',cortexvert,'faces',cortexface,'EdgeColor','none','FaceColor',[0.4 0.4 .4],'FaceAlpha',0.08); %plots a mesh surface for the cortex
    material dull;
    colormap HSV
    cm_cmp     = colormap;
    colour_arr = [49 37 7 19 43 55 13 1 62];
    s_factor   = 3;
    monocrom   = 0;
    
    for l=1:length(cmp{1,i}.area1.number)
        if monocrom
            edge_colour = cmp{1,i}.colour;
        else
            if cmp{1,i}.wei_norm(cmp{1,i}.area1.number(l),cmp{1,i}.area2.number(l))<0
                edge_colour = [0 0.34 1];  % blue; decrease in connectivity; SC1>SC2
            else
                edge_colour = [1 0.18 0];  % red; increase in connectivity; SC1<SC2
            end
        end
        edge_width  = cmp{1,i}.wei_norm(cmp{1,i}.area1.number(l),cmp{1,i}.area2.number(l))*s_factor;
        [X1,Y1,Z1]  = cylinder1(MNI_coord(cmp{1,i}.area1.number(l),:),MNI_coord(cmp{1,i}.area2.number(l),:),edge_width,100);
        wei_conn_lr = surf(X1,Y1,Z1,'FaceColor',edge_colour,'EdgeColor','none');
        sphere2(MNI_coord(cmp{1,i}.area1.number(l),1),MNI_coord(cmp{1,i}.area1.number(l),2),MNI_coord(cmp{1,i}.area1.number(l),3),edge_width,edge_colour,1); %text(MNI_coord(i,1)-2,MNI_coord(i,2)-3,MNI_coord(i,3),int2str(i),'FontWeight','bold','FontSize',12);
        sphere2(MNI_coord(cmp{1,i}.area2.number(l),1),MNI_coord(cmp{1,i}.area2.number(l),2),MNI_coord(cmp{1,i}.area2.number(l),3),edge_width,edge_colour,1); %text(MNI_coord(j,1)-2,MNI_coord(j,2)-6,MNI_coord(j,3),int2str(j),'FontWeight','bold','FontSize',12);
    end
    
    object_handles = findall(gcf,'Type','surface');
    rotate(object_handles,[1 0 0],-90);hold on;
    rotate(p,[1 0 0],-90);hold on;

    axis equal
    axis off
    camlight;
    print(gcf,strcat(out_dir,'/net',num2str(i),'_coronal.png'),'-dpng','-r300');
    
end

% Sagittal
for i=1:size(cmp,2)
    figure('units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto','color','white'); hold on;
    p=patch('vertices',cortexvert,'faces',cortexface,'EdgeColor','none','FaceColor',[0.4 0.4 .4],'FaceAlpha',0.08); %plots a mesh surface for the cortex
    material dull;
    colormap HSV
    cm_cmp     = colormap;
    colour_arr = [49 37 7 19 43 55 13 1 62];
    s_factor   = 3;
    monocrom   = 0;
    
    for l=1:length(cmp{1,i}.area1.number)
        if monocrom
            edge_colour = cmp{1,i}.colour;
        else
            if cmp{1,i}.wei_norm(cmp{1,i}.area1.number(l),cmp{1,i}.area2.number(l))<0
                edge_colour = [0 0.34 1];  % blue; decrease in connectivity; SC1>SC2
            else
                edge_colour = [1 0.18 0];  % red; increase in connectivity; SC1<SC2
            end
        end
        edge_width  = cmp{1,i}.wei_norm(cmp{1,i}.area1.number(l),cmp{1,i}.area2.number(l))*s_factor;
        [X1,Y1,Z1]  = cylinder1(MNI_coord(cmp{1,i}.area1.number(l),:),MNI_coord(cmp{1,i}.area2.number(l),:),edge_width,100);
        wei_conn_lr = surf(X1,Y1,Z1,'FaceColor',edge_colour,'EdgeColor','none');
        sphere2(MNI_coord(cmp{1,i}.area1.number(l),1),MNI_coord(cmp{1,i}.area1.number(l),2),MNI_coord(cmp{1,i}.area1.number(l),3),edge_width,edge_colour,1); %text(MNI_coord(i,1)-2,MNI_coord(i,2)-3,MNI_coord(i,3),int2str(i),'FontWeight','bold','FontSize',12);
        sphere2(MNI_coord(cmp{1,i}.area2.number(l),1),MNI_coord(cmp{1,i}.area2.number(l),2),MNI_coord(cmp{1,i}.area2.number(l),3),edge_width,edge_colour,1); %text(MNI_coord(j,1)-2,MNI_coord(j,2)-6,MNI_coord(j,3),int2str(j),'FontWeight','bold','FontSize',12);
    end
    
    object_handles = findall(gcf,'Type','surface');
    rotate(object_handles,[0 1 0],90);hold on;
    rotate(object_handles,[0 0 1],90)
    rotate(p,[0 1 0],90);hold on;
    rotate(p,[0 0 1],90)

    axis equal
    axis off
    camlight;
    print(gcf,strcat(out_dir,'/net',num2str(i),'_sagittal.png'),'-dpng','-r300');
    
end

end