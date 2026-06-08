% deafMovementBootstrap.m
%
% We will bootstrap with replacement 11 brains from each group, average the
% contrast maps, smooth them, identify the peak in each group, and
% calculate the distance between those two vertices. We will also bootstrap
% a second control group (drawing with replacement from the controls
% again) another smoothed map to see how far that peak is from the first
% control peak (e.g., trying to estiamte how stable the control peak is
% relative to the distance this peak is to the deaf group). 
%
% This calls a function (calculateSurfaceLabelMedoid.m) that
% will find the "center" vertex within a label. It will accept the full
% path to a label or the column vector of the label indices in a label. Be
% sure the label is from the same surface/overlay you are loading. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define path where category contrast surface overlays reside
dataDir = '';
fsdir = getenv('SUBJECTS_DIR');

% Add path to where this code and calculateSurfaceLabelMedoid.m live
code_dir = '';
addpath(genpath(code_dir))
addpath(genpath('/Applications/freesurfer/7.2.0/matlab'))


%%%%%%%%%%%%%%%% Variables to define %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of bootstraps
strap_num = 500;

% Directory of the fsaverage-space single subject contrast overlays
% Organized as /single_sub/contrast/D or C or H
surface_path = "";

% Which contrast?
hand = false;

if hand
    contrast = "hand-v-nonlimb";
    threshold = 90; % Use for hand contrast
    hemi='lh';
    mask_label_path = fullfile(fsdir,"fsaverage","label","lh.hand_movement_mask.label"); 
else
    contrast = "face-v-all";
    threshold = 80; % Use for face contrast (lowered a little to ensure a single cluster
    hemi='rh';
    mask_label_path = fullfile(fsdir,"fsaverage","label","rh.mFus_movement_mask.label");
end


% Directory where temp files will be written during bootstrapping
boot_dir = "";


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List of subs in each groups
subs_c = {"", "", "", "", "", "", "", "", "", "", ""};
subs_d = {"", "", "", "", "", "", "", "", "", "", ""};
subs_h = {"", "", "", "", "", "", "", "", "", "", ""};

% Matrix in which we'll store label distances for deaf (D) and control (C)
distances_DtoC = nan(1, strap_num);
distances_CtoC = nan(1, strap_num);
distances_HtoC = nan(1, strap_num);
distances_DtoH = nan(1, strap_num);

% Vector in which we'll store angle formed between C1->D/H and MFS1->MFS2 as
% well as C1->C2 and MFS
movement_angles_D_MFS = nan(strap_num,1);
movement_angles_C_MFS = nan(strap_num,1);
movement_angles_H_MFS = nan(strap_num,1);

% We can also keep track of whether D is closer to C1 than C2 which would
% imply that the C2 location is more posterior than C1 and therefore moved
% in the opposite direction that D does in which case we can make this
% pathlength (C1 to C2) negative. 
distances_CtoC_neg = nan(1, strap_num);
distances_HtoC_neg = nan(1, strap_num);

% MRIread does not seem to accept a full path, so we must change dir
cd(boot_dir); 

% Bootstrap loop 

for b = 1:strap_num
    % Draw with replacement 11 deaf, 11 control group1, and 11 control 2
    boot_d = subs_d(randi(11,[1,11]));
    boot_c1 = subs_c(randi(11,[1,11]));
    boot_c2 = subs_c(randi(11,[1,11]));
    boot_h = subs_h(randi(11,[1,11]));

    % Perform mri_concat and save each group's mean map to boot_dir
    cmd_c1 = strcat("!mri_concat --o ", boot_dir, "/c1.mgh --mean --i ");
    for sub=1:11
        cmd_c1 = strcat(cmd_c1, surface_path, "/", contrast, "/C/", boot_c1{sub}, "_fsavg_", contrast, "_", hemi, ".mgh ");
    end
    eval(cmd_c1); clear cmd_c1;

    cmd_c2 = strcat("!mri_concat --o ", boot_dir, "/c2.mgh --mean --i ");
    for sub=1:11
        cmd_c2 = strcat(cmd_c2, surface_path, "/", contrast, "/C/", boot_c2{sub}, "_fsavg_", contrast, "_", hemi, ".mgh ");
    end
    eval(cmd_c2); clear cmd_c2;

    cmd_d = strcat("!mri_concat --o ", boot_dir, "/d.mgh --mean --i ");
    for sub=1:11
        cmd_d = strcat(cmd_d, surface_path, "/", contrast, "/D/", boot_d{sub}, "_fsavg_", contrast, "_", hemi, ".mgh ");
    end
    eval(cmd_d); clear cmd_d;

    cmd_h = strcat("!mri_concat --o ", boot_dir, "/h.mgh --mean --i ");
    for sub=1:11
        cmd_h = strcat(cmd_h, surface_path, "/", contrast, "/H/", boot_h{sub}, "_fsavg_", contrast, "_", hemi, ".mgh ");
    end
    eval(cmd_h); clear cmd_h;

    % Smooth each overlay map with mri_surf2surf
    cmd = strcat("!mri_surf2surf --sval ", boot_dir, "/c1.mgh --hemi ", hemi, " --fwhm 7 --s fsaverage --tval ", boot_dir, "/c1_sm.mgh");
    eval(cmd); clear cmd;
    cmd = strcat("!mri_surf2surf --sval ", boot_dir, "/c2.mgh --hemi ", hemi, " --fwhm 7 --s fsaverage --tval ", boot_dir, "/c2_sm.mgh");
    eval(cmd); clear cmd;
    cmd = strcat("!mri_surf2surf --sval ", boot_dir, "/d.mgh --hemi ", hemi, " --fwhm 7 --s fsaverage --tval ", boot_dir, "/d_sm.mgh");
    eval(cmd); clear cmd;
    cmd = strcat("!mri_surf2surf --sval ", boot_dir, "/h.mgh --hemi ", hemi, " --fwhm 7 --s fsaverage --tval ", boot_dir, "/h_sm.mgh");
    eval(cmd); clear cmd;

    % Load each overlay and mask data
    mri = MRIread('c1_sm.mgh'); contrast_c1 = mri.vol; clear mri;
    mri = MRIread('c2_sm.mgh'); contrast_c2 = mri.vol; clear mri;
    mri = MRIread('d_sm.mgh'); contrast_d = mri.vol; clear mri;
    mri = MRIread('h_sm.mgh'); contrast_h = mri.vol; clear mri;

    [label_vertices, ~] = read_label([], mask_label_path);
    label_vertices = label_vertices(:,1) + 1;
    contrast_mask_c1 = nan(size(contrast_c1)); contrast_mask_c1(label_vertices) = contrast_c1(label_vertices);
    contrast_mask_c2 = nan(size(contrast_c2)); contrast_mask_c2(label_vertices) = contrast_c2(label_vertices);
    contrast_mask_d = nan(size(contrast_d)); contrast_mask_d(label_vertices) = contrast_d(label_vertices);
    contrast_mask_h = nan(size(contrast_h)); contrast_mask_h(label_vertices) = contrast_h(label_vertices);
    clear contrast_c1 contrast_c2 contrast_d contrast_h;

    % Identify the threshold percentile vertices in each group
    percentile_c1 = prctile(contrast_mask_c1(label_vertices), threshold);
    percentile_c2 = prctile(contrast_mask_c2(label_vertices), threshold);
    percentile_d = prctile(contrast_mask_d(label_vertices), threshold);
    percentile_h = prctile(contrast_mask_h(label_vertices), threshold);
    verts_c1 = find(contrast_mask_c1>=percentile_c1);
    verts_c2 = find(contrast_mask_c2>=percentile_c2);
    verts_d = find(contrast_mask_d>=percentile_d);
    verts_h = find(contrast_mask_h>=percentile_h);
    clear percentile_c1 percentile_c2 percentile_d percentile_h

    % Pass vertex indices to calculateSurfaceLabelMedoid
    % The verts being entered into the medoid function are 1-indexed
    surf_reference_path = fullfile(fsdir,'fsaverage','surf',strcat(hemi,".white"));
    medoid_c1 = calculateSurfaceLabelMedoid(surf_reference_path, verts_c1);
    medoid_c2 = calculateSurfaceLabelMedoid(surf_reference_path, verts_c2);
    medoid_d = calculateSurfaceLabelMedoid(surf_reference_path, verts_d);
    medoid_h = calculateSurfaceLabelMedoid(surf_reference_path, verts_h);

    % Get distance between D & C1, and C1 & C2 by calling the
    % calculateVertexDistance.m function which expects the path to a
    % surface (we'll use the white surface), a start vertex on that
    % surface, and an end vertex
    surf_for_dist = fullfile(fsdir,"fsaverage","surf",strcat(hemi,".white"));
    [path, pathlength_dc1] = calculateVertexDistance(surf_for_dist, medoid_c1, medoid_d);
    [path, pathlength_c1c2] = calculateVertexDistance(surf_for_dist, medoid_c1, medoid_c2);
    [path, pathlength_dc2] = calculateVertexDistance(surf_for_dist, medoid_d, medoid_c2);
    [path, pathlength_hc1] = calculateVertexDistance(surf_for_dist, medoid_c1, medoid_h);
    [path, pathlength_dh] = calculateVertexDistance(surf_for_dist, medoid_d, medoid_h);

    distances_DtoC(b) = pathlength_dc1;
    distances_HtoC(b) = pathlength_hc1;
    distances_CtoC(b) = pathlength_c1c2;
    distances_DtoH(b) = pathlength_dh;
    distances_CtoC_neg(b) = pathlength_c1c2; % may flip negative below
    distances_HtoC_neg(b) = pathlength_hc1; % may flip negative below


    if pathlength_dc1 < pathlength_dc2
        % If C2 location is on other side of C1 we'll count this as a
        % negative distance instead. Otherwise it's left positive
        distances_CtoC_neg(b) = -1 * pathlength_c1c2;
    end

    % And if D->C1 & H->C1 are less than D->H then H moved away from C1 in
    % a direction distinct from D. Otherwise, if it moved towards the D
    % coordinate then D->H would be less than D->C1
    if mean([pathlength_dc1 pathlength_hc1]) < pathlength_dh
        distances_HtoC_neg(b) = -1 * pathlength_hc1;
    end

    % Let's get the angle of movement to compare how the deaf ROI moved
    % relative to controls. We will use the vector formed by the posterior
    % and anterior points of the MFS as a relative reference, so the angles
    % formed by C1->D and C1->C2 are the angle between these two vectors
    % and the anterior-pointing MFS vector

    if strmatch('rh',hemi)
        refVector_vert1 = 40409; % posterior MFS on fsaverage brain
        refVector_vert2 = 13297; % anterior MFS on fsaverage brain
    elseif strmatch('lh',hemi)
        refVector_vert1 = 98939; % posterior MFS on fsaverage brain
        refVector_vert2 = 41735; % anterior MFS on fsaverage brain
    end

    surf_reference_path_angle = fullfile(fsdir,'fsaverage','surf',strcat(hemi,".inflated"));

    % Angles of  vectors relative to MFS1->MFS2
    [angle_D_MFS, angle_C1C2_MFS] = calculateMovementAngle_MFSref(surf_reference_path_angle, medoid_c1, medoid_d, medoid_c2, refVector_vert1, refVector_vert2);
    movement_angles_D_MFS(b) = angle_D_MFS;
    movement_angles_C_MFS(b) = angle_C1C2_MFS;
    clear angle_C1C2_MFS angle_D_MFS
    % Now for hearing
    [angle_H_MFS, angle_C1C2_MFS] = calculateMovementAngle_MFSref(surf_reference_path_angle, medoid_c1, medoid_h, medoid_c2, refVector_vert1, refVector_vert2);
    movement_angles_H_MFS(b) = angle_H_MFS;
    clear angle_H_MFS

    % Wipe the boot_dir to save space (be careful, the mgh part is important!)
    %!rm ./*.mgh
    
end


%% Save the data!
cd(code_dir)
if hand
    save('distances_bootstrap_values_handOTS_neg_andAngles.mat','distances_DtoC','distances_CtoC','distances_HtoC','distances_DtoH','distances_CtoC_neg','distances_HtoC_neg', 'movement_angles_D_MFS','movement_angles_C_MFS','movement_angles_H_MFS')
else
    save('distances_bootstrap_values_faceMFUS_neg_andAngles.mat','distances_DtoC','distances_CtoC','distances_HtoC','distances_DtoH','distances_CtoC_neg','distances_HtoC_neg', 'movement_angles_D_MFS','movement_angles_C_MFS','movement_angles_H_MFS')
end

% To visualize, open deafMovementQuantify.m 