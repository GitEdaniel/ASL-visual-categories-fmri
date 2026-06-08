% deafMovementBootstrap.m
%
% This will produce a bootstrapped contrast map with replacement
% and then compute distance and angle of movement from the grill-spector
% atlas ROI. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% List of subs in each group
subs_c = {"", "", "", "", "", "", "", "", "", "", ""};
subs_d = {"", "", "", "", "", "", "", "", "", "", ""};
subs_h = {"", "", "", "", "", "", "", "", "", "", ""};

% Matrix in which we'll store label distances to reference label
distances_DtoC = nan(1, strap_num);
distances_CtoC = nan(1, strap_num);
distances_HtoC = nan(1, strap_num);

% Vector in which we'll store angle formed between C1->D/H and MFS1->MFS2 as
% well as C1->C2 and MFS
movement_angles_D_MFS = nan(strap_num,1);
movement_angles_C_MFS = nan(strap_num,1);
movement_angles_H_MFS = nan(strap_num,1);

% MRIread does not seem to accept a full path, so let's cd
cd(boot_dir); 


% Bootstrap loop 
for b = 1:strap_num
    % Draw with replacement 11 deaf, 11 control group1, and 11 control 2
    boot_d = subs_d(randi(11,[1,11]));
    boot_c1 = subs_c(randi(11,[1,11]));
    boot_h = subs_h(randi(11,[1,11]));

    % Perform mri_concat and save each group's mean map to boot_dir
    cmd_c1 = strcat("!mri_concat --o ", boot_dir, "/c1.mgh --mean --i ");
    for sub=1:11
        cmd_c1 = strcat(cmd_c1, surface_path, "/", contrast, "/C/", boot_c1{sub}, "_fsavg_", contrast, "_", hemi, ".mgh ");
    end
    eval(cmd_c1); clear cmd_c1;

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
    cmd = strcat("!mri_surf2surf --sval ", boot_dir, "/d.mgh --hemi ", hemi, " --fwhm 7 --s fsaverage --tval ", boot_dir, "/d_sm.mgh");
    eval(cmd); clear cmd;
    cmd = strcat("!mri_surf2surf --sval ", boot_dir, "/h.mgh --hemi ", hemi, " --fwhm 7 --s fsaverage --tval ", boot_dir, "/h_sm.mgh");
    eval(cmd); clear cmd;


    % Load each overlay and mask data
    mri = MRIread('c1_sm.mgh'); contrast_c1 = mri.vol; clear mri;
    mri = MRIread('d_sm.mgh'); contrast_d = mri.vol; clear mri;
    mri = MRIread('h_sm.mgh'); contrast_h = mri.vol; clear mri;

    [label_vertices, ~] = read_label([], mask_label_path);
    label_vertices = label_vertices(:,1) + 1;
    contrast_mask_c1 = nan(size(contrast_c1)); contrast_mask_c1(label_vertices) = contrast_c1(label_vertices);
    contrast_mask_d = nan(size(contrast_d)); contrast_mask_d(label_vertices) = contrast_d(label_vertices);
    contrast_mask_h = nan(size(contrast_h)); contrast_mask_h(label_vertices) = contrast_h(label_vertices);
    clear contrast_c1 contrast_c2 contrast_d contrast_h;

    % Identify the threshold percentile vertices in each group
    percentile_c1 = prctile(contrast_mask_c1(label_vertices), threshold);
    percentile_d = prctile(contrast_mask_d(label_vertices), threshold);
    percentile_h = prctile(contrast_mask_h(label_vertices), threshold);
    verts_c1 = find(contrast_mask_c1>=percentile_c1);
    verts_d = find(contrast_mask_d>=percentile_d);
    verts_h = find(contrast_mask_h>=percentile_h);
    clear percentile_c1 percentile_c2 percentile_d percentile_h

    % Pass vertex indices to calculateSurfaceLabelMedoid
    % The verts being entered into the medoid function are 1-indexed
    surf_reference_path = fullfile(fsdir,'fsaverage','surf',strcat(hemi,".white"));
    medoid_c1 = calculateSurfaceLabelMedoid(surf_reference_path, verts_c1);
    medoid_d = calculateSurfaceLabelMedoid(surf_reference_path, verts_d);
    medoid_h = calculateSurfaceLabelMedoid(surf_reference_path, verts_h);

    % Load reference label and get medoid of that label
if hand
    % label_path = '/Applications/freesurfer/7.2.0/subjects/fsaverage/label/MPM_lh_OTS.label';
    % Label has holes in it which medoid code doesn't like so we derived
    % centroid from Freeview manually
    medoid_ref = 101953;
else
    % label_path = '/Applications/freesurfer/7.2.0/subjects/fsaverage/label/MPM_rh_mFus.label';
    medoid_ref = 118450;
end

    surf_for_dist = fullfile(fsdir,"fsaverage","surf",strcat(hemi,".white"));
    [path, pathlength_d] = calculateVertexDistance(surf_for_dist, medoid_ref, medoid_d);
    [path, pathlength_c1] = calculateVertexDistance(surf_for_dist, medoid_ref, medoid_c1);
    [path, pathlength_h] = calculateVertexDistance(surf_for_dist, medoid_ref, medoid_h);

    distances_DtoC(b) = pathlength_d;
    distances_HtoC(b) = pathlength_h;
    distances_CtoC(b) = pathlength_c1;

    % Let's get the angle of movement to compare how the deaf ROI moved
    % relative to controls. We will use the vector formed by the posterior
    % and anterior points of the MFS as a relative reference, so the angles
    % formed by C1->D and C1->C2 are the angle between these two vectors
    % and the anterior-pointing MFS vector
    
    if hand
        refVector_vert1 = 98939; % posterior MFS on fsaverage brain
        refVector_vert2 = 41735; % anterior MFS on fsaverage brain
    else
        refVector_vert1 = 40409; % posterior MFS on fsaverage brain
        refVector_vert2 = 13297; % anterior MFS on fsaverage brain
    end
        surf_reference_path_angle = fullfile(fsdir,'fsaverage','surf',strcat(hemi,".inflated"));
        
        % Angle of two vectors relative to MFS1->MFS2
        [angle_D_MFS, angle_C1C2_MFS] = calculateMovementAngle_MFSref(surf_reference_path_angle, medoid_ref, medoid_d, medoid_c1, refVector_vert1, refVector_vert2);
        movement_angles_D_MFS(b) = angle_D_MFS;
        movement_angles_C_MFS(b) = angle_C1C2_MFS;

        % Now for hearing
        % Angles of two vectors relative to MFS1->MFS2
        [angle_H_MFS, angle_C1C2_MFS] = calculateMovementAngle_MFSref(surf_reference_path_angle, medoid_ref, medoid_h, medoid_c1, refVector_vert1, refVector_vert2);
        movement_angles_H_MFS(b) = angle_H_MFS;

    % Wipe the boot_dir to save space
    !rm ./*.mgh
    
    progressbar(b/strap_num)
end


%% Save the data!
cd(code_dir)
if hand
    save('distances_andAngles_refLabelBootstrap_handOTS.mat','distances_DtoC','distances_CtoC','distances_HtoC', 'movement_angles_D_MFS','movement_angles_C_MFS','movement_angles_H_MFS')
else
    save('distances_andAngles_refLabelBootstrap_faceMFUS.mat','distances_DtoC','distances_CtoC','distances_HtoC', 'movement_angles_D_MFS','movement_angles_C_MFS','movement_angles_H_MFS')
end
    
% To visualize, open deafMovementQuantify_refLabelBootstrapped.m 