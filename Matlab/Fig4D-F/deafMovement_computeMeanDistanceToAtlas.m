% deafMovementBootstrap.m
%
% We will load each group, then average the
% contrast maps, smooth them, identify the peak in each group, and
% calculate the distance between those two vertices. We will also bootstrap
% a second control group (drawing with replacement from the controls
% again) another smoothe map to see how far that peak is from the first
% control peak (e.g., trying to estiamte how stable the control peak is
% relative to the distance this peak is to the deaf group). 
%
% This calls a function (calculateSurfaceLabelMedoid.m) that
% will find the "center" vertex within a label. It will accept the full
% path to a label or the column vector of the label indices in a label. Be
% sure the label is from the same surface/overlay you are loading. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataDir = '';
fsdir = getenv('SUBJECTS_DIR');
% Add path to where this code and calculateSurfaceLabelMedoid.m live
code_dir = '';
addpath(genpath(code_dir))
addpath(genpath('/Applications/freesurfer/7.2.0/matlab'))


%%%%%%%%%%%%%%%% Variables to define %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Matrix in which we'll store label distances for individuals to atlas ROI
distances_D = nan(1, 11);
distances_C = nan(1, 11);
distances_H = nan(1, 11);


% Vector in which we'll store angle formed between outside label and each
% group using MFS to give us signed reference angle
movement_angles_D_MFS = nan(11,1);
movement_angles_C_MFS = nan(11,1);
movement_angles_H_MFS = nan(11,1);

% Same as above for the mean ROI derived in each group
movement_angleMean_D_MFS = [];
movement_angleMean_C_MFS = [];
movement_angleMean_H_MFS = [];


% MRIread does not seem to accept a full path, so let's cd
cd(boot_dir); 


% Perform mri_concat and save each group's mean map to boot_dir
cmd_c1 = strcat("!mri_concat --o ", boot_dir, "/c1.mgh --mean --i ");
for sub=1:11
    cmd_c1 = strcat(cmd_c1, surface_path, "/", contrast, "/C/", subs_c{sub}, "_fsavg_", contrast, "_", hemi, ".mgh ");
end
eval(cmd_c1); clear cmd_c1;

cmd_d = strcat("!mri_concat --o ", boot_dir, "/d.mgh --mean --i ");
for sub=1:11
    cmd_d = strcat(cmd_d, surface_path, "/", contrast, "/D/", subs_d{sub}, "_fsavg_", contrast, "_", hemi, ".mgh ");
end
eval(cmd_d); clear cmd_d;

cmd_h = strcat("!mri_concat --o ", boot_dir, "/h.mgh --mean --i ");
for sub=1:11
    cmd_h = strcat(cmd_h, surface_path, "/", contrast, "/H/", subs_h{sub}, "_fsavg_", contrast, "_", hemi, ".mgh ");
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
clear contrast_c1 contrast_d contrast_h;

% Identify the threshold percentile vertices in each group
percentile_c1 = prctile(contrast_mask_c1(label_vertices), threshold);
percentile_d = prctile(contrast_mask_d(label_vertices), threshold);
percentile_h = prctile(contrast_mask_h(label_vertices), threshold);
verts_c1 = find(contrast_mask_c1>=percentile_c1);
verts_d = find(contrast_mask_d>=percentile_d);
verts_h = find(contrast_mask_h>=percentile_h);
clear percentile_c1 percentile_d percentile_h

% Pass vertex indices to calculateSurfaceLabelMedoid
% The verts being entered into the medoid function are 1-indexed
surf_reference_path = fullfile(fsdir,'fsaverage','surf',strcat(hemi,".white"));
medoid_c1 = calculateSurfaceLabelMedoid(surf_reference_path, verts_c1);
medoid_d = calculateSurfaceLabelMedoid(surf_reference_path, verts_d);
medoid_h = calculateSurfaceLabelMedoid(surf_reference_path, verts_h);

% Load reference label and get medoid of that label
if hand
    % Label has holes in it which medoid code doesn't like so we derived
    % Rosenke atlas centroid from Freeview 7.2.0 manually in each roi
    medoid_ref = 101953;
else
    medoid_ref = 118450;
end

% Get distance between ref and each group's ROI by calling the
% calculateVertexDistance.m function which expects the path to a
% surface (we'll use the white surface), a start vertex on that
% surface, and an end vertex
surf_for_dist = fullfile(fsdir,"fsaverage","surf",strcat(hemi,".white"));
[path, pathlength_d] = calculateVertexDistance(surf_for_dist, medoid_ref, medoid_d);
[path, pathlength_c1] = calculateVertexDistance(surf_for_dist, medoid_ref, medoid_c1);
[path, pathlength_h] = calculateVertexDistance(surf_for_dist, medoid_ref, medoid_h);
pathlength_mean_d = pathlength_d;
pathlength_mean_c = pathlength_c1;
pathlength_mean_h = pathlength_h; 

% Let's get the angle of movement 
if hand
    refVector_vert1 = 98939; % posterior MFS on fsaverage brain
    refVector_vert2 = 41735; % anterior MFS on fsaverage brain
else
    refVector_vert1 = 40409; % posterior MFS on fsaverage brain
    refVector_vert2 = 13297; % anterior MFS on fsaverage brain
end
surf_reference_path_angle = fullfile(fsdir,'fsaverage','surf',strcat(hemi,".inflated"));

% First we'll do the D and C groups (H is next section)
% Just angle between ref->D and ref->C relative to MFS1->MFS2
[angle_D_MFS, angle_C1_MFS] = calculateMovementAngle_MFSref(surf_reference_path_angle, medoid_ref, medoid_d, medoid_c1, refVector_vert1, refVector_vert2);
movement_angleMean_D_MFS = angle_D_MFS;
movement_angleMean_C_MFS = angle_C1_MFS;

% Now for hearing
[angle_H_MFS, angle_C1_MFS] = calculateMovementAngle_MFSref(surf_reference_path_angle, medoid_ref, medoid_h, medoid_c1, refVector_vert1, refVector_vert2);
movement_angleMean_H_MFS = angle_H_MFS;

clear medoid_c1 medoid_d medoid_h angle_D_MFS angle_C1_MFS angle_H_MFS

%% Save the data!
cd(code_dir)
if hand
    save('distances_andAngles_reflabel_handOTS.mat', 'movement_angleMean_C_MFS','movement_angleMean_D_MFS','movement_angleMean_H_MFS','pathlength_mean_d','pathlength_mean_c','pathlength_mean_h')
else
    save('distances_andAngles_reflabel_faceMFUS.mat', 'movement_angleMean_C_MFS','movement_angleMean_D_MFS','movement_angleMean_H_MFS','pathlength_mean_d','pathlength_mean_c','pathlength_mean_h')
end
% To visualize, open deafMovementQuantify_refLabel.m 