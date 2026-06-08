% deafMovementQuantify_refLabelBootstrapped.m
%
% This will load the results saved out from two scripts: (1)
% deafMovement_computeMeanDistanceToAtlas and (2)
% deafMovement_compareToAtlasBootstrap. The former derives the vector
% between the visfAtlas reference label and the mean face/hand ROI for each
% group. The latter bootstraps these face/hand ROIS to derive distributions
% from which we can calculate Cohen's D. So those two scripts will need
% to be run in order to visualize all the figures of this code. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data, where art thou
code_dir = '';

face = true;
hand = false;

%% Vectors visualizing displacement to visfAtlas
% Loads data from deafMovement_computeMeanDistanceToAtlas

if face
    load(fullfile(code_dir, 'distances_andAngles_reflabel_faceMFUS.mat'))
elseif hand
    load(fullfile(code_dir, 'distances_andAngles_reflabel_handOTS.mat'))
end

% Derive vector components (x,y) in each group
if hand
    x_d = pathlength_mean_d .* cosd(movement_angleMean_D_MFS+90);
    y_d = pathlength_mean_d .* sind(movement_angleMean_D_MFS+90);
    x_c = pathlength_mean_c .* cosd(movement_angleMean_C_MFS+90);
    y_c = pathlength_mean_c .* sind(movement_angleMean_C_MFS+90);
    x_h = pathlength_mean_h .* cosd(movement_angleMean_H_MFS+90);
    y_h = pathlength_mean_h .* sind(movement_angleMean_H_MFS+90);
else
    % Reverse rotation for other hemisphere
    x_d = pathlength_mean_d .* sind(movement_angleMean_D_MFS-90);
    y_d = pathlength_mean_d .* cosd(movement_angleMean_D_MFS-90);
    x_c = pathlength_mean_c .* sind(movement_angleMean_C_MFS-90);
    y_c = pathlength_mean_c .* cosd(movement_angleMean_C_MFS-90);
    x_h = pathlength_mean_h .* sind(movement_angleMean_H_MFS-90);
    y_h = pathlength_mean_h .* cosd(movement_angleMean_H_MFS-90);
end

% Plot arrows from origin
f=figure('color','w','position',[73 18 547 659]);
quiver(zeros(size(x_d)), zeros(size(y_d)), x_d, y_d, 0,'color',[0.8 0.2 0.4],'linewidth',5);
hold on;
quiver(zeros(size(x_c)), zeros(size(y_c)), x_c, y_c, 0,'color',[0.2 0.2 0.8],'linewidth',5);
quiver(zeros(size(x_h)), zeros(size(y_h)), x_h, y_h, 0,'color',[224, 191, 45]./255,'linewidth',5);
grid on;
set(gca,'box','off','fontsize',14,'tickdir','out','linewidth',2)
xlabel('Medial Distance [mm]','fontsize',16)
ylabel('Anterior-Posterior Distance [mm]','fontsize',16)
if face
    title('mFus reference shift','FontWeight','bold','fontsize',22)
else
    title('OTS reference shift','FontWeight','bold','fontsize',22)
end

if face
    exportgraphics(gcf, 'figure_mFus_vectorsMean_refLabel.pdf', 'ContentType', 'vector');
elseif hand
    exportgraphics(gcf, 'figure_OTS_vectorsMean_refLabel.pdf', 'ContentType', 'vector');
end


%% Histograms and Cohen's D
% Loads data from deafMovement_compareToAtlasBootstrap
if face
    load(fullfile(code_dir, 'distances_andAngles_refLabelBootstrap_faceMFUS.mat'))
elseif hand
    load(fullfile(code_dir, 'distances_andAngles_refLabelBootstrap_handOTS.mat'))
end

% Histograms

f2=figure('color','white','position',[67 1 476 553]);
% Compute kernel density for D to refLabel
[f, xi] = ksdensity(distances_DtoC);
hold on
% Filled density
h = fill(xi, f, [199, 123, 179]./255, ...
    'FaceAlpha', 0.35, ...
    'EdgeColor', [0 0 0], ...
    'LineWidth', 1.5);
% Compute kernel density for C to refLabel
[f, xi] = ksdensity(distances_CtoC);
hold on
% Filled density
h = fill(xi, f, [88, 132, 209]./255, ...
    'FaceAlpha', 0.35, ...
    'EdgeColor', [0 0 0], ...
    'LineWidth', 1.5);
% Compute kernel density for H to refLabel
[f, xi] = ksdensity(distances_HtoC);
hold on
% Filled density
h = fill(xi, f, [224, 191, 45]./255, ...
    'FaceAlpha', 0.35, ...
    'EdgeColor', [0 0 0], ...
    'LineWidth', 1.5);
% Clean styling
box off
set(gca, ...
    'FontName', 'Helvetica', ...
    'FontSize', 15, ...
    'LineWidth', 2)
xlabel('Distance [mm]','FontSize',15)
ylabel('Density','FontSize',15)
% Optional: remove top/right border
set(gca,'TickDir','out')

if hand
    effect = meanEffectSize(distances_DtoC(1:50),distances_HtoC(1:50),Effect="cohen");
    effect = meanEffectSize(distances_DtoC(1:50),distances_CtoC(1:50),Effect="cohen");
    text(7.98,0.275,'Deaf-Hearing CohensD = 1.01','fontsize',12)
    text(7.98,0.25,'Deaf-Control CohensD = 0.67','fontsize',12)
    text(11,0.22,'Deaf','FontSize',18,'FontWeight','bold','color',[199, 123, 179]./255)
    text(11,0.195,'Typical','FontSize',18,'FontWeight','bold','color',[88, 132, 209]./255)
    text(11,0.170,'Hearing','FontSize',18,'FontWeight','bold','color',[224, 191, 45]./255)
    title('OTS-hands shift','fontsize',22)
else
    effect = meanEffectSize(distances_DtoC(1:50),distances_HtoC(1:50),Effect="cohen");
    effect = meanEffectSize(distances_DtoC(1:50),distances_CtoC(1:50),Effect="cohen");
    text(7,0.6,'Deaf-Hearing CohensD = 1.6','fontsize',12)
    text(7,0.55,'Deaf-Control CohensD = 1.23','fontsize',12)
    text(11,0.5,'Deaf','FontSize',18,'FontWeight','bold','color',[199, 123, 179]./255)
    text(11,0.45,'Typical','FontSize',18,'FontWeight','bold','color',[88, 132, 209]./255)
    text(11,0.4,'Hearing','FontSize',18,'FontWeight','bold','color',[224, 191, 45]./255)
    title('mFus-faces shift','fontsize',22)
end

if hand
    exportgraphics(gcf, 'figure_OTS_shift_refLabelBootstrapped.pdf', 'ContentType', 'vector');
else
    exportgraphics(gcf, 'figure_mFus_shift_refLabelBootstrapped.pdf', 'ContentType', 'vector');
end
