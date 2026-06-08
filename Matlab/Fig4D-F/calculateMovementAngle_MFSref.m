function [angle_D_MFS, angle_C1C2_MFS] = calculateMovementAngle_MFSref(surface_path, originVertex, terminusVertex1, terminusVertex2, refVector_vert1, refVector_vert2)
%
% This code was written to better understand how a given functional ROI is
% moving in deaf participants compared to control. This was written
% specifically to synergize with the code deafMovementBootstrap. It takes
% in 3 vertices. The first is the origin vertex from which two vectors
% will be derived, for this origin vertex we will use the medoid of the
% control group bootstrap C1. The next is the medoid of the deaf group D,
% and the next is the medoid of the second control group bootstrap C2. We
% are interested in calculating the angle formed by vector C1-D and C1-C2.
% We will also rotate the vectors so that they align with a reference
% vector formed by the endpoints of the MFS, so refVector_vert1 is
% posterior MFS and refVector_vert2 is anterior MFS. We do this so that we
% can visualize the distribution of movement of C1-D and show that it is
% distinct from the spatial variation of bootstrapped vectors C1-C2. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C1 = originVertex;
D = terminusVertex1;
C2 = terminusVertex2;
MFS1 = refVector_vert1;
MFS2 = refVector_vert2;

% Load the surface to get the vertex coordinates
[vertices, surfFaces] = read_surf(surface_path);

D = vertices(D,:);
C1 = vertices(C1,:);
C2 = vertices(C2,:);
MFS1 = vertices(MFS1,:);
MFS2 = vertices(MFS2,:);

% Define vectors
vector_C1D  = D   - C1;
vector_C1C2 = C2  - C1;
% and reference vector
vector_MFS = MFS2 - MFS1;

% Compute angle using dot product
angle_DC = acos( dot(vector_C1D, vector_C1C2) / ...
                (norm(vector_C1D) * norm(vector_C1C2)) );
angle_DC = rad2deg(angle_DC);

% angle between DC and MFS
angle_D_MFS = acos( dot(vector_C1D, vector_MFS) / ...
                (norm(vector_C1D) * norm(vector_MFS)) );
angle_D_MFS = rad2deg(angle_D_MFS);

% angle between C1C2 and MFS
angle_C1C2_MFS = acos( dot(vector_C1C2, vector_MFS) / ...
                (norm(vector_C1C2) * norm(vector_MFS)) );
angle_C1C2_MFS = rad2deg(angle_C1C2_MFS);

% Below if is if you'd like to rotate the vectors so that C1D becomes
% aligned with MFS.

% % Normalize both vectors
% v1 = vector_C1D  / norm(vector_C1D);
% v2 = vector_MFS  / norm(vector_MFS);
% 
% % compute rotation axis
% rot_axis = cross(v1, v2);
% rot_axis = rot_axis / norm(rot_axis);
% 
%  % get rotation angle
% rot_angle = acos(dot(v1, v2));
% 
% % Build rotation matrix
% K = [     0        -rot_axis(3)  rot_axis(2);
%       rot_axis(3)       0       -rot_axis(1);
%      -rot_axis(2)  rot_axis(1)       0     ];
% 
% R = eye(3) + sin(rot_angle)*K + ...
%     (1 - cos(rot_angle))*(K*K);
% 
% % apply rotation
% vector_C1D_rot  = (R * vector_C1D')';
% vector_C1C2_rot = (R * vector_C1C2')';
% 
% % compute angle between rotate vectors
% angle_DC_rot = acos( dot(vector_C1D_rot, vector_C1C2_rot) / ...
%                     (norm(vector_C1D_rot) * norm(vector_C1C2_rot)) );
% 
% angle_DC_rot = rad2deg(angle_DC_rot);
% 
% vec1_rot = vector_C1D_rot;
% vec2_rot = vector_C1C2_rot;
% angle = angle_DC_rot;




