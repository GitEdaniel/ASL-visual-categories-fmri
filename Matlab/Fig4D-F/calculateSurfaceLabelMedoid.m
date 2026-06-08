function medoid_vertex = calculateSurfaceLabelMedoid(surface_path, label)
% Needs a surface and a "label" as either a full parth to a surface label
% that was defined on the surface you are inputting, or the list of
% vertices on that surface in which you want to derive the medoid. 
% NOTE: This code returns a vertex that is 1-indexed, because it will add 1
% to the label vertex indices when loading a label, or it assumes you've
% already added 1 to them if you're feeding it a list of vertices

surf_file  = surface_path;
%label_file =  label_path; 

% Check if label is a string
if ischar(label)
    label_file = label;  % Assign label_file if label is a string

    [label_vertices, ~] = read_label([], label_file);
    label_vertices = label_vertices(:,1) + 1;
end

% Check if it is a vector of indices, convert to column if necessary
if iscolumn(label) 
    label_vertices = label(:);  % Ensure label_vertices is a column vector
elseif isrow(label)
    label_vertices = label(:)';  
end


% Load surface
[vertices, faces] = read_surf(surf_file);
faces = faces + 1;   % convert from 0-based (FreeSurfer) to 1-based (MATLAB)


% Build full surface edge list
nV = size(vertices,1);

edges = [faces(:,[1 2]);
         faces(:,[2 3]);
         faces(:,[3 1])];

edges = unique(sort(edges,2),'rows');

% Get edge lengths
edge_lengths = sqrt(sum((vertices(edges(:,1),:) - ...
                         vertices(edges(:,2),:)).^2,2));

% Keep only edges that are fully inside label
label_mask = false(nV,1);
label_mask(label_vertices) = true;

keep = label_mask(edges(:,1)) & label_mask(edges(:,2));

edges_label = edges(keep,:);
edge_lengths_label = edge_lengths(keep);


% Build graph over label vertices and re-index label vertices to compact numbering
label_list = label_vertices(:);
nL = length(label_list);

map = zeros(nV,1);
map(label_list) = 1:nL;

edges_reindexed = [map(edges_label(:,1)), map(edges_label(:,2))];

G = graph(edges_reindexed(:,1), ...
          edges_reindexed(:,2), ...
          edge_lengths_label, ...
          nL);

% Compute geodesic medoid
total_dist = zeros(nL,1);

for k = 1:nL
    d = distances(G, k);
    total_dist(k) = sum(d);
end

[~, idx] = min(total_dist);

medoid_vertex = label_list(idx);
