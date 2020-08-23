addpath ..;
path_handling;

% This is the regularization smoothing parameter
alpha = 1e-4;



% Load mesh
[V,F] = readOBJ('sodabottle_remeshed_hr.obj');
vertexBd = unique(outline(F));


% Do CR killing vectors
% construct matrices
[crL, E, oE] = cr_vector_laplacian(V,F);
[crtL] = cr_twisted_vector_laplacian(V,F,E,oE);
m = size(crL,1) / 2;
aE = all_edges(F);
posEdges = find(oE(:)>0);
eE = nan(m,2); eE(E(posEdges),:) = aE(posEdges,:);
edgeMps = 0.5 * (V(eE(:,1),:) + V(eE(:,2),:));
edgeN = indexed_per_edge_normals(V, F, E);

%boundary edges
bd = bdry_edge_indices(E);
bc = zeros(size(bd,1),1);
bc(1:size(bd,1)) = 1;
bd = [bd;m+bd];
bc((size(bd,1)/2+1):size(bd,1)) = 0;

% Construct killing matrix, solve bvp
killingEnergy = crL + crtL + alpha*crL;
z = min_quad_with_fixed(killingEnergy,zeros(size(killingEnergy,1),1), ...
    bd, bc);
% Plot vectors
crVec = cr_plot_vec(V,F,E,oE,z);

% Average vector magnitudes onto vertices
crMagVert = average_per_edge_onto_vertices(eE, normrow(crVec));
crScal = 1 / max(normrow(crVec));
crMagVert = crMagVert * crScal;
crVec = crVec * crScal;

% Do rotating around cylindrical axis field
rotateVecs = [edgeMps(:,1:2)*[0 1; -1 0], ...
    zeros(size(edgeMps,1),1)];
rotateMag = normrow(rotateVecs);
rotateMagVert = average_per_edge_onto_vertices(eE, normrow(rotateMag));
%rotateScal = 1 / mean(rotateMagVert(vertexBd));
rotateScal = 1 / max(rotateMag);
rotateMagVert = rotateMagVert * rotateScal;
rotateVecs = rotateVecs * rotateScal;

% Plot
verts = {V, V};
faces = {F, F};
names = {'CR', 'rotate'};
plotPts = {edgeMps, edgeMps};
vecs = {crVec, rotateVecs};
mags = {crMagVert, rotateMagVert};

for i=1:numel(plotPts)
    % Plot
    clf;
    t = tsurf(faces{i},verts{i}, 'CData',mags{i});
    shading interp;
    caxis([0.3, 1]);
    cm = flipud(cbrewer('RdGy', 500));
    colormap(cm);
    colorbar;
    axis equal;
    axis off;
    t.EdgeColor = 'none';
    set(t,fphong);
    shad = struct();
    shad.SpecularStrength = 0.3;
    shad.DiffuseStrength = 0.6;
    shad.AmbientStrength = 0.2;
    set(t,shad);
    camproj('persp');
    view(0, 20.764);
    l = light('Position',[0.5 -0.5 2], 'Style', 'Infinite');
    s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
    delete(l)
    add_lights_o();
    apply_ambient_occlusion(t,'Factor',2,'AddLights',false,'SoftLighting',false);
    
    hold on;
    rng(89107);
    [snapInds] = sample_and_snap(plotPts{i},V,F,0.02);
    offsetN = repmat((campos-camtarget) / norm(campos-camtarget), size(plotPts{i},1), 1);
    offsetPlotPts = plotPts{i} + 5*offsetN;
    qu = vec_plot(offsetPlotPts(snapInds,:), vecs{i}(snapInds,:), ...
        'Scale',5, 'LogScale',false, ...
        'ArrowRadius',0.1, 'Color',[0.015686 0.35294 0.55294]);
    caxis([0.3, 1]);
    
    %set(gcf, 'Position',  [0, 0, 3840, 2160]);
    %imwrite(myaa('raw'), ['killingfield-' names{i} '.png']);
    saveas(gcf, ['killingfield-' names{i} '.png']);
end

