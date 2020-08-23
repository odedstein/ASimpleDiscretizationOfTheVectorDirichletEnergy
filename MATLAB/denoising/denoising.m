addpath ..;
path_handling;


% Load mesh
[V,F] = readOBJ('bunny.obj');
orign = size(V,1);
R = axisangle_to_rotmatrix([1 0 0],pi/2);
V = V*R';
[V,F,S] = loop(V,F,2);

% Pick input VF
z0 = @(x) ((x(:,3)<=0.36723)*1 + ...
    (x(:,3)>0.36723).*cos(4*(x(:,3)-0.36723))) .* ...
    [0.2*x(:,2)+cos(2*(x(:,2)+0.1)), 0.5*(x(:,1)+0.3), 0.3*cos(pi*x(:,3))];
rng(0);
vertNoise = 1 + 0.45 * S * ...
    (rand(orign,1)-0.5) + 0.6*(rand(size(V,1),1)-0.5);
vertNoiseRot = 2*pi/3 * (rand(size(V,1),1) - 0.5);

% CR
t = 2e-5;
[crL, crM, E, oE] = cr_vector_laplacian_and_mass(V,F);
m = size(crL,1) / 2;
aE = all_edges(F);
posEdges = find(oE(:)>0);
eE = nan(m,2); eE(E(posEdges),:) = aE(posEdges,:);
noise = 0.5*(vertNoise(eE(:,1),:) + vertNoise(eE(:,2),:));
noiseRot = 0.5*(vertNoiseRot(eE(:,1),:) + vertNoiseRot(eE(:,2),:));
% Deal with pts and vecs
edgeMps = 0.5 * (V(eE(:,1),:) + V(eE(:,2),:));
edgeN = indexed_per_edge_normals(V,F,E);
vecs = z0(edgeMps);
for j=1:size(vecs,1)
    vecs(j,:) = vecs(j,:) * ...
        axisangle_to_rotmatrix(edgeN(j,:),noiseRot(j))' * ...
        noise(j);
end
% find para and perp values at bdry
[paraExact, perpExact] = project_onto_cr_space(z0(edgeMps), (1:m)', V, F, E, oE);
[para, perp] = project_onto_cr_space(vecs, (1:m)', V, F, E, oE);
% Solve actual system
crZ = [para;perp];
Q = (crM + t*crL);
for iter=1:15
    crZ = Q \ (crM * crZ);
end
% Construct vecs to plot
crExactVec = cr_plot_vec(V,F,E,oE,[paraExact;perpExact]);
crNoisyVec = cr_plot_vec(V,F,E,oE,[para;perp]);
crPlotVec = cr_plot_vec(V,F,E,oE,crZ);



names = {'exact', 'noisy', 'smoothed'};
vecs = {crExactVec, crNoisyVec, crPlotVec};

for i=1:numel(names)
    % Plot
    fprintf('Plotting %s\n', names{i});
    clf;
    crMagVert = average_per_edge_onto_vertices(eE,normrow(vecs{i}));
    t = tsurf(F,V, 'CData',crMagVert);
    cm = cbrewer('YlOrBr', 500);
    colormap(cm);
    colorbar;
    axis equal;
    axis off;
    t.EdgeColor = 'none';
    
    set(t,fphong);
    shad = struct();
    shad.SpecularStrength = 0.2;
    shad.DiffuseStrength = 0.4;
    shad.AmbientStrength = 0.6;
    set(t,shad);
    camproj('persp');
    view(-7.9885, 24.512);
    l = light('Position',[-0.5 -0.5 2], 'Style', 'Infinite');
    s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
    delete(l)
    add_lights_o();
    apply_ambient_occlusion(t,'Factor',1,'AddLights',false,'SoftLighting',false);
    
    hold on;
    rng(8910761);
    [snapInds] = sample_and_snap(edgeMps,V,F,150);
    %offsetN = edgeN;
    offsetN = repmat((campos-camtarget) / norm(campos-camtarget), ...
        size(edgeMps,1), 1);
    offsetPlotPts = edgeMps + 0.1*offsetN;
    qu = vec_plot(offsetPlotPts(snapInds,:), vecs{i}(snapInds,:), ...
        'Scale',0.08, 'LogScale',false, ...
        'ArrowRadius',0.1, 'Color',[0.015686 0.35294 0.55294]); %[0.0078431 0.21961 0.3451]);
    caxis([min(normrow(crExactVec)), max(normrow(crExactVec))]);
    set(gcf, 'Position',  [0, 0, 3840, 2160]);
    drawnow;
    
    %imwrite(myaa('raw'), ['denoising-' names{i} '.png']);
    saveas(gcf, ['denoising-' names{i} '.png']);
end
