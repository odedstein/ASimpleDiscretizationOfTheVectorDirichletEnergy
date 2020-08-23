addpath ..;
path_handling;


% Load mesh
[V,F] = readOBJ('venus_hr.obj');
[V,Io] = remove_unreferenced(V,F);
F = Io(F);

% Pick input VF
mn = mean(V);
rotAtVec = [0 1 0];
z0 = @(x) ([-(x(:,2)-mn(2)), (x(:,1)-mn(1)), zeros(size(x,1),1)]);
rng(21);
vertNoise = 1 + 0.5 * (rand(size(V,1),1)-0.5);
vertNoiseRot = 2*pi/3 * (rand(size(V,1),1) - 0.5);

% CR
t = 6e-3;
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

% Construct exact solution
exactVecs = z0(edgeMps);
notChangedMuchVecs = (1:m)';
for j=1:2
    exactVecs(:,3) = exactVecs(:,3) * 0.3;
    exactVecs = exactVecs - dot(exactVecs,edgeN,2).*edgeN;
    divFac = 1 + normrow(exactVecs);
    exactVecs = exactVecs ./ divFac;
end
exactVecs(isnan(exactVecs)) = 0;

% Add noise
vecs = exactVecs;
for j=1:size(vecs,1)
    vecs(j,:) = vecs(j,:) * ...
        axisangle_to_rotmatrix(edgeN(j,:),noiseRot(j))' * ...
        noise(j);
end

% find para and perp values at bdry
[para, perp] = project_onto_cr_space(vecs, (1:m)', V, F, E, oE);
% Solve actual system
crZ = [para;perp];
Q = (crM + t*crL);
for iter=1:100
    crZ = Q \ (crM * crZ);
end
% Construct vecs to plot
crNoisyVec = cr_plot_vec(V,F,E,oE,[para;perp]);
crPlotVec = cr_plot_vec(V,F,E,oE,crZ);


names = {'noisy', 'smoothed'};
vecs = {crNoisyVec, crPlotVec};

for i=1:numel(names)
    % Plot
    clf;
    col = [0.8 0.8 0.8];
    t = tsurf(F,V, 'CData',ones(size(V,1),1));
    axis equal;
    axis off;
    t.EdgeColor = 'none';
    set(t,fphong,'FaceVertexCData',repmat(col,size(V,1),1));
    shad = struct();
    shad.SpecularStrength = 0.1;
    shad.DiffuseStrength = 0.5;
    shad.AmbientStrength = 0.6;
    set(t,shad);
    camproj('persp');
    view(-137.53, 6.2523);
    l = light('Position',[-0.5 0.5 2], 'Style', 'Infinite');
    s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
    delete(l)
    add_lights_o();
    apply_ambient_occlusion(t,'Factor',1.3,'AddLights',false,'SoftLighting',false);
    caxis([0 1]);
    
    hold on;
    rng(8291237);
    [snapInds] = sample_and_snap(edgeMps,V,F,0.2);
    %offsetN = edgeN;
    offsetN = repmat((campos-camtarget) / norm(campos-camtarget), m, 1);
    edgePlotPts = edgeMps + 1*offsetN;
    qu = vec_plot(edgePlotPts(snapInds,:), vecs{i}(snapInds,:), ...
        'Scale',2, 'LogScale',false, ...
        'ArrowRadius',0.1, 'Color',[0.015686 0.35294 0.55294]);
    
    
    %set(gcf, 'Position',  [0, 0, 3840, 2160]);
    drawnow;
    
    %imwrite(myaa('raw'), ['teaser-smooth-' names{i} '.png']);
    saveas(gcf, ['teaser-smooth-' names{i} '.png']);
end