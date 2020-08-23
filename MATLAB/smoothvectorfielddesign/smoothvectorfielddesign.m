addpath ..;
path_handling;


% Load mesh
[V,F] = readOBJ('koala.obj');

% Pick initial vertices and vectors
pts = [-3.48029, -1.46583, 3.84537; ...
    -3.48029, 1.46583, 3.84537; ...
    2.38831, 0, 3.42307];
vecs = [-1, 0, 0; ...
    -1, 0, 0; ...
    -1, 0, 0];

% CR
[crL, E, oE] = cr_vector_laplacian(V,F);
m = size(crL,1) / 2;
aE = all_edges(F);
posEdges = find(oE(:)>0);
eE = nan(m,2); eE(E(posEdges),:) = aE(posEdges,:);
% Deal with pts and vecs
edgeMps = 0.5 * (V(eE(:,1),:) + V(eE(:,2),:));
crInd = pts_closest_to(edgeMps, pts);
% find para and perp values at bdry
[para, perp] = project_onto_cr_space(vecs, crInd, V, F, E, oE);
% Solve actual system
crZ = min_quad_with_fixed(crL, zeros(size(crL,1),1), ...
    [crInd; crInd+m], [para;perp]);
% Construct vecs to plot
crPlotVec = cr_plot_vec(V,F,E,oE,crZ);



names = {'none', 'CR'};
plotPts = {edgeMps(crInd,:), edgeMps};
fixed = {(1:size(crInd,1))', crInd};
vecs = {crPlotVec(crInd,:), crPlotVec};

for i=1:numel(plotPts)
    % Resize vecs to a log scale
    if(strcmp('none', names{i}))
        szby = (log(0.1+normrow(vecs{2})) - min(log(0.1+normrow(vecs{2}))));
        resizedVec = vecs{i} ./ normrow(vecs{i}) .* szby(crInd);
    else
        resizedVec = vecs{i} ./ normrow(vecs{i}) .* ...
            (log(0.1+normrow(vecs{i})) - min(log(0.1+normrow(vecs{i}))));
    end
    
    % Plot
    clf;
    if(strcmp('none', names{i}))
        t = tsurf(F,V, 'CData', 0.1*ones(size(V,1),1));
    else
        crMagVert = average_per_edge_onto_vertices(eE,normrow(resizedVec));
        t = tsurf(F,V, 'CData',crMagVert);
    end
    cm = cbrewer('YlGn', 500);
    colormap(cm);
    colorbar;
    caxis([min(normrow(crPlotVec)), max(normrow(crPlotVec))]);
    axis equal;
    axis off;
    t.EdgeColor = 'none';
    set(t,fphong);
    shad = struct();
    shad.SpecularStrength = 0.2;
    shad.DiffuseStrength = 0.3;
    shad.AmbientStrength = 0.7;
    set(t,shad);
    camproj('persp');
    view(-42.438, 20.764);
    l = light('Position',[-0.5 -0.5 2], 'Style', 'Infinite');
    s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
    delete(l)
    add_lights_o();
    apply_ambient_occlusion(t,'Factor',1,'AddLights',false,'SoftLighting',false);
    
    if(strcmp('none', names{i}))
        snapInds = fixed{i};
        presentInds = fixed{i};
    else
        rng(891072);
        [snapInds, presentInds] = sample_and_snap(plotPts{i},V,F,10,plotPts{i}(fixed{i},:));
    end
    
    
    offsetN = repmat((campos-camtarget) / norm(campos-camtarget), size(plotPts{i},1), 1);
    offsetPlotPts = plotPts{i} + 0.2*offsetN;
    offsetPlotPts(fixed{i},:) = offsetPlotPts(fixed{i},:) + 0.2*offsetN(fixed{i},:);
    hold on;
    qu = vec_plot(offsetPlotPts(snapInds,:), resizedVec(snapInds,:), ...
        'Scale',0.5, 'LogScale',false, ...
        'ArrowRadius',0.1, 'Color',[0.015686 0.35294 0.55294], ...
        'Highlight',presentInds, 'HighlightColor',[0.87059 0.17647 0.14902]);
    
    %set(gcf,'Pos',[1442 0 1680 1050]);
    saveas(gcf, ['smoothvectorfielddesign-' names{i} '.png']);
end
