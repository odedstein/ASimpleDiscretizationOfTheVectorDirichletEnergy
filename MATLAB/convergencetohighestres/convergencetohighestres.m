addpath ..;
path_handling;

% These are the numbers of times we will subdivise
ns = 6;

% These are the bdry conditions and the exact solution of the problem
bdrycond = @(x) cos(atan2(x(:,2),x(:,1))).^2 .* ...
    [[-x(:,2), x(:,1)], zeros(size(x,1),1)];

% These are the meshes we care about
meshes = {'spot.obj', 'spot-unequalresolns.obj'};

for meshid=1:numel(meshes)
    % Load mesh
    meshname = meshes{meshid};
    [V,F] = readOBJ(meshname);
    
    
    % Plot raw mesh
    clf;
    t = tsurf(F,V, 'CData',ones(size(V,1),1));
    
    t.EdgeColor = '#111111';
    t.FaceColor = '#D0D0D0';
    camproj('persp');
    axis equal;
    axis off;
    l = light('Position',[-1 -1 2], 'Style', 'Infinite');
    s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
    apply_ambient_occlusion(t,'Factor',1,'AddLights',false,'SoftLighting',true);
    view(-21.792,24.549);
    delete(l)
    
    set(gcf, 'Position',  [0, 0, 1440, 900]);
    %imwrite(myaa('raw'), ['spot_mesh_' meshname '.png']);
    saveas(gcf, ['spot_mesh_' meshname '.png']);
    
    
    % Set up variables we care about
    hs = nan(1,ns-1);
    perVertCr = cell(1,ns);
    S = cell(1,ns-1); %Subdivision mat to highest-res
    J = cell(1,ns-1);
    
    % Do convergence experiment
    for i=1:ns
        fprintf('doing %d / %d\n', i, ns);
        
        % Construct spot mesh, fix bdry, update subdivision matrix
        [V,F,Ss,Jj] = loop(V,F);
        for j=1:(i-1)
            S{j} = Ss * S{j};
            oldJ = J{j};
            J{j} = oldJ(Jj);
        end
        b = unique(outline(F));
        V(b,3) = 0;
        V(b,1:2) = V(b,1:2) ./ normrow(V(b,1:2));
        
        if(i~=ns)
            S{i} = speye(size(V,1), size(V,1));
            J{i} = (1:size(F,1))';
        end
        
        
        % CR
        [crL, E, oE] = cr_vector_laplacian(V,F);
        m = size(crL,1) / 2;
        aE = all_edges(F);
        posEdges = find(oE(:)>0);
        eE = nan(m,2); eE(E(posEdges),:) = aE(posEdges,:);
        % Deal with boundary conditions
        edgeMps = 0.5 * (V(eE(:,1),:) + V(eE(:,2),:));
        bd = bdry_edge_indices(E);
        % find para and perp values
        [para, perp] = project_onto_cr_space(bdrycond(edgeMps(bd,:)), ...
            bd, V, F, E, oE);
        bc = [para; perp];
        bd = [bd; bd+m];
        % Solve actual system
        crZ = min_quad_with_fixed(crL, zeros(size(crL,1),1), bd, bc);
        % Construct vecs to plot
        crPlotVec = cr_plot_vec(V,F,E,oE,crZ);
        % Average plot vecs onto vertices
        perVertCr{i} = [average_per_edge_onto_vertices(eE, crPlotVec(:,1)), ...
            average_per_edge_onto_vertices(eE, crPlotVec(:,2)), ...
            average_per_edge_onto_vertices(eE, crPlotVec(:,3))];
        
        % Save max edge length
        if(i~=ns)
            hs(i) = max(normrow(V(eE(:,1),:) - V(eE(:,2),:)));
        end
    end
    
    
    % Compute actual errors vs highest-res resolution using S
    l2errsCr = nan(1,ns-1);
    M = massmatrix(V,F,'full');
    for i=1:(ns-1)
        %Error calc by backporting single vector from hr solution onto old
        % solution, then integrating on HR mesh the LR per-vert fcts
        Mloc = S{i}' * M * S{i};
        endCr = perVertCr{end};
        crErr = normrow(perVertCr{i} - endCr(1:size(perVertCr{i},1),:));
        l2errsCr(i) = sqrt(crErr' * Mloc * crErr);
    end
    
    
    
    crEnd = perVertCr{end};
    crEndPlot = log(1+1000*normrow(crEnd)) .* crEnd ./ normrow(crEnd);
    crEndPlot(isnan(crEndPlot)) = 0;
    
    
    
    % Plot
    clf;
    t = tsurf(F,V, 'CData',normrow(crEndPlot));
    cm = cbrewer('PuBu',200);
    colormap(cm);
    colorbar;
    shading interp;
    axis equal;
    axis off;
    t.EdgeColor = 'none';
    shad = struct();
    shad.SpecularStrength = 0.1;
    shad.DiffuseStrength = 0.5;
    shad.AmbientStrength = 0.6;
    set(t,shad);
    camproj('persp');
    view(-50.691, 17.036);
    l = light('Position',[-1 -1 2], 'Style', 'Infinite');
    s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
    delete(l)
    add_lights_o();
    apply_ambient_occlusion(t,'Factor',1,'AddLights',false,'SoftLighting',false);
    
    hold on;
    rng(865752);
    [snapInds] = sample_and_snap(V,V,F,5);
    offsetN = repmat((campos-camtarget) / norm(campos-camtarget), size(V,1), 1);
    VPlotPts = V + 0.2*offsetN;
    qu = vec_plot(VPlotPts(snapInds,:), crEndPlot(snapInds,:), ...
        'Scale',0.03, 'LogScale',false, ...
        'ArrowRadius',0.1, 'Color',[0 0 0]);
    
    %set(gcf, 'Position',  [0, 0, 3840, 2160])
    %imwrite(myaa('raw'), ['spot_fct_' meshname '.png']);
    saveas(gcf, ['spot_fct_' meshname '.png']);
    
    % Error plot
    figure;
    loglog(hs, l2errsCr, '-or', ...
        hs, hs, '--k', ...
        hs, hs.^2, '--k');
    set(gca, 'xdir', 'reverse' )
    xlabel('h');
    ylabel('L2 error');
    title(['Bdry value problem, convergence to highest-res solution of each algo.', ...
        'vecs are log scaled']);
    legend('CR', 'h', 'h^2');
    saveas(gcf,['convergencetohighestres-' meshname '.eps'],'epsc');
end
