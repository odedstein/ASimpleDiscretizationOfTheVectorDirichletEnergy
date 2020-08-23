addpath ..;
path_handling;

% These are the numbers of boundaries we will have
%bdns = floor(2.^linspace(1,6,7)); %For figures in paper
bdns = floor(2.^linspace(3,6.5,8));  %For rebuttal performance evaln.

% These are the bdry conditions and the exact solution of the problem
bdrycond = @(x) 2*(x(:,3) - 0.5) .* ...
    [[-x(:,2), x(:,1)], zeros(size(x,1),1)];
% Solution 
soln = @(x) 2*(x(:,3) - 0.5) .* ...
    [[-x(:,2), x(:,1)], zeros(size(x,1),1)];


% Set up variables we care about
hs = nan(1,numel(bdns));
nverts = nan(1,numel(bdns));
l2errsCr = nan(1,numel(bdns));
soltimeCr = nan(1,numel(bdns));

% Do convergence experiment
for i=1:numel(bdns)
    fprintf('doing %d / %d\n', i, numel(bdns));
    
    % Construct cylinder
    [X,Y,Z] = cylinder(ones(1,2*bdns(i)), 6*bdns(i));
    [V,F] = surf_to_mesh(X,Y,Z);
    b = unique(outline(F));
    
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
    tic;
    crZ = min_quad_with_fixed(crL, zeros(size(crL,1),1), bd, bc);
    soltimeCr(i) = toc;
    % Construct vecs to plot
    crPlotVec = cr_plot_vec(V,F,E,oE,crZ);
    % Average vector magnitudes onto vertices
    crMag = average_per_edge_onto_vertices(eE, normrow(crPlotVec));
    
    % Measure L2 error
    M = crouzeix_raviart_massmatrix(V,F);
    edgeErr = normrow(crPlotVec - soln(edgeMps));
    l2errsCr(i) = sqrt(edgeErr'*M*edgeErr);
    
    % Save max edge length
    hs(i) = max(normrow(V(eE(:,1),:) - V(eE(:,2),:)));
    
    % Save verts
    nverts(i) = size(V,1);
    
end


% Plot actual image
clf;
solnVec = soln(V);
t = tsurf(F,V, 'CData',normrow(solnVec));
cm = cbrewer('BuGn',200);
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
view(0, 11.213);
l = light('Position',[-1 -1 2], 'Style', 'Infinite');
s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
delete(l)
add_lights_o();
apply_ambient_occlusion(t,'Factor',1,'AddLights',false,'SoftLighting',false);

hold on;
rng(946612456);
snapInds = sample_and_snap(V,V,F,80);
snapInds = snapInds(V(snapInds,2)<0);
offsetN = repmat((campos-camtarget) / norm(campos-camtarget), size(V,1), 1);
VPlotPts = V + 0.1*offsetN;
qu = vec_plot(VPlotPts(snapInds,:), solnVec(snapInds,:), ...
    'Scale',0.06, 'LogScale',false, ...
    'ArrowRadius',0.1, 'Color',[0 0 0]);

%set(gcf, 'Position',  [0, 0, 3840, 2160])
%imwrite(myaa('raw'), 'cylinder_fct.png');
saveas(gcf, 'cylinder_fct.png');

% Error plot
figure;
loglog(hs, l2errsCr, '-or', ...
    hs, hs, '--k', ...
    hs, hs.^2, '--k');
set(gca, 'xdir', 'reverse' )
xlabel('h');
ylabel('L2 error');
title('Bdry value problem with known solution, convergence.');
legend('CR', 'h', 'h^2');
saveas(gcf,'curvedexactconvergence-hs.eps','epsc');

clf;
loglog(nverts, l2errsCr, '-or', ...
    nverts, 1./sqrt(nverts), '--k', ...
    nverts, 1./nverts, '--k');
xlabel('nverts');
ylabel('L2 error');
title('Bdry value problem with known solution, convergence.');
legend('CR', '1/sqrt(n)', '1/n');
saveas(gcf, 'curvedexactconvergence-ns.eps', 'epsc');

clf;
loglog(soltimeCr, l2errsCr, '-or', ...
    soltimeCr, 1./sqrt(soltimeCr), '--k', ...
    soltimeCr, 1./soltimeCr, '--k');
xlabel('computation time t in secs');
ylabel('L2 error');
title('Bdry value problem with known solution, convergence.');
legend('CR', '1/sqrt(t)', '1/t');
saveas(gcf, 'curvedexactconvergence-comptime.eps', 'epsc');
