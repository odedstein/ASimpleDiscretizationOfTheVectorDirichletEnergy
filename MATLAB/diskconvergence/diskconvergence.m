addpath ..;
path_handling;

% These are the numbers of boundaries we will have
%bdns = 6 * 2.^(1:7); %For paper figure
bdns = 6 * 2.^linspace(4,7.5,8); %For performance evaluation

% These are the bdry conditions and the exact solution of the problem
bdrycond = @(x) [sin( 2*atan2(x(:,2),x(:,1)) ).^2, ...
    cos( 2*atan2(x(:,2),x(:,1)) ).^2, ...
    zeros(size(x,1),1)];
% Solution in mathematica script
soln = @(x) 0.5*[1-x(:,1).^4+6*(x(:,1).^2).*(x(:,2).^2) - x(:,2).^4, ...
    1+x(:,1).^4-6*(x(:,1).^2).*(x(:,2).^2) + x(:,2).^4 ...
    zeros(size(x,1),1)];


% Set up variables we care about
hs = nan(1,numel(bdns));
nverts = nan(1,numel(bdns));
l2errsCr = nan(1,numel(bdns));
l2errsKnop = nan(1,numel(bdns));
l2errsSharp = nan(1,numel(bdns));

% Do convergence experiment
for i=1:numel(bdns)
    fprintf('doing %d / %d\n', i, numel(bdns));
    
    % Construct disk
    [V,F,b] = disk(bdns(i));
    V = [V zeros(size(V,1),1)];
    
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
    % Average vector magnitudes onto vertices
    crMag = average_per_edge_onto_vertices(eE, normrow(crPlotVec));
    
    % Measure L2 error
    M = crouzeix_raviart_massmatrix(V,F);
    edgeErr = normrow(crPlotVec - soln(edgeMps));
    l2errsCr(i) = sqrt(edgeErr'*M*edgeErr);
    
    
    % Save max edge length & verts
    hs(i) = max(normrow(V(eE(:,1),:) - V(eE(:,2),:)));
    nverts(i) = size(V,1);
end


% Plot actual image
clf;
solnVec = soln(V);

tsurf(F,V, 'CData',normrow(solnVec) );
cm = parula(5000);
colormap(cm);
colorbar;
c = caxis;
shading interp;
axis equal;
axis off;

hold on;
rng(946616456);
snapInds = sample_and_snap(V,V,F,40);
offsetN = repmat((campos-camtarget) / norm(campos-camtarget), size(V,1), 1);
VPlotPts = V + 0.1*offsetN;
qu = vec_plot(VPlotPts(snapInds,:), solnVec(snapInds,:), ...
    'Scale',0.15, 'LogScale',false, ...
    'ArrowRadius',0.1, 'Color',[0 0 0]);
caxis(c);

%set(gcf, 'Position',  [0, 0, 3840, 2160])
%imwrite(myaa('raw'), 'disk_fct.png');
saveas(gcf, 'disk_fct.png');


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
saveas(gcf, 'diskconvergence-hs.eps', 'epsc');

clf;
loglog(nverts, l2errsCr, '-or', ...
    nverts, 1./sqrt(nverts), '--k', ...
    nverts, 1./nverts, '--k');
xlabel('nverts');
ylabel('L2 error');
title('Bdry value problem with known solution, convergence.');
legend('CR', '1/sqrt(n)', '1/n');
saveas(gcf, 'diskconvergence-ns.eps', 'epsc');
