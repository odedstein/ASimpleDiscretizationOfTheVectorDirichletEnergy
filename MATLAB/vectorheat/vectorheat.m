addpath ..;
path_handling;

% Load model
[V,F] = readOBJ('ramskull.obj');
tcr = mean(edge_lengths(V,F), 'all'); %average edge length;

% Pick initial vertex and vector
pt = [0.029439 0.554511 0.330999];
v = [0.030781 0.614479 0.361481] - pt; v = v/norm(v);

% Do CR vector heat
% construct matrices
[crL, crM, E, oE] = cr_vector_laplacian_and_mass(V,F);
m = size(crL,1) / 2;
aE = all_edges(F);
posEdges = find(oE(:)>0);
eE = nan(m,2); eE(E(posEdges),:) = aE(posEdges,:);
cL = crouzeix_raviart_cotmatrix(V,F);
cM = crouzeix_raviart_massmatrix(V,F);
% find closest point to pt in edges by index
edgeMps = 0.5 * (V(eE(:,1),:) + V(eE(:,2),:));
crInd = pts_closest_to(edgeMps, pt);
% find para and perp values
[para, perp] = project_onto_cr_space(v, crInd, V, F, E, oE);
% Perform the actual vector heat algorithm
Y0 = zeros(2*m,1);
Y0([crInd,crInd+m]) = [para,perp];
Yt = min_quad_with_fixed(0.5*(crM + tcr*crL), -crM*Y0, ...
    [crInd;crInd+m], [para;perp]);
u0 = zeros(m,1);
u0(crInd) = sqrt(para.^2+perp.^2);
ut = min_quad_with_fixed(0.5*(cM + tcr*cL), -cM*u0, crInd, u0(crInd));
phi0 = zeros(m,1);
phi0(crInd) = 1;
phit = min_quad_with_fixed(0.5*(cM + tcr*cL), -cM*phi0, crInd, 1);
Xt = repmat(ut ./ ...
    (phit.*sqrt(Yt(1:m,:).^2 + Yt((m+1):end,:).^2)), [2 1]) ...
    .* Yt;
% Construct vecs to plot
crVec = cr_plot_vec(V,F,E,oE,Xt);

% Do scalar vector heat for color
L = -cotmatrix(V,F);
M = massmatrix(V,F,'barycentric');
indInV = pts_closest_to(V, pt);
phi0 = zeros(size(V,1),1);
phi0(indInV) = 1;
phit = min_quad_with_fixed(0.5*(M + tcr*L), -M*phi0, indInV, 1);
utscal = phit;
gradut = grad(V,F)*utscal;
scalarX = -gradut ./ normrow(gradut);
scalarX(isnan(scalarX)) = 0;
scalarHeat = min_quad_with_fixed(-0.5*L, div(V,F)*scalarX, indInV, 1);


names = {'raw', 'CR'};
plotPts = {pt, edgeMps};
vecs = {v, crVec};
fixeds = {1, crInd};

for i=1:numel(plotPts)
    % Plot
    clf;
    t = tsurf(F,V, 'CData',scalarHeat);
    shading interp;
    cm = hot(300);
    cm = cm(1:floor(5*end/6),:).^4;
    colormap(cm);
    colorbar;
    axis equal;
    axis off;
    t.EdgeColor = 'none';
    set(t,fphong);
    shad = struct();
    shad.SpecularStrength = 0.4;
    shad.DiffuseStrength = 0.8;
    shad.AmbientStrength = 0.3;
    set(t,shad);
    c = caxis();
    camproj('persp');
    view(58.535, 35.394);
    camroll(25);
    l = light('Position',[0.5 0.5 2], 'Style', 'Infinite');
    s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
    delete(l)
    add_lights_o();
    apply_ambient_occlusion(t,'Factor',1,'AddLights',false,'SoftLighting',false);
    
    hold on;
    rng(129758);
    [snapInds, presentInds] = sample_and_snap(plotPts{i},V,F,40,plotPts{i}(fixeds{i},:));
    offsetN = repmat((campos-camtarget) / norm(campos-camtarget), size(plotPts{i},1), 1);
    offsetPlotPts = plotPts{i} + 0.03*offsetN;
    qu = vec_plot(offsetPlotPts(snapInds,:), vecs{i}(snapInds,:), ...
        'Scale',0.1, 'LogScale',false, ...
        'ArrowRadius',0.1, 'Color',[0.015686 0.35294 0.55294], ...
        'Highlight',presentInds, ...
        'HighlightColor',[0.87059 0.17647 0.14902]);
    caxis(c);
    
    %set(gcf, 'Position',  [0, 0, 3840, 2160]);
    %imwrite(myaa('raw'), ['vectorheat-' names{i} '.png']);
    saveas(gcf,  ['vectorheat-' names{i} '.png']);

end



