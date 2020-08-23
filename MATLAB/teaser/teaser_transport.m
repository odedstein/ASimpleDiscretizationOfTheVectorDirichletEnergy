addpath ..;
path_handling;



% Load goat head
[V,F] = readOBJ('venus_hr.obj');
t = mean(edge_lengths(V,F), 'all').^2; %average edge length;

% Pick initial vertex and vector
pt = [-4.7504, 18.368, 41.73];
v = [0 0 -1];


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
% normalize v (this is ok, we are not comparing to anything)
edgeN = indexed_per_edge_normals(V,F,E);
v = v - dot(v,edgeN(crInd,:),2) .* edgeN(crInd,:);
v = v ./ normrow(v);
% find para and perp values
[para, perp] = project_onto_cr_space(v, crInd, V, F, E, oE);
% Perform the actual vector heat algorithm
Y0 = zeros(2*m,1);
Y0([crInd,crInd+m]) = [para,perp];
Yt = min_quad_with_fixed(0.5*(crM + t*crL), -crM*Y0, ...
    [crInd;crInd+m], [para;perp]);
u0 = zeros(m,1);
u0(crInd) = sqrt(para.^2+perp.^2);
ut = min_quad_with_fixed(0.5*(cM + t*cL), -cM*u0, crInd, u0(crInd));
phi0 = zeros(m,1);
phi0(crInd) = 1;
phit = min_quad_with_fixed(0.5*(cM + t*cL), -cM*phi0, crInd, 1);
Xt = repmat(ut ./ ...
    (phit.*sqrt(Yt(1:m,:).^2 + Yt((m+1):end,:).^2)), [2 1]) ...
    .* Yt;
% Construct vecs to plot
crVec = cr_plot_vec(V,F,E,oE,Xt);

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
rng(8917344);
[snapInds, presentInds] = sample_and_snap(edgeMps,V,F,0.17,edgeMps(crInd,:));
[az,el] = view;
viewmat = viewmtx(az,el,25);
offsetN = repmat((campos-camtarget) / norm(campos-camtarget), m, 1);
edgePlotPts = edgeMps + 1*offsetN;
qu = vec_plot(edgePlotPts(snapInds,:), crVec(snapInds,:), ...
    'Scale',2, 'LogScale',false, ...
    'ArrowRadius',0.1, 'Color',[0.015686 0.35294 0.55294], ...
    'Highlight',presentInds, 'HighlightColor',[0.87059 0.17647 0.14902]);


%set(gcf, 'Position',  [0, 0, 3840, 2160])
%imwrite(myaa('raw'), 'teaser_transport.png');
saveas(gcf, 'teaser_transport.png');