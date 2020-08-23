addpath ..;
path_handling;


% Load mesh
[V,F] = readOBJ('mountain.obj');
smL = -cotmatrix(V,F);
smM = massmatrix(V,F);
V = (smM + 1e-2*smL) \ (smM*V);
[V,F] = loop(V, F, 1);
F = flip_ears(V, F);
orign = size(V,1);
R = axisangle_to_rotmatrix([1 0 0],pi/2);
V = V*R';

% CR
[crL, crM, E, oE] = cr_vector_laplacian_and_mass(V,F);
m = size(crL,1) / 2;
aE = all_edges(F);
posEdges = find(oE(:)>0);
eE = nan(m,2); eE(E(posEdges),:) = aE(posEdges,:);
edgeMps = 0.5 * (V(eE(:,1),:) + V(eE(:,2),:));
bdry = true(m,1);
bdry(E(oE<1)) = false;
interior = ~bdry;
bdry2 = [bdry; bdry];
interior2 = [interior; interior];

% Pick RHS, make sure avgs to zero
[~,maxpt] = max(V(:,3));
mp = V(maxpt,:);
analytic_f = @(x) [sin(0.05*normrow(x(:,1:2)-mp(:,1:2))), ...
    sin(0.063*normrow(x(:,1:2)-mp(:,1:2))), ...
    zeros(size(x,1),1)];
[para_f, perp_f] = project_onto_cr_space(analytic_f(edgeMps), ...
    (1:m)', V, F, E, oE);
f = [para_f; perp_f];
f = (crM + 0.5*crL) \ (crM*f);
f = f - sum(crM*f)/sum(diag(crM));
fm = crM * f;


% Solve Neumann problem
uneu = crL \ fm;

% Solve Dirichlet problem
usub = crL(interior2,interior2) \ fm(interior2);
udir(interior2) = usub;
udir(bdry2) = 0;


% Construct vecs to plot
crRhsVec = cr_plot_vec(V,F,E,oE,f);
crDirichletVec = cr_plot_vec(V,F,E,oE,udir);
crNeumannVec = cr_plot_vec(V,F,E,oE,uneu);

dncax = [min(min(udir),min(uneu)), max(max(udir),max(uneu))];
dirrat = (max(f)-min(f)) / (max(udir)-min(udir));
neurat = (max(f)-min(f)) / (max(uneu)-min(uneu));
frat = 1;



names = {'rhs', 'Dirichlet', 'Neumann'};
vecs = {crRhsVec, crDirichletVec, crNeumannVec};
rats = {frat, dirrat, neurat};

for iso=[true false]
    for i=1:numel(names)
        if(iso)
            ncols = 20;
        else
            ncols = 500;
        end
        
        % Plot
        fprintf('Plotting %s\n', names{i});
        clf;
        crMagVert = average_per_edge_onto_vertices(eE,normrow(vecs{i}));
        t = tsurf(F,V, 'CData',crMagVert);
        if(i==1)
            cm = cbrewer('Reds', ncols);
        elseif(i==2)
            cm = cbrewer('YlGn', 1.5*ncols);
            cm = cm(1:ncols, :);
        else
            cm = cbrewer('PuOr', 2*ncols);
            cm = cm(1:ncols, :);
        end
        colormap(cm);
        colorbar;
        axis equal;
        axis off;
        t.EdgeColor = 'none';
        
        if(iso)
            add_isolines(t,'LineWidth',2);
        end
        
        set(t,fphong);
        shad = struct();
        shad.SpecularStrength = 0.2;
        shad.DiffuseStrength = 0.4;
        shad.AmbientStrength = 0.6;
        set(t,shad);
        camproj('persp');
        view(-46.552, 30.891);
        l = light('Position',[-0.5 -0.5 2], 'Style', 'Infinite');
        s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
        delete(l)
        add_lights_o();
        apply_ambient_occlusion(t,'Factor',1,'AddLights',false,'SoftLighting',false);
        
        hold on;
        rng(3910721);
        [snapInds] = sample_and_snap(edgeMps,V,F,0.1);
        offsetN = repmat((campos-camtarget) / norm(campos-camtarget), ...
            size(edgeMps,1), 1);
        offsetPlotPts = edgeMps + 1.1*offsetN;
        qu = vec_plot(offsetPlotPts(snapInds,:), vecs{i}(snapInds,:), ...
            'Scale',rats{i}*2.7, 'LogScale',false, ...
            'ArrowRadius',0.1, 'Color',[0.015686 0.35294 0.55294]); %[0.0078431 0.21961 0.3451]);
        hold off;
        
        %set(gcf, 'Position',  [0, 0, 1920, 1080]);
        drawnow;
        
        namestr = sprintf('bdryconditions-%s-iso-%d.png', names{i}, iso);
        print(gcf, namestr,'-dpng','-r300');
    end
end
