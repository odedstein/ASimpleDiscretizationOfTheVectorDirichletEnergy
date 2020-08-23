addpath '..';
path_handling;

meshes = {'mushroom.obj', 'brucewick.obj', 'torusex.obj'};

for meshnum = 1:numel(meshes)
    meshname = meshes{meshnum};
    
    [V,F] = readOBJ(meshname);
    switch(meshnum)
        case 1
            V = (massmatrix(V,F) - 2e-2*cotmatrix(V,F)) \ (massmatrix(V,F)*V);
            [V,F] = loop(V, F, 2);
        case 2
            [V,F] = loop(V, F, 3);
        case 3
            [V,F] = loop(V, F, 1);
    end
    m = size(F,1);
    
    [S, E, oE] = cr_shape_operator(V,F);
    k1s = nan(m,1);
    k2s = nan(m,1);
    kns = nan(m,1);
    u1s = nan(m,3);
    u2s = nan(m,3);
    uns = nan(m,3);
    for f=1:m
        Sl = reshape(S(f,:), 3, 3);
        [evecs, evals] = eig(Sl);
        [~,I] = sort(abs(diag(evals)));
        kns(f) = evals(I(1),I(1));
        k1s(f) = evals(I(2),I(2));
        k2s(f) = evals(I(3),I(3));
        uns(f,:) = evecs(:,I(1));
        u1s(f,:) = evecs(:,I(2));
        u2s(f,:) = evecs(:,I(3));
    end
    
    %Principal curvatures are unreliable at umbilic points and at
    % points with zero Gaussian curvature
    unreliable = abs(k1s-k2s)<1e-6 | abs(k1s)<1e-6 | abs(k2s)<1e-6;
    u1s(unreliable,:) = 0;
    u2s(unreliable,:) = 0;
    
    N = normals(V,F);
    N = N ./ normrow(N);
    
    gauss = k1s .* k2s;
    
    % Plot
    clf;
    plotf = sign(gauss).*log(1+abs(gauss));
    t = tsurf(F,V, 'CData',plotf);
    cm = flipud(cbrewer('RdBu',500));
    colormap(cm);
    colorbar;
    shading interp;
    axis equal;
    axis off;
    sortedg = sort(plotf);
    mmxu = max(abs(sortedg(floor(0.05*end))), abs(sortedg(floor(0.95*end))));
    caxis([-mmxu, mmxu]);
    t.EdgeColor = 'none';
    shad = struct();
    shad.SpecularStrength = 0.2;
    shad.DiffuseStrength = 0.6;
    shad.AmbientStrength = 0.5;
    set(t,shad);
    camproj('persp');
    view(-33.641, 16.329);
    l = light('Position',[-1 -1 2], 'Style', 'Infinite');
    s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
    delete(l)
    add_lights_o();
    apply_ambient_occlusion(t,'Factor',1,'AddLights',false,'SoftLighting',false);
    
    hold on;
    switch(meshnum)
        case 1
            scale = 0.09;
            sampleres = 10;
        case 2
            scale = 0.045;
            sampleres = 30;
        case 3
            scale = 0.03;
            sampleres = 100;
    end
    rng(865532);
    B = barycenter(V,F);
    [snapInds] = sample_and_snap(B,V,F,sampleres);
    offsetN = repmat((campos-camtarget) / norm(campos-camtarget), size(B,1), 1);
    VPlotPts = B + 0.05*offsetN;
    col1 = [0.3 0.3 0.36];
    col2 = [0.7 0.7 0.8];
    qu1 = vec_plot(VPlotPts(snapInds,:), u1s(snapInds,:), ...
        'Scale',scale, 'LogScale',false, ...
        'ArrowRadius',0.3, 'Color',col1, 'Head',false);
    qu15 = vec_plot(VPlotPts(snapInds,:), -u1s(snapInds,:), ...
        'Scale',scale, 'LogScale',false, ...
        'ArrowRadius',0.3, 'Color',col1, 'Head',false);
    qu2 = vec_plot(VPlotPts(snapInds,:), u2s(snapInds,:), ...
        'Scale',scale, 'LogScale',false, ...
        'ArrowRadius',0.3, 'Color',col2, 'Head',false);
    qu25 = vec_plot(VPlotPts(snapInds,:), -u2s(snapInds,:), ...
        'Scale',scale, 'LogScale',false, ...
        'ArrowRadius',0.3, 'Color',col2, 'Head',false);
    
    title('Principal curvature directions. Color is Gauss curvature, log scale.');
    
    set(gcf, 'Position',  [0, 0, 1280, 800])
    name = [meshname '_principalcurvs.png'];
    %imwrite(myaa('raw'), name);
    saveas(gcf, name);
    
end
