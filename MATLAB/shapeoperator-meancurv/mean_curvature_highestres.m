addpath '..';
path_handling;

nrefs = 8;
mcCr = nan(1,nrefs);
mcCot = nan(1,nrefs);
hs = nan(1,nrefs);

for ref=1:nrefs
    fprintf('Doing ref %d\n', ref);
    
    [V, F] = readOBJ('fish.obj');
    [V, F] = loop(V, F, ref-1);
    
    [S, E, oE] = cr_shape_operator(V,F);
    
    m = max(E(:));
    aE = all_edges(F);
    posEdges = find(oE(:)>0);
    eE = nan(m,2); eE(E(posEdges),:) = aE(posEdges,:);
    edgeMps = 0.5 * (V(eE(:,1),:) + V(eE(:,2),:));
    hs(ref) = min(normrow(V(eE(:,1),:) - V(eE(:,2),:)));
    
    perFaceMeanCurv = 0.5 * abs(S(:,1) + S(:,5) + S(:,9));
    integratedMeanCurv = 0.5 * doublearea(V,F)' * perFaceMeanCurv;
    mcCr(ref) = integratedMeanCurv;
    
    L = cotmatrix(V,F);
    perFaceIntegratedMeanCurvCot = normrow(L*V);
    integratedMeanCurvCot = 0.5 * sum(perFaceIntegratedMeanCurvCot);
    mcCot(ref) = integratedMeanCurvCot;
    
    % Plot
    if (ref==nrefs-1)
        clf;
        plotf = perFaceMeanCurv;
        t = tsurf(F,V, 'FaceVertexCData', plotf);
        cm = cbrewer('PuBuGn', 500);
        colormap(cm);
        colorbar;
        axis equal;
        axis off;
        t.EdgeColor = 'none'; %t.EdgeColor = '#222222';
        
        sortedC = log(1+sort(plotf));
        caxis([sortedC(floor(0.05*end)), sortedC(floor(0.95*end))]);
        
        shad = struct();
        shad.SpecularStrength = 0.2;
        shad.DiffuseStrength = 0.4;
        shad.AmbientStrength = 0.6;
        set(t,shad);
        camproj('persp');
        view(138.39, 22.17);
        l = light('Position',[1 1 2], 'Style', 'Infinite');
        s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
        delete(l)
        add_lights_o();
        apply_ambient_occlusion(t,'Factor',1,'AddLights',false,'SoftLighting',false);
        
        title('log scale local curvature squared');
        
        set(gcf, 'Position',  [0, 0, 1280, 800]);
        %imwrite(myaa('raw'), 'mean_curvature_highestres.png');
        saveas(gcf, 'mean_curvature_highestres.png');
    end
    
end

errCr = abs(mcCr(1:(end-1)) - mcCr(end));
errCot = abs(mcCot(1:(end-1)) - mcCot(end));
hsm = hs(1:(end-1));

figure;
loglog(hsm, errCr, '-or', ...
    hsm, errCot, '-og', ...
    hsm, hsm, '--k', ...
    hsm, hsm.^2, '--k');
set(gca, 'Xdir', 'reverse');
xlabel('edge length h');
ylabel('Total mean curvature error');
title('Computing the total mean curvature of a fish');
legend('CR', 'Cot laplacian', 'h', 'h^2');

saveas(gcf, 'mean_curvature_highestres_asmyp.eps', 'epsc');
