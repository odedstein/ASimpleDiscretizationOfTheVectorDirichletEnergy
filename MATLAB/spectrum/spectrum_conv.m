addpath ..;
path_handling;

refs = 7;

errs_cr = nan(1,refs);
hs = nan(1,refs);

for meshind = 1:refs
    fprintf('Computing ref %d\n', meshind);
    
    % Load mesh
    [V,F] = readOBJ('sphere-lo.obj');
    [V,F] = loop(V,F,meshind-1);
    V = V ./ normrow(V);
    
    % CR
    [crL, crM, E, oE] = cr_vector_laplacian_and_mass(V,F);
    m = size(crL,1) / 2;
    aE = all_edges(F);
    posEdges = find(oE(:)>0);
    eE = nan(m,2); eE(E(posEdges),:) = aE(posEdges,:);
    [~, eDs] = eigs(crL, crM, 3, 'smallestabs');
    errs_cr(meshind) = abs(1-eDs(1,1));
    hs(meshind) = max(normrow(V(eE(:,1),:) - V(eE(:,2),:)));
    
    if(meshind ~= 4)
        continue;
    end
    
    % Plot
    clf;
    t = tsurf(F,V, 'CData',ones(size(V,1),1));
    cm = cbrewer('PuBu', 500);
    colormap(cm);
    axis equal;
    axis off;
    t.EdgeColor = 'none';
    
    set(t,fphong);
    shad = struct();
    shad.SpecularStrength = 0.2;
    shad.DiffuseStrength = 0.4;
    shad.AmbientStrength = 0.6;
    set(t,shad);
    camproj('persp');
    view(-7.9885, 24.512);
    l = light('Position',[-0.5 -0.5 2], 'Style', 'Infinite');
    s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
    delete(l)
    add_lights_o();
    apply_ambient_occlusion(t,'Factor',1,'AddLights',false,'SoftLighting',false);
    set(gcf, 'Position',  [0, 0, 1440, 900]);
    
    %imwrite(myaa('raw'), 'spectrum-surf-conv.png');
    saveas(gcf, 'spectrum-surf-conv.png');
end

clf;
loglog(hs, errs_cr, '-or', ...
    hs, hs, '--k', ...
    hs, hs.^2, '--k');
set(gca, 'xdir', 'reverse');
xlabel('h');
ylabel('eigenvalue error');
title('eigenvalue of sphere, convergence.');
legend('CR', 'h', 'h^2');
saveas(gcf, 'spectrum-surf-conv-plot.eps', 'epsc');
    
