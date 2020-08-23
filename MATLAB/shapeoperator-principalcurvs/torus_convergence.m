addpath '..';
path_handling;

resols = 4 * 2.^linspace(1,7,8);
R = 0.5;
r = 0.1;

integratedk1 = nan(1,numel(resols));
integratedk2 = nan(1,numel(resols));
integratedgauss = nan(1,numel(resols));
hs = nan(1,numel(resols));

for resol=1:numel(resols)
    fprintf('Doing resol %d\n', resol);
    
    uzh = resols(resol);
    [V,F] = torus(R, r, floor(uzh*3.4), floor(uzh), 0);
    m = size(F,1);
    
    [S, E, oE] = cr_shape_operator(V,F);
    k1s = nan(m,1);
    k2s = nan(m,1);
    u1s = nan(m,3);
    u2s = nan(m,3);
    for f=1:m
        Sl = reshape(S(f,:), 3, 3);
        [evecs, evals] = eig(Sl);
        [~,I] = sort(abs(diag(evals)));
        k1s(f) = evals(I(2),I(2));
        k2s(f) = evals(I(3),I(3));
        u1s(f,:) = evecs(:,I(2));
        u2s(f,:) = evecs(:,I(3));
    end
    gauss = k1s .* k2s;
    
    [k1ex, k2ex, u1ex, u2ex] = ...
        torus_principal_curvatures(R, r, barycenter(V,F));
    gaussex = k1ex .* k2ex;
    A = 0.5*doublearea(V,F);
    
    %Principal curvatures are unreliable at umbilic points and at
    % points with zero Gaussian curvature
    unreliable = abs(k1ex-k2ex)<1e-6 | abs(k1ex)<1e-6 | abs(k2ex)<1e-6;
    u1s(unreliable,:) = 0;
    u2s(unreliable,:) = 0;
    u1ex(unreliable,:) = 0;
    u2ex(unreliable,:) = 0;
    
    k1err = abs(k1s - k1ex);
    k2err = abs(k2s - k2ex);
    gausserr = abs(gauss - gaussex);
    u1err = min(normrow(u1s-u1ex), normrow(u1s+u1ex));
    u2err = min(normrow(u2s-u2ex), normrow(u2s+u2ex));
    
    integratedk1(resol) = sum(A.*abs(k1s));
    integratedk2(resol) = sum(A.*abs(k2s));
    integratedgauss(resol) = sum(A.*abs(gauss));
    
    %Mesh stats
    m = max(E(:));
    aE = all_edges(F);
    posEdges = find(oE(:)>0);
    eE = nan(m,2); eE(E(posEdges),:) = aE(posEdges,:);
    edgeMps = 0.5 * (V(eE(:,1),:) + V(eE(:,2),:));
    hs(resol) = min(normrow(V(eE(:,1),:) - V(eE(:,2),:)));
    
    %Plot the Gauss curvature
    if(resol==numel(resols)-2)
        clf;
        plotf = gauss;
        t = tsurf(F,V, 'FaceVertexCData', plotf);
        cm = flipud(cbrewer('RdBu',500));
        colormap(cm);
        colorbar;
        axis equal;
        axis off;
        t.EdgeColor = 'none';
        
        sortedC = sort(plotf);
        ahl = max(abs(sortedC(floor(0.05*end))), ...
            abs(sortedC(floor(0.95*end))));
        caxis([-ahl, ahl]);
        
        shad = struct();
        shad.SpecularStrength = 0.2;
        shad.DiffuseStrength = 0.4;
        shad.AmbientStrength = 0.6;
        set(t,shad);
        camproj('persp');
        view(148.45, 22.88);
        l = light('Position',[1 1 2], 'Style', 'Infinite');
        s = add_shadow(t,l,'Color',[0.5 0.5 0.5],'Fade','infinite');
        delete(l)
        add_lights_o();
        apply_ambient_occlusion(t,'Factor',1,'AddLights',false,'SoftLighting',false);
        
        title('gauss curvature');
        
        set(gcf, 'Position',  [0, 0, 1280, 800]);
        name = sprintf('torus_gauss_plot.png');
        %imwrite(myaa('raw'), name);
        saveas(gcf, name);
    end
end


%Integrated error
clf;
loglog(hs, abs(integratedk1-sum(A.*abs(k1ex))), '-or', ...
    hs, abs(integratedk2-sum(A.*abs(k2ex))), '-og', ...
    hs, abs(integratedgauss-sum(A.*abs(gaussex))), '-ob', ...
    hs, hs, '--k', ...
    hs, hs.^2, '--k');
set(gca, 'Xdir', 'reverse');
xlabel('edge length h');
ylabel('Error, integrated');
title('Integrated error torus');
legend('k1', 'k2', 'gausscurvature', 'h', 'h^2');
saveas(gcf, ...
    sprintf('torusproperties_int.eps'), ...
    'epsc');

fprintf('\n');

