addpath ..;
path_handling;

% These are the numbers of boundaries we will have
bdns = 6 * 2.^linspace(1,4,6);
bdnsIrreg = 6 * 2.^linspace(2,7,6);

% These are the bdry conditions and the exact solution of the problem
bdrycond = @(x) [sin( 2*atan2(x(:,2),x(:,1)) ).^2, ...
    cos( 2*atan2(x(:,2),x(:,1)) ).^2, ...
    zeros(size(x,1),1)];
% Solution from mathematica
soln = @(x) 0.5*[1-x(:,1).^4+6*(x(:,1).^2).*(x(:,2).^2) - x(:,2).^4, ...
    1+x(:,1).^4-6*(x(:,1).^2).*(x(:,2).^2) + x(:,2).^4 ...
    zeros(size(x,1),1)];


% Set up variables we care about
hsReg = nan(1,numel(bdns));
l2errsCrReg = nan(1,numel(bdns));
regratioReg = nan(1,numel(bdns));
hsIrreg = nan(1,numel(bdns));
l2errsCrIrreg = nan(1,numel(bdns));
regratioIrreg = nan(1,numel(bdns));

% Do convergence experiment
regs = {'reg', 'irreg'};
for i=1:numel(bdns)
    for regind = 1:numel(regs)
        regtype = regs{regind};
        fprintf('doing %d / %d\n', i, numel(bdns));
        
        % Construct disk
        if(strcmp(regtype,'reg'))
            [V,F,b] = disk(bdns(i));
            regratioReg(i) = max(circumradius(V,F)./inradius(V,F));
        elseif(strcmp(regtype,'irreg'))
            [V,F,b] = disk(bdnsIrreg(i),true);
            regratioIrreg(i) = max(circumradius(V,F)./inradius(V,F));
        end
        V = [V zeros(size(V,1),1)];
        if(any(is_vertex_nonmanifold(F)))
            error('nonmanifold mesh');
        end
        
        if((i==1&&strcmp(regtype,'irreg')) || i==2 || ...
                (i==3&&strcmp(regtype,'reg')))
                figure;
                tsurf(F,V, 'LineWidth',3);
                alpha 0;
                hold on;
                axis equal;
                axis off;
            if(strcmp(regtype,'reg'))
                % Plot both reg and irreg mesh
                title('reg. mesh');
            elseif(strcmp(regtype,'irreg'))
                title('irreg. mesh');
            end
            
            %set(gcf, 'Position',  [0, 0, 3840, 2160])
            %imwrite(myaa('raw'),['wireframe-' regtype '-' int2str(i) '.png']);
            saveas(gcf, ['wireframe-' regtype '-' int2str(i) '.png']);
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
        % Average vector magnitudes onto vertices
        crMag = average_per_edge_onto_vertices(eE, normrow(crPlotVec));
        
        % Measure L2 error
        M = crouzeix_raviart_massmatrix(V,F);
        edgeErr = normrow(crPlotVec - soln(edgeMps));
        if(strcmp(regtype,'reg'))
            l2errsCrReg(i) = sqrt(edgeErr'*M*edgeErr);
        elseif(strcmp(regtype,'irreg'))
            l2errsCrIrreg(i) = sqrt(edgeErr'*M*edgeErr);
        end
        
        
        % Save max edge length
        if(strcmp(regtype,'reg'))
            hsReg(i) = max(normrow(V(eE(:,1),:) - V(eE(:,2),:)));
        elseif(strcmp(regtype,'irreg'))
            hsIrreg(i) = max(normrow(V(eE(:,1),:) - V(eE(:,2),:)));
        end
        
    end
end


% Error plot
figure;
loglog(hsReg, l2errsCrReg, '-or', ...
    hsIrreg, l2errsCrIrreg, '-xr', ...
    hsReg, hsReg, '--k', ...
    hsReg, hsReg.^2, '--k');
set(gca, 'xdir', 'reverse' )
xlabel('h');
ylabel('L2 error');
title('Bdry value problem with known solution, convergence.');
legend('CR (reg)', ...
    'CR (irreg)', ...
    'h', 'h^2');
saveas(gcf, 'triangleregularity-convergence-plot.eps', 'epsc');

