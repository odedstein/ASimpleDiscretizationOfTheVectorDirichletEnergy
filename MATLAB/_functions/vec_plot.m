function t = vec_plot(pts, vecs, varargin)
%VEC_PLOT Simple wrapper around quiver. Plots vectors vecs at pts;

v = 1;
scale = 1;
autoScale = 'on';
color = [1 1 1];
arrowRadius = 0.1;
logScale = false;
highlight = [];
highlightColor = [1 0 0];
head = true;
while v<=numel(varargin) && ischar(varargin{v})
    switch varargin{v}
        case 'Scale'
            v = v+1;
            assert(v<=numel(varargin));
            scale = varargin{v};
            autoScale = 'off';
        case 'Color'
            v = v+1;
            assert(v<=numel(varargin));
            color = varargin{v};
        case 'ArrowRadius'
            v = v+1;
            assert(v<=numel(varargin));
            arrowRadius = varargin{v};
        case 'LogScale'
            v = v+1;
            assert(v<=numel(varargin));
            logScale = varargin{v};
        case 'Highlight'
            v = v+1;
            assert(v<=numel(varargin));
            highlight = varargin{v};
        case 'HighlightColor'
            v = v+1;
            assert(v<=numel(varargin));
            highlightColor = varargin{v};
        case 'Head'
            v = v+1;
            assert(v<=numel(varargin));
            head = varargin{v};
        otherwise
            break;
    end
    v = v+1;
end

vecsNorm = normrow(vecs);
if(logScale)
    vecs = vecs ./ vecsNorm .* log(1+vecsNorm);
    vecsNorm = normrow(vecs);
end
%if(strcmp(autoScale,'on'))
%    vecs = vecs ./ max(vecsNorm);
%end

% t = quiver3(pts(:,1),pts(:,2),pts(:,3), ...
%     scale*vecs(:,1),scale*vecs(:,2),scale*vecs(:,3), ...
%     'MarkerEdgeColor', color, ...
%     'LineWidth', arrowRadius, ...
%     'MarkerSize', arrowHeadRadius, ...
%     'AutoScale', autoScale);

t = cell(size(pts,1),1);
shad = {};
shad.SpecularStrength = 0.05;
shad.DiffuseStrength = 0.05;
shad.AmbientStrength = 0.9;
for i=1:size(pts,1)
    if(sum(highlight==i)>0)
        plotCol = highlightColor;
    else
        plotCol = color;
    end
    if(head)
        t{i} = arrow3D(pts(i,:), scale*vecs(i,:), plotCol, 0.7, ...
            arrowRadius);
    else
        t{i} = arrow3D(pts(i,:), scale*vecs(i,:), plotCol, 1, arrowRadius);
    end
    set(t{i}, shad);
end

end

