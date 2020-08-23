function t = scatter_plot(pts, moreargs)
%SCATTER_PLOT Simple wrapper around scatter. Plots scattered data pts;

t = scatter3(pts(:,1), pts(:,2), pts(:,3));
if(nargin>2)
    set(t, moreargs{:})
end

end

