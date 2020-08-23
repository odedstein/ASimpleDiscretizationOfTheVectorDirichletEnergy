function avgd = average_per_face_onto_vertices(V,F,z)
%AVERADE_PER_FACE_ONTO_VERTICES Given faces F, average the
%per-face scalar data z (column vector) onto the vertices v

n = size(V,1);

A = 0.5*doublearea(V,F);

is = [F(:,1); F(:,2); F(:,3)];
js = ones(size(is));
ks = [A; A; A];
perVertA = full(sparse(is, js, ks, n, 1));

ks = [A.*z; A.*z; A.*z];
avgd = full(sparse(is, js, ks, n, 1));
avgd = avgd ./ perVertA;

