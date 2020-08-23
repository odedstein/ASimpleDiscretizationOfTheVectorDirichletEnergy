function z = project_onto_complex_pervertex_space(vec,X,N)
%PROJECT_ONTO_COMPLEX_PERVERTEX_SPACE Given 3D embedded vectors vec at
%points with a normal N and a tangent space basis X, project vec into the
%tangent space of each point parametrized by C

tangentVec = vec - dot(vec,N,2).*N;
tangentVecNorms = normrow(tangentVec);
tangentVecNormalized = tangentVec ./ tangentVecNorms;

tangX = X - dot(X,N,2).*N;
tangX = tangX ./ normrow(tangX);

vecDot = dot(tangX, tangentVecNormalized, 2);
vecCross = dot(N, cross(tangX, tangentVecNormalized, 2), 2);
vecAngle = atan2(vecCross, vecDot);

z = tangentVecNorms .* (cos(vecAngle) + 1i*sin(vecAngle));
z(tangentVecNorms < 1e-15) = 0;

end

