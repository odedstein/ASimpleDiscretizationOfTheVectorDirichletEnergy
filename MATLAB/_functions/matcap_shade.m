function matcap_shade(mc,tsh,V,F)
%MATCAP_SHADE
% Given a matcap image mc = im2double(imread('image.png')), apply matcap
%  shading to the current figure and replace the figure with the shaded
%  image.
% mc: matcap image
% tsh: matlab trisurf handle
% V,F plotted mesh underlying the tsh handle.


N = per_vertex_normals(V,F);
pos = get(gcf, 'Position');
posl = pos;
imscale = 4;
posl(3) = imscale*posl(3);
posl(4) = imscale*posl(4);
figure(gcf);
axis off;
[AZ,EL] = view;
M = eye(3,4)*viewmtx(AZ,EL)'*eye(4,3);

%% THIS IS SIMPLER BUT QUANTIZES NORMALS
%tsh.FaceVertexCData = 0.5+0.5*N*M;
%I = im2double(getfield(getframe(gcf),'cdata'));
%tsh.FaceVertexCData = repmat([1 0 1], size(V,1),1);
%C = getfield(getframe(gcf),'cdata');
%A = ((C(:,:,1)==255) & (C(:,:,2)==0) & (C(:,:,3)==255) );
%PN = (I.*A)*2-1;

%% MORE COMPLICATED
[I,A,FI,B] = shader(tsh);
[II,~,IV] = find(FI(:));
N = N*M;
B1 = B(:,:,1);B2 = B(:,:,2);B3 = B(:,:,3);
PN = zeros(numel(B1),3);
PN(II,:) = normalizerow( ...
    N(F(IV,1),:).*B1(II) + N(F(IV,2),:).*B2(II) + N(F(IV,3),:).*B3(II));
PN = reshape(PN,size(B)).*A;

%% Plot
[Xmc,Ymc] = meshgrid(linspace(-1,1,size(mc,2)),linspace(-1,1,size(mc,1)));
Imc = [];
Imc(:,:,1) = interp2(Xmc,Ymc,mc(:,:,1),PN(:,:,1),-PN(:,:,2)).*A + I(:,:,1).*~A;
Imc(:,:,2) = interp2(Xmc,Ymc,mc(:,:,2),PN(:,:,1),-PN(:,:,2)).*A + I(:,:,2).*~A;
Imc(:,:,3) = interp2(Xmc,Ymc,mc(:,:,3),PN(:,:,1),-PN(:,:,2)).*A + I(:,:,3).*~A;

clf;
set(gcf,'Visible','Off')
set(gca,'Position',[0 0 1 1]);
set(gcf, 'Position',  posl)

imshow(Imc);
%imshow(Imc.*((0.5+0.5*I).*A+~A))

idat = im2double(getfield(getframe(gcf),'cdata'));
actimg = imresize(idat, 1/imscale, 'Antialiasing',true, 'Method','bicubic');

clf;
set(gcf,'Visible','On')
set(gca,'Position',[0 0 1 1]);
set(gcf, 'Position',  pos)
imshow(actimg);

end

