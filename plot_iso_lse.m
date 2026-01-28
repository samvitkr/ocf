close all
clear
Nx=640;
Ny=512;
ret=180;
Lx=2*pi/0.25;
Ly=2*pi*3/2;
x=[0:Nx-1]*Lx/Nx-Lx/2;
y=[0:Ny-1]*Ly/Ny-Ly/2;
mp=matfile('../data/mean_profiles.mat');
z=mp.zz;
[X,Y,Z]=(meshgrid(x,y,z));
X=single(X)*ret; 
Y=single(Y)*ret;
Z=single(Z)*ret;
k=64;


fq2=sprintf('../data/fieldsq2neg_lse_k_%03d.mat',k);
mq2=matfile(fq2);
fign=sprintf('isoq2pos_lse_k_%03d.fig',k);

% l2tot=fftshift(fftshift(mq2.lambda2tot,1),2);
% l2tot = permute(l2tot, [2 1 3]);

l2=fftshift(fftshift(mq2.lambda2,1),2);
l2 = permute(l2, [2 1 3]);

fig=figure;

%subplot(1,2,1)
hold on
p = patch(isosurface(X, Y, Z, l2, -0.001));
isonormals(X, Y, Z, l2, p); % Recalculate normals for better lighting
scatter3(0,0,z(k)*ret,100,'green','filled')
hold off
set(p, 'FaceColor', 'cyan', 'EdgeColor', 'none');
daspect([1 1 1]);
view(3);
axis tight;
camlight;
lighting gouraud;
xlabel('x^+')
ylabel('y^+')
zlabel('z^+')

% subplot(1,2,2)
% hold on
% p = patch(isosurface(X, Y, Z, l2, -1));
% isonormals(X, Y, Z, l2, p); % Recalculate normals for better lighting
% scatter3(0,0,z(k)*ret,100,'green','filled')
% hold off
% set(p, 'FaceColor', 'cyan', 'EdgeColor', 'none');
% daspect([1 1 1]);
% view(3);
% axis tight;
% camlight;
% lighting gouraud;
% xlabel('x^+')
% ylabel('y^+')
% zlabel('z^+')
% 
% saveas(fig,fign)