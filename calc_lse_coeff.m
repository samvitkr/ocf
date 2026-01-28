clear all
clc

Nx=640;
Ny=512;
Nz=192;
kcond=64;

fn=sprintf('../data/velgrad_corr_k_%03d.mat',kcond);
%m=matfile(fn);
m=matfile(fn);
R1=zeros(Nx,Ny,Nz,14,'single');
R2=zeros(Nx,Ny,Nz,14,'single');

R1(:,:,:,1) 	= ifft2(m.phiuu,'symmetric');
R1(:,:,:,2) 	= ifft2(m.phiuv,'symmetric');
R1(:,:,:,3) 	= ifft2(m.phiuw,'symmetric');
R1(:,:,:,4) 	= ifft2(m.phiup,'symmetric');
R1(:,:,:,5) 	= ifft2(m.phiuf,'symmetric');
R1(:,:,:,6)	= ifft2(m.phiududx,'symmetric');
R1(:,:,:,7)	= ifft2(m.phiudvdx,'symmetric');
R1(:,:,:,8)	= ifft2(m.phiudwdx,'symmetric');
R1(:,:,:,9)	= ifft2(m.phiududy,'symmetric');
R1(:,:,:,10)	= ifft2(m.phiudvdy,'symmetric');
R1(:,:,:,11)	= ifft2(m.phiudwdy,'symmetric');
R1(:,:,:,12)	= ifft2(m.phiududz,'symmetric');
R1(:,:,:,13)	= ifft2(m.phiudvdz,'symmetric');
R1(:,:,:,14)	= ifft2(m.phiudwdz,'symmetric');

R2(:,:,:,1)       = ifft2(m.phiwu,'symmetric');
R2(:,:,:,2)       = ifft2(m.phiwv,'symmetric');
R2(:,:,:,3)       = ifft2(m.phiww,'symmetric');
R2(:,:,:,4)       = ifft2(m.phiwp,'symmetric');
R2(:,:,:,5)       = ifft2(m.phiwf,'symmetric');
R2(:,:,:,6)       = ifft2(m.phiwdudx,'symmetric');
R2(:,:,:,7)       = ifft2(m.phiwdvdx,'symmetric');
R2(:,:,:,8)       = ifft2(m.phiwdwdx,'symmetric');
R2(:,:,:,9)       = ifft2(m.phiwdudy,'symmetric');
R2(:,:,:,10)      = ifft2(m.phiwdvdy,'symmetric');
R2(:,:,:,11)      = ifft2(m.phiwdwdy,'symmetric');
R2(:,:,:,12)      = ifft2(m.phiwdudz,'symmetric');
R2(:,:,:,13)      = ifft2(m.phiwdvdz,'symmetric');
R2(:,:,:,14)      = ifft2(m.phiwdwdz,'symmetric');

clear m

L1=zeros(Nx,Ny,Nz,14,'single');
L2=zeros(Nx,Ny,Nz,14,'single');

uij=[R1(1,1,kcond,1),R1(1,1,kcond,3);
     R2(1,1,kcond,1),R2(1,1,kcond,3)];

M=inv(uij.');
L1 =single( R1.*M(1,1)+R2.*M(2,1));
L2 =single( R1.*M(1,2)+R2.*M(2,2));

fn=sprintf('../data/lse_coeff_k_%03d.mat',kcond);
mf=matfile(fn,"Writable",true);
mf.L1=real(L1);
mf.L2=real(L2);
