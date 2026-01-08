close all
clear
ret=180;
nu=1/ret;
tstart=300;
tend=1024;
step=2;
counter=0;
Nx=640;
Ny=512;
Nz=192;
jcond=64;

phiuu=zeros(Nx,Ny,Nz);
phivv=zeros(Nx,Ny,Nz);
phiuv=zeros(Nx,Ny,Nz);
phivu=zeros(Nx,Ny,Nz);
phiuw=zeros(Nx,Ny,Nz);
phivw=zeros(Nx,Ny,Nz);
phiup=zeros(Nx,Ny,Nz);
phivp=zeros(Nx,Ny,Nz);

phiududx=zeros(Nx,Ny,Nz);
phiudvdx=zeros(Nx,Ny,Nz);
phiudwdx=zeros(Nx,Ny,Nz);
phiududy=zeros(Nx,Ny,Nz);
phiudvdy=zeros(Nx,Ny,Nz);
phiudwdy=zeros(Nx,Ny,Nz);
phiududz=zeros(Nx,Ny,Nz);
phiudvdz=zeros(Nx,Ny,Nz);
phiudwdz=zeros(Nx,Ny,Nz);
phivdudx=zeros(Nx,Ny,Nz);
phivdvdx=zeros(Nx,Ny,Nz);
phivdwdx=zeros(Nx,Ny,Nz);
phivdudy=zeros(Nx,Ny,Nz);
phivdvdy=zeros(Nx,Ny,Nz);
phivdwdy=zeros(Nx,Ny,Nz);
phivdudz=zeros(Nx,Ny,Nz);
phivdvdz=zeros(Nx,Ny,Nz);
phivdwdz=zeros(Nx,Ny,Nz);

phiuvoz=zeros(Nx,Ny,Nz);
phivvoz=zeros(Nx,Ny,Nz);
phiuwoy=zeros(Nx,Ny,Nz);
phivwoy=zeros(Nx,Ny,Nz);

nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
for time=tstart:step:tend
        time
	fn=sprintf('Sol000%04d0000000.h5',tstep)
	fnmat=sprintf('gradflux000%04d0000000.mat',tstep)
	baseDir = fullfile(getenv('MSIPROJECT'), 'xuanx004', 'ocf','ocf180');
	fname   = fullfile(baseDir,fn);   % <-- edit if naming differs
	datadir=fullfile(getenv('MSIPROJECT'),'shared','kuma0458','open_channel_flow_180','data');
	fnamemat=fullfile(datadir,fnmat);
	% info=h5info(fname)
	fprintf('Reading %s\n', fname);
	
	% --- Read datasets (no leading slash) ---
	zz   = h5read(fname, '/zz');
	zw   = h5read(fname, '/zw');
	u    = h5read(fname, '/u');
	v    = h5read(fname, '/v');
	w    = h5read(fname, '/w');
	pp   = h5read(fname, '/pp');
	pex  = h5read(fname, '/pex');
	pey  = h5read(fname, '/pey');
	time = h5read(fname, '/time');



        fvel=sprintf("../data/velfields_%07d.mat",time);
        m=matfile(fvel);
	fvelg=sprintf("../data/velgrad_%07d.mat",time);
        mg=matfile(fvelg);
	ft=sprintf("../data/transferfields_%07d.mat",time);
        mt=matfile(ft);
	viscF=fft2(mt.visc(:,:,Ny/2+1:end))./(Nz*Nx);
%m=matfile('velfields_0050000.mat');
	ufj=m.uFourier(:,:,jcond);
	vfj=m.vFourier(:,:,jcond);
%	wfj=m.wFourier(:,:,jcond);
	ufj(1,1)=0;
	vfj(1,1)=0;
%	wfj(1,1)=0;

	phiuu=phiuu+conj(ufj).*m.uFourier(:,:,Ny/2+1:end);
	phivv=phivv+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
%	phiww=phiww+conj(wfj).*m.wFourier(:,:,Ny/2+1:end);
	phiuv=phiuv+conj(ufj).*m.vFourier(:,:,Ny/2+1:end);
	phivu=phivu+conj(vfj).*m.uFourier(:,:,Ny/2+1:end);
	phiuw=phiuw+conj(ufj).*m.wFourier(:,:,Ny/2+1:end);
%	phiwu=phiwu+conj(wfj).*m.uFourier(:,:,Ny/2+1:end);
	phivw=phivw+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
%	phiwv=phiwv+conj(wfj).*m.vFourier(:,:,Ny/2+1:end);

        phiududx=phiududx+conj(ufj).*fft2(mg.dudx(:,:,Ny/2+1:end))./(Nx*Nz);
	phiudvdx=phiudvdx+conj(ufj).*fft2(mg.dvdx(:,:,Ny/2+1:end))./(Nx*Nz);
	phiudwdx=phiudwdx+conj(ufj).*fft2(mg.dwdx(:,:,Ny/2+1:end))./(Nx*Nz);
	phiududy=phiududy+conj(ufj).*fft2(mg.dudy(:,:,Ny/2+1:end))./(Nx*Nz);
	phiudvdy=phiudvdy+conj(ufj).*fft2(mg.dvdy(:,:,Ny/2+1:end))./(Nx*Nz);
	phiudwdy=phiudwdy+conj(ufj).*fft2(mg.dwdy(:,:,Ny/2+1:end))./(Nx*Nz);
	phiududz=phiududz+conj(ufj).*fft2(mg.dudz(:,:,Ny/2+1:end))./(Nx*Nz);
	phiudvdz=phiudvdz+conj(ufj).*fft2(mg.dvdz(:,:,Ny/2+1:end))./(Nx*Nz);
	phiudwdz=phiudwdz+conj(ufj).*fft2(mg.dwdz(:,:,Ny/2+1:end))./(Nx*Nz);	
	phiufx=phiufx+conj(ufj).*viscF;

        phivdudx=phivdudx+conj(vfj).*fft2(mg.dudx(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdvdx=phivdvdx+conj(vfj).*fft2(mg.dvdx(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdwdx=phivdwdx+conj(vfj).*fft2(mg.dwdx(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdudy=phivdudy+conj(vfj).*fft2(mg.dudy(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdvdy=phivdvdy+conj(vfj).*fft2(mg.dvdy(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdwdy=phivdwdy+conj(vfj).*fft2(mg.dwdy(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdudz=phivdudz+conj(vfj).*fft2(mg.dudz(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdvdz=phivdvdz+conj(vfj).*fft2(mg.dvdz(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdwdz=phivdwdz+conj(vfj).*fft2(mg.dwdz(:,:,Ny/2+1:end))./(Nx*Nz);
	phivfx=phivfx+conj(vfj).*viscF;

end

phiuu=phiuu./nf;
phivv=phivv./nf;
%phiww=phiww./nf;
phiuv=phiuv./nf;
phivu=phivu./nf;
phiuw=phiuw./nf;
%phiwu=phiwu./nf;
phivw=phivw./nf;
%phiwv=phiwv./nf;

phiududx=phiududx./nf;
phiudvdx=phiudvdx./nf;
phiudwdx=phiudwdx./nf;
phiududy=phiududy./nf;
phiudvdy=phiudvdy./nf;
phiudwdy=phiudwdy./nf;
phiududz=phiududz./nf;
phiudvdz=phiudvdz./nf;
phiudwdz=phiudwdz./nf;
phiufx=phiufx./nf;

phivdudx=phivdudx./nf;
phivdvdx=phivdvdx./nf;
phivdwdx=phivdwdx./nf;
phivdudy=phivdudy./nf;
phivdvdy=phivdvdy./nf;
phivdwdy=phivdwdy./nf;
phivdudz=phivdudz./nf;
phivdvdz=phivdvdz./nf;
phivdwdz=phivdwdz./nf;
phivfx=phivfx./nf;

Ruu=ifft2(phiuu*(Nz*Nx),'symmetric');
Rvv=ifft2(phivv*(Nz*Nx),'symmetric');
%Rww=ifft2(phiww*(Nz*Nx),'symmetric');
Ruv=ifft2(phiuv*(Nz*Nx),'symmetric');
Rvu=ifft2(phivu*(Nz*Nx),'symmetric');
Ruw=ifft2(phiuw*(Nz*Nx),'symmetric');
%Rwu=ifft2(phiwu*(Nz*Nx),'symmetric');
Rvw=ifft2(phivw*(Nz*Nx),'symmetric');
%Rwv=ifft2(phiwv*(Nz*Nx),'symmetric');

Rududx=ifft2(phiududx*(Nz*Nx),'symmetric');
Rudvdx=ifft2(phiudvdx*(Nz*Nx),'symmetric');
Rudwdx=ifft2(phiudwdx*(Nz*Nx),'symmetric');
Rududy=ifft2(phiududy*(Nz*Nx),'symmetric');
Rudvdy=ifft2(phiudvdy*(Nz*Nx),'symmetric');
Rudwdy=ifft2(phiudwdy*(Nz*Nx),'symmetric');
Rududz=ifft2(phiududz*(Nz*Nx),'symmetric');
Rudvdz=ifft2(phiudvdz*(Nz*Nx),'symmetric');
Rudwdz=ifft2(phiudwdz*(Nz*Nx),'symmetric');

Rvdudx=ifft2(phivdudx*(Nz*Nx),'symmetric');
Rvdvdx=ifft2(phivdvdx*(Nz*Nx),'symmetric');
Rvdwdx=ifft2(phivdwdx*(Nz*Nx),'symmetric');
Rvdudy=ifft2(phivdudy*(Nz*Nx),'symmetric');
Rvdvdy=ifft2(phivdvdy*(Nz*Nx),'symmetric');
Rvdwdy=ifft2(phivdwdy*(Nz*Nx),'symmetric');
Rvdudz=ifft2(phivdudz*(Nz*Nx),'symmetric');
Rvdvdz=ifft2(phivdvdz*(Nz*Nx),'symmetric');
Rvdwdz=ifft2(phivdwdz*(Nz*Nx),'symmetric');

Rufx=ifft2(phiufx*(Nz*Nx),'symmetric');
Rvfx=ifft2(phivfx*(Nz*Nx),'symmetric');

fn=sprintf('../data/velgrad_corr_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.Ruu=Ruu;
mf.Rvv=Rvv;
%mf.Rww=Rww;
mf.Ruv=Ruv;
mf.Rvu=Rvu;
mf.Ruw=Ruw;
%mf.Rwu=Rwu;
mf.Rvw=Rvw;
%mf.Rwv=Rwv;
%mf.yCheb=yCheb(Ny/2+1:end);
%mf.j=jc;

mf.Rududx=Rududx;
mf.Rudvdx=Rudvdx;
mf.Rudwdx=Rudwdx;
mf.Rududy=Rududy;
mf.Rudvdy=Rudvdy;
mf.Rudwdy=Rudwdy;
mf.Rududz=Rududz;
mf.Rudvdz=Rudvdz;
mf.Rudwdz=Rudwdz;

mf.Rvdudx=Rvdudx;
mf.Rvdvdx=Rvdvdx;
mf.Rvdwdx=Rvdwdx;
mf.Rvdudy=Rvdudy;
mf.Rvdvdy=Rvdvdy;
mf.Rvdwdy=Rvdwdy;
mf.Rvdudz=Rvdudz;
mf.Rvdvdz=Rvdvdz;
mf.Rvdwdz=Rvdwdz;

mf.Rufx=Rufx;
mf.Rvfx=Rvfx;
