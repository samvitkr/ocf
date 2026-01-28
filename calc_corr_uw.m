close all
clear
ret=180;
nu=1/ret;
tstart=300;
%tend=310;
tend=1024;
step=2;
counter=0;
Nx=640;
Ny=512;
Nz=192;
kcond=100;

phiuu=complex(zeros(Nx,Ny,Nz,'single'));
phiuv=complex(zeros(Nx,Ny,Nz,'single'));
phiuw=complex(zeros(Nx,Ny,Nz,'single'));
phiup=complex(zeros(Nx,Ny,Nz,'single'));
phiuf=complex(zeros(Nx,Ny,Nz,'single'));

phiwu=complex(zeros(Nx,Ny,Nz,'single'));
phiwv=complex(zeros(Nx,Ny,Nz,'single'));
phiww=complex(zeros(Nx,Ny,Nz,'single'));
phiwp=complex(zeros(Nx,Ny,Nz,'single'));
phiwf=complex(zeros(Nx,Ny,Nz,'single'));

phiududx=complex(zeros(Nx,Ny,Nz,'single'));
phiudvdx=complex(zeros(Nx,Ny,Nz,'single'));
phiudwdx=complex(zeros(Nx,Ny,Nz,'single'));
phiududy=complex(zeros(Nx,Ny,Nz,'single'));
phiudvdy=complex(zeros(Nx,Ny,Nz,'single'));
phiudwdy=complex(zeros(Nx,Ny,Nz,'single'));
phiududz=complex(zeros(Nx,Ny,Nz,'single'));
phiudvdz=complex(zeros(Nx,Ny,Nz,'single'));
phiudwdz=complex(zeros(Nx,Ny,Nz,'single'));

phiwdudx=complex(zeros(Nx,Ny,Nz,'single'));
phiwdvdx=complex(zeros(Nx,Ny,Nz,'single'));
phiwdwdx=complex(zeros(Nx,Ny,Nz,'single'));
phiwdudy=complex(zeros(Nx,Ny,Nz,'single'));
phiwdvdy=complex(zeros(Nx,Ny,Nz,'single'));
phiwdwdy=complex(zeros(Nx,Ny,Nz,'single'));
phiwdudz=complex(zeros(Nx,Ny,Nz,'single'));
phiwdvdz=complex(zeros(Nx,Ny,Nz,'single'));
phiwdwdz=complex(zeros(Nx,Ny,Nz,'single'));


nf=(tend-tstart)/step+1;
%load('lambda_stats.mat')

for tstep=tstart:step:tend
        tstep
	fn=sprintf('Sol000%04d0000000.h5',tstep)
	fnmat=sprintf('gradflux000%04d0000000.mat',tstep)
	baseDir = fullfile(getenv('MSIPROJECT'), 'xuanx004', 'ocf','ocf180');
	fname   = fullfile(baseDir,fn);   % <-- edit if naming differs
	datadir=fullfile(getenv('MSIPROJECT'),'shared','kuma0458','open_channel_flow_180','data');
	fnamemat=fullfile(datadir,fnmat);
	load(fnamemat);
	% info=h5info(fname)
	fprintf('Reading %s\n', fname);
	
	% --- Read datasets (no leading slash) ---
	u    = h5read(fname, '/u');
	v    = h5read(fname, '/v');
	p   = h5read(fname, '/pp');

	ufj=fft2(u(:,:,kcond));
	wfj=fft2(wc(:,:,kcond));
	ufj(1,1)=0;
	wfj(1,1)=0;

	phiuu=phiuu+conj(ufj).*fft2(u)./(Nx*Ny);
	phiuv=phiuv+conj(ufj).*fft2(v)./(Nx*Ny);
        phiuw=phiuw+conj(ufj).*fft2(wc)./(Nx*Ny);	
	phiup=phiup+conj(ufj).*fft2(p)./(Nx*Ny);
	phiuf=phiuf+conj(ufj).*fft2(visc)./(Nx*Ny);

	phiwu=phiwu+conj(wfj).*fft2(u)./(Nx*Ny);
	phiwv=phiwv+conj(wfj).*fft2(v)./(Nx*Ny);
        phiww=phiww+conj(wfj).*fft2(wc)./(Nx*Ny);
        phiwp=phiwp+conj(wfj).*fft2(p)./(Nx*Ny);
	phiwf=phiwf+conj(wfj).*fft2(visc)./(Nx*Ny);

        phiududx=phiududx+conj(ufj).*fft2(dudx)./(Nx*Ny);
	phiudvdx=phiudvdx+conj(ufj).*fft2(dvdx)./(Nx*Ny);
	phiudwdx=phiudwdx+conj(ufj).*fft2(dwdx)./(Nx*Ny);
	phiududy=phiududy+conj(ufj).*fft2(dudy)./(Nx*Ny);
	phiudvdy=phiudvdy+conj(ufj).*fft2(dvdy)./(Nx*Ny);
	phiudwdy=phiudwdy+conj(ufj).*fft2(dwdy)./(Nx*Ny);
	phiududz=phiududz+conj(ufj).*fft2(dudz)./(Nx*Ny);
	phiudvdz=phiudvdz+conj(ufj).*fft2(dvdz)./(Nx*Ny);
	phiudwdz=phiudwdz+conj(ufj).*fft2(dwdz)./(Nx*Ny);	

	phiwdudx=phiwdudx+conj(wfj).*fft2(dudx)./(Nx*Ny);
        phiwdvdx=phiwdvdx+conj(wfj).*fft2(dvdx)./(Nx*Ny);
        phiwdwdx=phiwdwdx+conj(wfj).*fft2(dwdx)./(Nx*Ny);
        phiwdudy=phiwdudy+conj(wfj).*fft2(dudy)./(Nx*Ny);
        phiwdvdy=phiwdvdy+conj(wfj).*fft2(dvdy)./(Nx*Ny);
        phiwdwdy=phiwdwdy+conj(wfj).*fft2(dwdy)./(Nx*Ny);
        phiwdudz=phiwdudz+conj(wfj).*fft2(dudz)./(Nx*Ny);
        phiwdvdz=phiwdvdz+conj(wfj).*fft2(dvdz)./(Nx*Ny);
        phiwdwdz=phiwdwdz+conj(wfj).*fft2(dwdz)./(Nx*Ny);

end

phiuu=phiuu./nf;
phiuv=phiuv./nf;
phiuw=phiuw./nf;
phiup=phiup./nf;
phiuf=phiuf./nf;

phiwu=phiwu./nf;
phiwv=phiwv./nf;
phiww=phiww./nf;
phiwp=phiwp./nf;
phiwf=phiwf./nf;

phiududx=phiududx./nf;
phiudvdx=phiudvdx./nf;
phiudwdx=phiudwdx./nf;
phiududy=phiududy./nf;
phiudvdy=phiudvdy./nf;
phiudwdy=phiudwdy./nf;
phiududz=phiududz./nf;
phiudvdz=phiudvdz./nf;
phiudwdz=phiudwdz./nf;

phiwdudx=phiwdudx./nf;
phiwdvdx=phiwdvdx./nf;
phiwdwdx=phiwdwdx./nf;
phiwdudy=phiwdudy./nf;
phiwdvdy=phiwdvdy./nf;
phiwdwdy=phiwdwdy./nf;
phiwdudz=phiwdudz./nf;
phiwdvdz=phiwdvdz./nf;
phiwdwdz=phiwdwdz./nf;

%Ruu=ifft2(phiuu,'symmetric');
%Rvv=ifft2(phivv,'symmetric');
%Ruv=ifft2(phiuv,'symmetric');
%Rvu=ifft2(phivu,'symmetric');
%Ruw=ifft2(phiuw,'symmetric');
%Rvw=ifft2(phivw,'symmetric');
%
%Rududx=ifft2(phiududx,'symmetric');
%Rudvdx=ifft2(phiudvdx,'symmetric');
%Rudwdx=ifft2(phiudwdx,'symmetric');
%Rududy=ifft2(phiududy,'symmetric');
%Rudvdy=ifft2(phiudvdy,'symmetric');
%Rudwdy=ifft2(phiudwdy,'symmetric');
%Rududz=ifft2(phiududz,'symmetric');
%Rudvdz=ifft2(phiudvdz,'symmetric');
%Rudwdz=ifft2(phiudwdz,'symmetric');
%
%Rvdudx=ifft2(phivdudx,'symmetric');
%Rvdvdx=ifft2(phivdvdx,'symmetric');
%Rvdwdx=ifft2(phivdwdx,'symmetric');
%Rvdudy=ifft2(phivdudy,'symmetric');
%Rvdvdy=ifft2(phivdvdy,'symmetric');
%Rvdwdy=ifft2(phivdwdy,'symmetric');
%Rvdudz=ifft2(phivdudz,'symmetric');
%Rvdvdz=ifft2(phivdvdz,'symmetric');
%Rvdwdz=ifft2(phivdwdz,'symmetric');
%
%Rufx=ifft2(phiufx*(Nz*Nx),'symmetric');
%Rvfx=ifft2(phivfx*(Nz*Nx),'symmetric');

fn=sprintf('../data/velgrad_corr_k_%03d.mat',kcond);
mf=matfile(fn,"Writable",true);
mf.phiuu=phiuu;
mf.phiuv=phiuv;
mf.phiuw=phiuw;
mf.phiup=phiup;
mf.phiuf=phiuf;

mf.phiwu=phiwu;
mf.phiwv=phiwv;
mf.phiww=phiww;
mf.phiwp=phiwp;
mf.phiwf=phiwf;

mf.phiududx=phiududx;
mf.phiudvdx=phiudvdx;
mf.phiudwdx=phiudwdx;
mf.phiududy=phiududy;
mf.phiudvdy=phiudvdy;
mf.phiudwdy=phiudwdy;
mf.phiududz=phiududz;
mf.phiudvdz=phiudvdz;
mf.phiudwdz=phiudwdz;

mf.phiwdudx=phiwdudx;
mf.phiwdvdx=phiwdvdx;
mf.phiwdwdx=phiwdwdx;
mf.phiwdudy=phiwdudy;
mf.phiwdvdy=phiwdvdy;
mf.phiwdwdy=phiwdwdy;
mf.phiwdudz=phiwdudz;
mf.phiwdvdz=phiwdvdz;
mf.phiwdwdz=phiwdwdz;

