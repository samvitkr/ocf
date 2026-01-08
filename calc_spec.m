ret=180;
nu=1/ret;
tstep=400;
fn=sprintf('Sol000%04d0000000.h5',tstep)
baseDir = fullfile(getenv('MSIPROJECT'), 'xuanx004', 'ocf','ocf180');
fname   = fullfile(baseDir,fn);   % <-- edit if naming differs
datadir=fullfile('..','data');
% info=h5info(fname)
fprintf('Reading %s\n', fname);
u    = h5read(fname, '/u');
pex  = h5read(fname, '/pex');
pey  = h5read(fname, '/pey');
[Nx,Ny,Nz] = size(u);
kx=pex*[0:Nx/2-1,0,-Nx/2+1:-1];
ky=pey*[0:Ny/2-1,0,-Ny/2+1:-1];
dkx=kx(2)-kx(1);
dky=ky(2)-ky(1);
[Kx,Ky]=meshgrid(kx,ky);
Kx=Kx';
Ky=Ky';
load('../data/spectra.mat')
phivoz=phivoz./(dkx.*dky);
phiwoy=phiwoy./(dkx.*dky);

psc1=phivoz;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
pscm=real(pscm(1:Nx/2,1:Ny/2,:));
phivozm=pscm;

psc1=phiwoy;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
pscm=real(pscm(1:Nx/2,1:Ny/2,:));
phiwoym=pscm;

m=matfile('spectra2d.mat','Writable',true)
m.phivozm=phivozm;
m.phiwoym=phiwoym;
m.Kx =Kx;
m.Ky=Ky;
m.kx=kx;
m.ky=ky;
