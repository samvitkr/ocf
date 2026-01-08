%% Read timestep 10 from OCF data and save into a .mat file
ret=180;
nu=1/ret;
tstart=300;
tend=1024;
step=2;
counter=0;
Nx=640;
Ny=512;
Nz=192;
phivoz = zeros(Nx,Ny,Nz);
phiwoy = zeros(Nx,Ny,Nz);
vozm = zeros(Nz,1);
woym = zeros(Nz,1);
viscm = zeros(Nz,1);
um = zeros(Nz,1);

uum = zeros(Nz,1);
vvm = zeros(Nz,1);
wwm = zeros(Nz,1);

uvm = zeros(Nz,1);
uwm = zeros(Nz,1);
vwm = zeros(Nz,1);

dudzm = zeros(Nz,1);

for tstep=tstart:step:tend
	counter=counter+1;
	fn=sprintf('Sol000%04d0000000.h5',tstep)
	fnmat=sprintf('gradflux000%04d0000000.mat',tstep)
	baseDir = fullfile(getenv('MSIPROJECT'), 'xuanx004', 'ocf','ocf180');
	fname   = fullfile(baseDir,fn);   % <-- edit if naming differs
	datadir=fullfile(getenv('MSIPROJECT'),'shared','kuma0458','open_channel_flow_180','data');
    	fnamemat=fullfile(datadir,fnmat);
	load(fnamemat)
	% info=h5info(fname)
	fprintf('Reading %s\n', fname);
	
	% --- Read datasets (no leading slash) ---
	zz   = h5read(fname, '/zz');
	u    = h5read(fname, '/u');
	v    = h5read(fname, '/v');
	pp   = h5read(fname, '/pp');
	pex  = h5read(fname, '/pex');
	pey  = h5read(fname, '/pey');
	time = h5read(fname, '/time');
	%dz   = h5read(fname, '/dz');
	%dzw   = h5read(fname, '/dzw');
	toc
	%% wave numbers and fft
	%[Nx,Ny,Nz] = size(u)
	
	kx=pex*[0:Nx/2-1,0,-Nx/2+1:-1]';
	ky=pey*[0:Ny/2-1,0,-Ny/2+1:-1]';
	

	fv=fft2(v);
	fw=fft2(wc);
	foz=fft2(dvdx-dudy);
	foy=fft2(dudz-dwdx);

	phivoz=phivoz+	fv.*(conj(foz));
	phiwoy=phiwoy+	fw.*(conj(foy));

	um= um+	squeeze(mean(u,[1 2]));

        uum=uum+squeeze(mean(u.*u,[1 2]));
        vvm=vvm+squeeze(mean(v.*v,[1 2]));
	wwm=wwm+squeeze(mean(wc.*wc,[1 2]));

        uwm=uwm+squeeze(mean(u.*wc,[1 2]));	
	uvm=uvm+squeeze(mean(u.*v,[1 2]));
        vwm=vwm+squeeze(mean(v.*wc,[1 2]));

	dudzm=	dudzm+	squeeze(mean(dudz,[1 2]));
	vozm=	vozm+	squeeze(mean(voz,[1 2]));
	woym=	woym+	squeeze(mean(woy,[1 2]));
	viscm=	viscm+	squeeze(mean(visc,[1 2]));
end


m=matfile('../data/spectra.mat','Writable',true)

m.phivoz=phivoz./(counter*Nx*Ny);
m.phiwoy=phiwoy./(counter*Nx*Ny);
m.zz=zz;
m.kx=kx;
m.ky=ky;

mp=matfile('../data/mean_profiles.mat','Writable',true)

mp.um=um./(counter);

mp.uum=uum./(counter);
mp.vvm=vvm./(counter);
mp.wwm=wwm./(counter);

mp.uvm=uvm./(counter);
mp.uwm=uwm./(counter);
mp.vwm=vwm./(counter);

mp.dudzm=um./(counter);
mp.vozm=vozm./(counter);
mp.woym=woym./(counter);
mp.viscm=viscm./(counter);
mp.zz=zz;
