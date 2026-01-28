clear
close all
Nx=640;
Ny=512;
Nz=192;

kcond=90;
tstart=300;
tend=1024;
step=2;

nf=(tend-tstart)/step+1;
nbin=256;
umin=-10;
umax=10;
wmin=-5;
wmax=5;
density=zeros(nbin,nbin);

%load('lambda_stats.mat')
for tstep=tstart:step:tend
        tstep
        fn=sprintf('Sol000%04d0000000.h5',tstep)
    % fnmat=sprintf('gradflux000%04d0000000.mat',tstep)
    baseDir = fullfile(getenv('MSIPROJECT'), 'xuanx004', 'ocf','ocf180');
    fname   = fullfile(baseDir,fn);   % <-- edit if naming differs
    % datadir=fullfile(getenv('MSIPROJECT'),'shared','kuma0458','open_channel_flow_180','data');
    % fnamemat=fullfile(datadir,fnmat);
    % info=h5info(fname)
    fprintf('Reading %s\n', fname);
    u    = h5read(fname, '/u');
    % v    = h5read(fname, '/v');
    w    = h5read(fname, '/w');
        % fvel=sprintf("../data/velfields_%07d.mat",time);
         % m=matfile(fvel);
         ul=reshape( u(:,:,kcond)-mean(u(:,:,kcond),"all") ,Nx*Ny,1 );
         wl=reshape( w(:,:,kcond) ,Nx*Ny,1 );
        uw=[ul,wl];
        [bandwidth,densityinst,ubin,wbin]=kde2d(uw,nbin,[umin wmin],[umax wmax]);
        density=density+densityinst;
end
density=density./nf;

fn=sprintf('../data/jointpdf_uw_k_%03d.mat',kcond);
mf=matfile(fn,"Writable",true);
mf.density=density;
mf.bandwidth=bandwidth;
mf.ubin=ubin;
mf.wbin=wbin;                                                                                                                                                         

q2 = (ubin<0).*(wbin>0);
q4 = (ubin>0).*(wbin<0);

mf.uq2 = sum(ubin.*density.*q2,'all')/sum(density.*q2,'all');
mf.wq2 = sum(wbin.*density.*q2,'all')/sum(density.*q2,'all');
mf.uq4 = sum(ubin.*density.*q4,'all')/sum(density.*q4,'all');
mf.wq4 = sum(wbin.*density.*q4,'all')/sum(density.*q4,'all');
