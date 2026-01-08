%% Read timestep 10 from OCF data and save into a .mat file
tic
ret=180;
nu=1/ret;
tstart=406;
tend=1024;
step=2;
for tstep=tstart:step:tend
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
    %dz   = h5read(fname, '/dz');
    %dzw   = h5read(fname, '/dzw');
    toc
    %% interpolate wc to cenn centers
    wc=permute(w,[3 1 2]);
    wc=interp1(zw(1:end-1),wc(1:end-1,:,:),zz);
    wc=permute(wc,[2 3 1]);
    toc
    %% wave numbers and fft
    [Nx,Ny,Nz] = size(u);
    kx=pex*[0:Nx/2-1,0,-Nx/2+1:-1]';
    ky=pey*[0:Ny/2-1,0,-Ny/2+1:-1]';
    %% x derivative
    fxu	=fft(u,[],1);
    fxv	=fft(v,[],1);
    fxw	=fft(wc,[],1);


    dudx	=ifft( fxu.*(1i.*kx),[],1,'symmetric');
    dvdx	=ifft( fxv.*(1i.*kx),[],1,'symmetric');
    dwdx	=ifft( fxw.*(1i.*kx),[],1,'symmetric');
    d2udx2	=ifft( fxu.*(-(kx).^2),[],1,'symmetric');
    toc
    clear fxu fxv fxw
    %% y derivative
    fyu     =fft(permute(u ,[2 1 3]),[],1);
    fyv     =fft(permute(v ,[2 1 3]),[],1);
    fyw     =fft(permute(wc,[2 1 3]),[],1);

    dudy    =permute(ifft( fyu.*(1i.*ky),[],1,'symmetric'),[2 1 3])	;
    dvdy    =permute(ifft( fyv.*(1i.*ky),[],1,'symmetric'),[2 1 3])	;
    dwdy    =permute(ifft( fyw.*(1i.*ky),[],1,'symmetric'),[2 1 3])	;
    d2udy2  =permute(ifft( fyu.*(-(ky).^2),[],1,'symmetric'),[2 1 3]);
    toc
    clear fyu fyv fyw



    %% finite difference for difference wall normal derivatives
    % dudz(i)= (1/dz(i+1) - 1/( dz(i+1)+dz(i) )*u(i+1) + (1/dz(i) -1/dz(i+1))*u(i) +
    % (1/(dz(i+1)+dz(i)) - 1/dz(i))*u(i-1);

    dudz=zeros(size(u));
    dvdz=zeros(size(u));
    dwdz=zeros(size(u));
    % dwcdz=zeros(size(u));
    d2udz2=zeros(size(u));
    dwdz(:,:,2:end) = permute(diff(permute(w,[3 2 1]),1,1)./diff(zw),[3 2 1]);
    dz=diff(zz);
    dudz(:,:,1)=(u(:,:,2)-u(:,:,1))./dz(1);
    dvdz(:,:,1)=(v(:,:,2)-v(:,:,1))./dz(1);



    %for i =1:length(dzz)-1
    %    a = 1/dz(i+1) - 1/( dz(i+1)+dz(i) );
    %    b = 1/dz(i) - 1/dz(i+1);
    %    c = 1/( dz(i+1)+dz(i)) - 1/dz(i);
    %
    %    a2=2/((dz(i+1)+dz(i))*dz(i+1));
    %    b2=-2/(dz(i)*dz(i+1));
    %    c2=2/((dz(i+1)+dz(i))*dz(i));
    %
    %    k=i+1;
    %    dudz(:,:,k)= a.*u(:,:,k+1) +b.*u(:,:,k) +c.*u(:,:,k-1);
    %    dvdz(:,:,k)= a.*v(:,:,k+1) +b.*v(:,:,k) +c.*v(:,:,k-1);
    %    % dwcdz(:,:,k)= a.*wc(:,:,k+1) +b.*wc(:,:,k) +c.*wc(:,:,k-1);
    %    d2udz2(:,:,k)= a2.*u(:,:,k+1) +b2.*u(:,:,k) +c2.*u(:,:,k-1);
    %end

    nz=Nz;
    i      = 1:(nz-2);
    dz_i   = dz(i);
    dz_ip1 = dz(i+1);
    sumdz  = dz_i + dz_ip1;

    a  = 1./dz_ip1 - 1./sumdz;
    b  = 1./dz_i   - 1./dz_ip1;
    c  = 1./sumdz  - 1./dz_i;

    % reshape for implicit expansion along 3rd dim
    a  = reshape(a, 1, 1, []);
    b  = reshape(b, 1, 1, []);
    c  = reshape(c, 1, 1, []);

    % Apply to u, v
    dudz(:,:,2:nz-1) =a.*u(:,:,3:nz)+b.*u(:,:,2:nz-1)+c.* u(:,:,1:nz-2);
    dvdz(:,:,2:nz-1) =a.*v(:,:,3:nz)+b.*v(:,:,2:nz-1)+c.* v(:,:,1:nz-2);
    % dwcdz(:,:,2:nz-1) =a.*wc(:,:,3:nz)+b.*wc(:,:,2:nz-1)+c.*wc(:,:,1:nz-2);
    % Second derivative coefficients
    a2 =  2 ./ (sumdz .* dz_ip1);
    b2 = -2 ./ (dz_i   .* dz_ip1);
    c2 =  2 ./ (sumdz .* dz_i);
    a2 = reshape(a2, 1, 1, []);
    b2 = reshape(b2, 1, 1, []);
    c2 = reshape(c2, 1, 1, []);
    d2udz2(:,:,2:nz-1) =a2.*u(:,:,3:nz)+b2.*u(:,:,2:nz-1)+c2.*u(:,:,1:nz-2);

    dudz(:,:,end)=(u(:,:,end)-u(:,:,end-1))./dz(end);
    dvdz(:,:,end)=(v(:,:,end)-v(:,:,end-1))./dz(end);
    d2udz2(:,:,1)=d2udz2(:,:,2);
    d2udz2(:,:,end)=d2udz2(:,:,end-1);
    %%
    wc     = single(wc);
    dudx   = single(dudx);
    dvdx   = single(dvdx);
    dwdx   = single(dwdx);
    d2udx2 = single(d2udx2);

    dudy   = single(dudy);
    dvdy   = single(dvdy);
    dwdy   = single(dwdy);
    d2udy2 = single(d2udy2);

    dudz   = single(dudz);
    dvdz   = single(dvdz);
    dwdz   = single(dwdz);
    d2udz2 = single(d2udz2);

    oz=dvdx-dudy;
    oy=dudz-dwdx;
    voz=single(v.*oz);
    woy=single(w.*oy);
    visc =single(nu*( d2udx2+d2udy2+d2udz2));

    toc
    % --- Save everything to MAT file ---
    % outFile =  ('ocf_test.mat');
    outFile=fnamemat;
    save(outFile,'zz','wc',...
        'dudx','dvdx','dwdx','d2udx2',...
        'dudy','dvdy','dwdy','d2udy2',...
        'dudz','dvdz','dwdz','d2udz2',...
        'voz','woy','visc','-v7.3')
    %save(outFile, 'zz','zw','u','v','w','wc','pp','pex','pey','time','-v7.3');
    % -v7.3 needed for large 3D arrays
    toc
    fprintf('computed gradients for %0.3f\n', time);
end
