clear 
Nx=640;
Ny=512;
Nz=192;
filpos=single(zeros(Nx,Ny,Nz));
filneg=single(zeros(Nx,Ny,Nz));
mt=matfile('../data/siren_filters.mat','Writable',true);
fp=filpos;
fp(Nx/2+1:end,Ny/2+1:end,:) = mt.mask_pos;
filpos = fp+flip(fp,1)+flip(fp,2)+flip(flip(fp,2),1);
mt.filpos=fftshift(filpos);

fn=filneg;
fn(Nx/2+1:end,Ny/2+1:end,:) = mt.mask_neg;
filneg = fn+flip(fn,1)+flip(fn,2)+flip(flip(fn,2),1);
mt.filneg=fftshift(filneg);
