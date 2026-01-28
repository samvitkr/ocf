Nx=2048;
Ny=512;
Nz=1536;
nproc=6;

tstart=31;
tend=32;
for time=tstart:tend
	S_11 = single(zeros(Ny,Nx,Nz));
	S_12 = single(zeros(Ny,Nx,Nz));
	S_13 = single(zeros(Ny,Nx,Nz));
	S_22 = single(zeros(Ny,Nx,Nz));
	S_23 = single(zeros(Ny,Nx,Nz));
	S_33 = single(zeros(Ny,Nx,Nz));
	
	O_21 = single(zeros(Ny,Nx,Nz));
	O_13 = single(zeros(Ny,Nx,Nz));
	O_32 = single(zeros(Ny,Nx,Nz));

	lambda2 = single(zeros(Ny,Nx,Nz));
	Q = single(zeros(Ny,Nx,Nz));
	
	fvgx=sprintf('velgradx_%03d.mat',time);
	fvgy=sprintf('velgrady_%03d.mat',time);
	fvgz=sprintf('velgradz_%03d.mat',time);
	fo=sprintf('vort_%03d',time);
	

	mvgx=matfile(fvgx);
	mvgy=matfile(fvgy);
	mvgz=matfile(fvgz);
	mo=matfile(fo);

	time
	tic
	S_11=mvgx.dudx;
	S_12=0.5*( mvgy.dudy + mvgx.dvdx );
	S_13=0.5*( mvgz.dudz + mvgx.dwdx );
	S_22=mvgy.dvdy;
	S_23=0.5*( mvgy.dwdy + mvgz.dvdz );
	S_33=mvgz.dwdz;
 		
	O_21 = 0.5*mo.omegaz(:,:,:);
	O_13 = 0.5*mo.omegay(:,:,:);
	O_32 = 0.5*mo.omegax(:,:,:);	

	toc

	O = zeros(3,3);
	S = zeros(3,3);

	for k =1:Nz
		for i =1:Nx
		for j =1:Ny
		S(1,1) = S_11(j,i,k);%mvg.dudx(i,j,k);
                S(1,2) = S_12(j,i,k);%0.5*( mvg.dudy(i,j,k) +mvg.dvdx(i,j,k) );
                S(1,3) = S_13(j,i,k);%0.5*( mvelgz.dudz(i,j,kstart+k)+mvg.dwdx(i,j,k));
                S(2,1) = S(1,2);
                S(2,2) = S_22(j,i,k);% mvg.dvdy(i,j,k);
                S(2,3) = S_23(j,i,k);%0.5*( mvelgz.dvdz(i,j,kstart+k)+mvg.dwdy(i,j,k));
                S(3,1) = S(1,3);%,ks 0.5*( mvelgz.dudz(i,j,kstart+k)+mvg.dwdx(i,j,k));
                S(3,2) = S(2,3);
                S(3,3) = S_33(j,i,k);%mvelgz.dwdz(i,j,kstart+k);

                O(1,3) = O_13(j,i,k);%0.5*mo.omega_y(i,j,kstart+k);
                O(2,1) = O_21(j,i,k);%0.5*mo.omega_z(i,j,kstart+k);
                O(3,2) = O_32(j,i,k);% 0.5*mo.omega_x(i,j,kstart+k);
                O(1,2) =-O(2,1);
                O(2,3) =-O(3,2);
                O(3,1) =-O(1,3);

                A = S*S + O*O;
		B = O*O';
		C = S*S';
		
                ll = sort(eig(A));
                lambda2(j,i,k) = single(ll(2));
		Q(j,i,k)= 0.5*(trace(B)-trace(C));
		end
		end
	end
	fl=sprintf("lambda_%03d",time);
	ml=matfile(fl,'Writable',true);
	ml.lambda2=single(lambda2);
	ml.Q=single(Q);

	toc
end

