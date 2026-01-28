Nx=640;
Ny=512;
Nz=192;
kcond=64;

	fn=sprintf('../data/fieldsq2neg_lse_k_%03d.mat',kcond);
	%fn=sprintf('../data/fieldsq4_lse_k_%03d.mat',kcond);
	mn=matfile(fn,"Writable",true);
	S_11 = single(zeros(Nx,Ny,Nz));
	S_12 = single(zeros(Nx,Ny,Nz));
	S_13 = single(zeros(Nx,Ny,Nz));
	S_22 = single(zeros(Nx,Ny,Nz));
	S_23 = single(zeros(Nx,Ny,Nz));
	S_33 = single(zeros(Nx,Ny,Nz));	
	O_21 = single(zeros(Nx,Ny,Nz));
	O_13 = single(zeros(Nx,Ny,Nz));
	O_32 = single(zeros(Nx,Ny,Nz));
	lambda2 = single(zeros(Nx,Ny,Nz));
	Q = single(zeros(Nx,Ny,Nz));
	
	tic
	S_11=		mn.lse(1:end,1:end,1:end,6);
	S_12=0.5*( mn.lse(1:end,1:end,1:end,9) + mn.lse(1:end,1:end,1:end,7) );
	S_13=0.5*( mn.lse(1:end,1:end,1:end,8) + mn.lse(1:end,1:end,1:end,12) );
	S_22=		mn.lse(1:end,1:end,1:end,10);
	S_23=0.5*( mn.lse(1:end,1:end,1:end,11)+mn.lse(1:end,1:end,1:end,13) );
	S_33=mn.lse(1:end,1:end,1:end,14);
 		
	O_21 = 0.5*(mn.lse(1:end,1:end,1:end,7)  - mn.lse(1:end,1:end,1:end,9)) ;
	O_13 = 0.5*(mn.lse(1:end,1:end,1:end,12) - mn.lse(1:end,1:end,1:end,8))  ;
	O_32 = 0.5*(mn.lse(1:end,1:end,1:end,11) - mn.lse(1:end,1:end,1:end,13))  ;	

	toc

	O = zeros(3,3);
	S = zeros(3,3);

	for k =1:Nz
		for i =1:Nx
		for j =1:Ny
		S(1,1) = S_11(i,j,k);%mvg.dudx(i,j,k);
                S(1,2) = S_12(i,j,k);%0.5*( mvg.dudy(i,j,k) +mvg.dvdx(i,j,k) );
                S(1,3) = S_13(i,j,k);%0.5*( mvelgz.dudz(i,j,kstart+k)+mvg.dwdx(i,j,k));
                S(2,1) = S(1,2);
                S(2,2) = S_22(i,j,k);% mvg.dvdy(i,j,k);
                S(2,3) = S_23(i,j,k);%0.5*( mvelgz.dvdz(i,j,kstart+k)+mvg.dwdy(i,j,k));
                S(3,1) = S(1,3);%,ks 0.5*( mvelgz.dudz(i,j,kstart+k)+mvg.dwdx(i,j,k));
                S(3,2) = S(2,3);
                S(3,3) = S_33(i,j,k);%mvelgz.dwdz(i,j,kstart+k);

                O(1,3) = O_13(i,j,k);%0.5*mo.omega_y(i,j,kstart+k);
                O(2,1) = O_21(i,j,k);%0.5*mo.omega_z(i,j,kstart+k);
                O(3,2) = O_32(i,j,k);% 0.5*mo.omega_x(i,j,kstart+k);
                O(1,2) =-O(2,1);
                O(2,3) =-O(3,2);
                O(3,1) =-O(1,3);

                A = S*S + O*O;
		B = O*O';
		C = S*S';
		
                ll = sort(eig(A));
                lambda2(i,j,k) = single(ll(2));
		Q(i,j,k)= 0.5*(trace(B)-trace(C));
		end
		end
	end
	mn.lambda2=single(lambda2);
	mn.Q=single(Q);

	toc

