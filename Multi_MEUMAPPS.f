! =========================================================================
!  program to mest solution of time dependent Ginzburg-Landau equation
! =========================================================================
! [MEUMAPPS-SS]
!  program to mest solution of time dependent Ginzburg-Landau equation
!
!  The code has been originally devloped by B. Rad and his colleagues
!  [Radhakrishnan, Gorti, and Babu, Metall. Mater. Trans. A (2016)]
!  B. Rad and Y. Song have updated to use multi-elements and multi-phases
!  Reference for the update: [Zhou, et al., Acta Mater. (2014)]
!
! =========================================================================

      program Multi_MEUMAPPS

! =========================================================================

      !$acc routine (k_space) gang

! =========================================================================

! use p3dfft and mpi libraries
      use p3dfft
! include modules
! ModGparams: include parameters
      use ModGparams
! ModFunc: include subroutines of f_trans, inv_trans, invert, and langevin
      use ModFunc
! ModInitials: include subroutines of eltensor, InitElasticMat, InitTransGrad, k_space, InitFields, Langevin and NucAdd
      use ModInitials
! ModElastic: include subroutines of eldis and elener
      use ModElastic

      implicit none
      include 'mpif.h'

! =========================================================================
! Local parameters are defined here.
! Most parameters are defined in ModGparam.

      integer i, j, k, l
      integer ii, jj, kk, i0, i1, j0
      integer ivar, jvar
      integer nphc, nphc1, nelc, pmax
      integer ivn, fvn
      integer ifreq

      real *8 pphiloc
      real *8 con_min, glob_min, con_max, glob_max, max_phi, phi_max
      real *8 term, term11, termch, tprc
      real *8 t_step, eldis_time, time(6)

      complex(p3dfft_type) term1, term2, term3
      double complex k_sq, k1_sq, k2_sq, k3_sq, k1k2, k2k3, k3k1

! Arrays =================================
! for phase field
      real *8, dimension(:,:,:,:), allocatable :: phi
      real *8, dimension(:,:,:,:), allocatable :: phi_noise

! for concentrations
      real *8, dimension(:), allocatable :: Cloc, Caviter
      real *8, dimension(:,:,:,:), allocatable :: Cons

! for energy related
      real *8, dimension(:), allocatable :: GibbsF

      real *8, dimension(:,:,:), allocatable :: el_en
      real *8, dimension(:,:,:,:), allocatable :: el_grad
      real *8, dimension(:,:,:,:), allocatable :: f_grad, f_diff

      real *8, dimension(:,:,:), allocatable :: u1, u2, u3
      real *8, dimension(:,:,:), allocatable :: grad_1_u1, grad_1_u2, grad_1_u3
      real *8, dimension(:,:,:), allocatable :: grad_2_u1, grad_2_u2, grad_2_u3
      real *8, dimension(:,:,:), allocatable :: grad_3_u1, grad_3_u2, grad_3_u3
      real *8, dimension(:,:,:,:), allocatable :: grad_x, grad_y, grad_z

! for FFT related / complex fields
      real *8, dimension(:,:,:,:), allocatable :: k_del_sq
      real *8, dimension(:,:,:,:), allocatable :: Fconx, Fcony, Fconz

      complex(p3dfft_type), dimension(:,:,:,:), allocatable :: dft_phi
      complex(p3dfft_type), dimension(:,:,:,:), allocatable :: dft_phi_noise

      complex(p3dfft_type), dimension(:,:,:,:), allocatable :: dft_Cons
      complex(p3dfft_type), dimension(:,:,:,:), allocatable :: dft_Fconx, dft_Fcony, dft_Fconz

      complex(p3dfft_type), dimension(:,:,:,:), allocatable :: dft_grad_x, dft_grad_y, dft_grad_z
      complex(p3dfft_type), dimension(:,:,:,:), allocatable :: dft_f_grad, dft_f_diff
      
      complex(p3dfft_type), dimension(:,:,:,:), allocatable :: dft_k_del_sq

! for additional
      real *8, dimension(:,:,:), allocatable :: dummy
      real *8, dimension(:,:,:), allocatable :: dummy2
      real *8, dimension(:,:,:), allocatable :: dummy3

      complex(p3dfft_type), dimension(:,:,:), allocatable :: dft_dummy
      complex(p3dfft_type), dimension(:,:,:), allocatable :: dft_dummy2
      complex(p3dfft_type), dimension(:,:,:), allocatable :: dft_dummy3

! =========================================================================

! =========================================================================
!     MPI Initializations
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
! =========================================================================

! =========================================================================
      initerr = 0
      e_mean = 0.d0

! ===== Read data from input file

      open(2, file='Input.txt', status='old')
      time = 0.d0
      time(1) = MPI_Wtime()
      eldis_time = 0.d0

! MPI sets
      read(2,*) buffchar
      read(2,*) ndim, dims(1), dims(2)
      if ( ndim .eq. 1 ) then
         dims(1) = 1
         dims(2) = nprocs
      end if

! simulation pararmeters
      read(2,*) buffchar
      read(2,*) Nx, Ny, Nz, del_N
      read(2,*) nrun, step, elastep0, Nt, t_step, noise_step
      read(2,*) b_mob, diff_coeff, m_star
      read(2,*) sd, noise
      read(2,*) ifreq

! nucleation
      read(2,*) buffchar
      read(2,*) nucstat, targv, vartot, nmaxppt, iseed
      read(2,*) var_lo, var_hi
      read(2,*) radppt(1), radppt(2), radppt(3)
      read(2,*) addrad(1), addrad(2), addrad(3)
      read(2,*) disc_ener

! total phases/elements
      read(2,*) buffchar
      read(2,*) numph, numel
      read(2,*) intfantiph
      read(2,*) e_mean(1), e_mean(2), e_mean(3)

! ==== allocate additional alloy
      numel = numel - 1
      zdim = numph * numel
      pmax = 2 * numel + 1
      ooNxyz = 1.d0/float(Nx*Ny*Nz)

      allocate ( Nvar(numph), GibbsF(numph), Wellph(numph), SumE(numph) )
      allocate ( Ctype(numph), CEla(6,6,numph), SEla(6,6,numph) )
      allocate ( Cav(numel), Cimax(numel), Cigrad(numel) )
      allocate ( Gparams(pmax, numph), EIntfph(3, numph), TStrph(3, numph) )

      TStrph = 0.d0
      CEla = 0.d0
      Gparams = 0.d0
      EIntfph = 0.d0

! =====
! ===== continue reading the file

! Concentrations:
      read(2,*) buffchar
      read(2,*) constat
      do i = 1, numel
        read(2,*) nelc, Cav(nelc), Cimax(nelc), Cigrad(nelc)
      end do

! phase 1: matrix
      read(2,*) buffchar
      read(2,*) phi_init
      read(2,*) nphc, Nvar(nphc), Ctype(nphc)
      read(2,*) CEla(1,1,nphc),CEla(1,2,nphc),CEla(4,4,nphc)
      read(2,*) (Gparams(j,nphc), j=1, pmax)

! phase 2+: precipitates
      do i = 2, numph
        read(2,*) buffchar
        read(2,*) nphc, Nvar(nphc), Ctype(nphc)

        if(Ctype(nphc).eq.1) then
        read(2,*) CEla(1,1,nphc),CEla(1,2,nphc),CEla(4,4,nphc)
        else
        read(2,*) CEla(1,1,nphc), CEla(1,2,nphc), CEla(1,3,nphc), CEla(2,2,nphc), CEla(2,3,nphc), CEla(3,3,nphc)
        read(2,*) CEla(4,4,nphc), CEla(4,5,nphc), CEla(4,6,nphc), CEla(5,5,nphc), CEla(5,6,nphc), CEla(6,6,nphc)
        end if

        read(2,*) TStrph(1, nphc), TStrph(2, nphc), TStrph(3, nphc)
        read(2,*) EIntfph(1,nphc), EIntfph(2,nphc), EIntfph(3,nphc)
        write(95,*) 'reading Gparams'
        call flush(95)
        read(2,*) (Gparams(j,nphc), j=1, pmax)

        if (i.ne.nphc) initerr=4
      end do

      close(2)

      ! initial check:
      if ((vartot + 1 - sum(Nvar)).ne.0) initerr = 1
      if (nelc.ne.numel) initerr = 5

! =========================================================================

! =========================================================================
!     Default value of free energy from thermocalc is J/mole
!     Convert to J/m^3 using molar volume of the alloy 

      do nphc = 1, numph
        do nelc = 1, numel

            Gparams(nelc, nphc) = Gparams(nelc, nphc) / 1.055e-5

        end do

        Gparams(nelc, nphc) = Gparams(nelc, nphc) / 1.055e-5

      end do

      iseed = iseed - myid

!     Increase the Gibbs Energy of the matrix by an amount equal to the
!     stored energy of deformation due to dislocations

      do nelc = 1, numel
         Gparams(nelc,1) = Gparams(nelc,1) + disc_ener 
      end do

! =========================================================================
!     Set up work structures for P3DFFT
!     -- note use of TRANSPOSED arrays
      call p3dfft_setup (dims, Nx, Ny, Nz, MPI_COMM_WORLD, Nx, Ny, Nz, .true.)

!     Get dimensions for the original array of real numbers, X-pencils
      call p3dfft_get_dims(ist, ien, isize, 1)

!     Get dimensions for the R2C-forward-transformed array of complex numbers
!     Z-pencils
      call p3dfft_get_dims(fst, fen, fsize, 2)
! =========================================================================

! =========================================================================
!     Determine local array sizes and allocate arrays
! Z matrix
      allocate ( Zintf(zdim, zdim), Zinv(zdim, zdim) )
      allocate ( Zr(zdim) )
      allocate ( Xintf(numel, numph) )

! concentrations
      allocate ( Cloc(numel), Caviter(numel), Cmobs(numel), MuC(numel) )
      allocate ( Cons(numel, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( dft_Cons(numel, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )

! phase-field
      allocate ( Pph(numph), Gph(numph) )
      allocate ( Ppr_phi(vartot) )

      allocate ( phi(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( phi_noise(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( P_phi(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )

      allocate ( dft_phi(vartot, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate ( dft_phi_noise(vartot, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      
! energy related
      allocate ( f_grad(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( f_diff(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( el_grad(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      
      allocate ( u1(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( u2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( u3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      
      allocate ( grad_1_u1(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( grad_1_u2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( grad_1_u3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )

      allocate ( grad_2_u1(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( grad_2_u2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( grad_2_u3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )

      allocate ( grad_3_u1(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( grad_3_u2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( grad_3_u3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      
      allocate ( grad_x(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( grad_y(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( grad_z(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )

      allocate ( dft_grad_x(vartot, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate ( dft_grad_y(vartot, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate ( dft_grad_z(vartot, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )

      allocate ( Fconx(numel, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( Fcony(numel, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( Fconz(numel, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
            
      allocate ( dft_Fconx(numel, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate ( dft_Fcony(numel, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate ( dft_Fconz(numel, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      
      allocate ( dft_f_grad(vartot, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate ( dft_f_diff(vartot, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      
      allocate ( el_en(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      
      allocate ( grad_coeff(3,3,vartot),  trans_ppt(6,vartot) )
      
! FFT related
      allocate ( kf(3, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate ( k_del_sq(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      
      allocate ( dft_k_del_sq(vartot, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate ( kf_sq(3, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      
! additional
      allocate ( dummy(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( dummy2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      allocate ( dummy3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3)) )
      
      allocate ( dft_dummy(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate ( dft_dummy2(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )
      allocate ( dft_dummy3(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3)) )

! =========================================================================
! start initialization
! =========================================================================
! file number: to save data
      write(file_num,'(i4.4)') myid+1
! =========================================================================
! =========================================================================
! Complete the elasticity matrix for the phase assuming homogeneous
! elasticity for all the cubic phases by applying cubic symmetry

      call InitElasticMat()

! =========================================================================
      trans_ppt = 0.d0
      grad_coeff = 0.d0

! calculate the transformation strain for all variants
! in the computational frame

      ! interface width
      b_width = 5.0 * del_N

      ! for interfacial energy
      SumE = sum(EIntfph, dim = 1)

      ivn = 0

      do nphc = 2, numph

        ! double well depth
        Wellph(nphc) = 6.0 * SumE(nphc) / b_width
        if(myid.eq.0) write(*,*) '\n Wellph=', nphc, Wellph(nphc)

        ! in the subroutine: determine trans_ppt and grad_coeff
        ivn = ivn + Nvar(nphc-1)

        call InitTransGrad(nphc, ivn)

      end do

      ! Anti-phase: consider only one value - can be updated later
      Wellph(1) = 6.0 * intfantiph / b_width
      if(myid.eq.0) write(*,*) '\n well_var: Wellph(1)=', Wellph(1)

! ===== check calculated values

      if(myid.eq.0) then
        do i = 1,vartot
            write(95,*) 'var=', i
            write(95,*) (trans_ppt(j,i), j=1,6)
        end do

        do i = 1,vartot
            do j = 1, 3
                write(95,*) 'var=', i
                write(95,*) (grad_coeff(j,k,i), k=1,3)
            end do
        end do
      end if
      call flush(95)

! =========================================================================
! Compute (non-dimensionalizing) gradient coefficients for various interfaces

      del_N_sq = del_N*del_N
      grad_coeff = grad_coeff/del_N_sq

      do nelc = 1, numel

        Cmobs(nelc) = diff_coeff / del_N_sq

      end do

      if(myid.eq.0) write(*,*) 'l_star=', l_star
      if(myid.eq.0) write(*,*) 'mob_con1: Cmobs(1) =', Cmobs(1)

!     Non-dimensionalize time_step

      if(myid.eq.0) write(*,*) 't_step=', t_step
      call flush(6)

! =========================================================================
! Set up k vectors in Fourier space
      call k_space()

! =====
      glob_min = 1.d0

! =====
! Set up z matrix here for entries that are constant with time

! initialization
      Zintf = 0.d0

! set up for multi phases/components

      do nelc = 1, numel

        i1 = nelc + numel + 1

        do nphc = 1, numph-1

            ! coordinates
            i0 = (numph - 1) * (nelc - 1) + nphc
            j0 = numph * (nelc - 1) + nphc
        
            nphc1 = nphc + 1

            Zintf(i0, j0)   = Gparams(nelc, nphc)
            Zintf(i0, j0+1) = - Gparams(nelc, nphc1)

            Zr(i0) = Gparams(nelc, nphc) * Gparams(i1, nphc) - Gparams(nelc, nphc1) * Gparams(i1, nphc1)

        end do
      end do
! =========================================================================

! =========================================================================
! Initialization: fields =====

      IF(nrun.eq.1) THEN

        ! Initialize the fields
        ! nucleation: subroutine

        call InitFields(nucstat, constat, phi, Cons, phi_noise)

! ========== end: nucleation cases

      ELSE

! ===== when reading initial fields from a file
        ! check file existance
        INQUIRE(file="data_save.0001", EXIST=file_exists)
        if(.NOT.file_exists) initerr = 3

        open(5,file='data_save.'//file_num, status='old', &
            form='unformatted')

        read(5) ist, ien
        read(5) Cons, phi, el_en
! for the previous type concentration field: below
!        read(5) Cons(1,:,:,:), Cons(2,:,:,:), phi, el_en
        read(5) u1, u2, u3
        read(5) grad_1_u1, grad_1_u2, grad_1_u3
        read(5) grad_2_u1, grad_2_u2, grad_2_u3
        read(5) grad_3_u1, grad_3_u2, grad_3_u3
        close (5)

        if (nucstat.ne.'none') then

            call NucAdd(nucstat, phi)

        end if

      END IF

! end initialization of fields
! =========================================================================

! === end: initial setting
! =========================================================================
! ===== save initial fileds

      write(dir_name,'(a5,i7.7)') 'step_',step-1

      if(myid.eq.0) then
        call system('mkdir -p '//dir_name)
        open(9,file=dir_name//'/data_extract_mult.in')
        write(9,*) Nx, Ny, Nz, nprocs, vartot, numel
        write(9,*) dims
        close(9)
      end if

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      open (4, file=dir_name//'/data_save.'//file_num, &
            form='unformatted')
      rewind(4)
      write(4) ist, ien
      write(4) Cons, phi, el_en
      write(4) u1, u2, u3
      write(4) grad_1_u1, grad_1_u2, grad_1_u3
      write(4) grad_2_u1, grad_2_u2, grad_2_u3
      write(4) grad_3_u1, grad_3_u2, grad_3_u3
      call flush(4)
      close(4)

! =========================================================================

! ===== Error check before proceeding the main loop

      if (initerr.eq.1) write(*,*) "Check variant numbers: Initial error ", initerr
      if (initerr.eq.2) write(*,*) "Check input inp file: Initial error ", initerr
      if (initerr.eq.3) write(*,*) "Check input field files: Initial error ", initerr
      if (initerr.eq.4) write(*,*) "Check variant number: Initial error ", initerr
      if (initerr.eq.5) write(*,*) "Check element number: Initial error ", initerr

      if( initerr .ne. 0) then
        
        STOP "Check the initial error "

      end if

! =====

! =========================================================================
! Align MPI processors before entering the main calculation
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

! =========================================================================

! ===== main loop starts

4     continue

! Introduce fluctuation in phi through the noise term

      if (step.le.noise_step) then

      DO k = ist(3), ien(3)
      DO j = ist(2), ien(2)
      DO i = ist(1), ien(1)

        do ivar = 1, vartot
          call langevin(iseed,sd,noise)
          phi_noise(ivar,i,j,k) = phi_noise(ivar,i,j,k) + noise
        end do

      END DO
      END DO
      END DO

      else

      phi_noise = 0.0

      end if

! ===== variants / precipitates - in Fourier space
      
      DO ivar = 1, vartot

        dummy(:,:,:) = phi(ivar,:,:,:)
        call f_trans(dummy(:,:,:), dft_dummy(:,:,:))


        do k = fst(3), fen(3)
        do j = fst(2), fen(2)
        do i = fst(1), fen(1)

            k1_sq = kf_sq(1,i,j,k)
            k2_sq = kf_sq(2,i,j,k)
            k3_sq = kf_sq(3,i,j,k)
            k1k2 = kf(1,i,j,k)*kf(2,i,j,k)
            k2k3 = kf(2,i,j,k)*kf(3,i,j,k)
            k3k1 = kf(3,i,j,k)*kf(1,i,j,k)

            dft_grad_x(ivar,i,j,k) = dft_dummy(i,j,k) * kf(1,i,j,k)
            dft_grad_y(ivar,i,j,k) = dft_dummy(i,j,k) * kf(2,i,j,k)
            dft_grad_z(ivar,i,j,k) = dft_dummy(i,j,k) * kf(3,i,j,k)

            dft_k_del_sq(ivar,i,j,k) = dft_dummy(i,j,k)* ( grad_coeff(1,1,ivar) * k1_sq &
                                + grad_coeff(2,2,ivar)*k2_sq + grad_coeff(3,3,ivar)*k3_sq &
                                + (grad_coeff(1,2,ivar)+grad_coeff(2,1,ivar)) * k1k2 &
                                + (grad_coeff(2,3,ivar)+grad_coeff(3,2,ivar)) * k2k3 &
                                + (grad_coeff(3,1,ivar)+grad_coeff(1,3,ivar)) * k3k1 )


        end do
        end do
        end do

        call inv_trans( dft_grad_x(ivar,:,:,:), grad_x(ivar,:,:,:) )

        call inv_trans(dft_grad_y(ivar,:,:,:), grad_y(ivar,:,:,:) )

        call inv_trans(dft_grad_z(ivar,:,:,:), grad_z(ivar,:,:,:) )

        call inv_trans(dft_k_del_sq(ivar,:,:,:), k_del_sq(ivar,:,:,:) )

      END DO

! ===== initilize the energy related values at the current iteration

      fch = 0.d0
      fw = 0.d0
      f_int = 0.0
      f_grad = 0.d0

      Fconx = 0.d0
      Fcony = 0.d0
      Fconz = 0.d0

! =====

! calculate the interpolating function p(phi)
! and grain boundary energy function g(phi)
! for all variants based on local phi

      do k = ist(3), ien(3)
      do j = ist(2), ien(2)
      do i = ist(1), ien(1)
    
        ! initialization
        Pph = 0.d0
        Gph = 0.d0

        ! variants
        ivn = 0
        do nphc = 2, numph
            ivn = ivn + Nvar(nphc - 1)
            fvn = ivn + Nvar(nphc) -1

            do ivar = ivn, fvn

                fi = phi(ivar,i,j,k)
                fi_sq = fi * fi

                pphiloc = fi_sq * fi * ( 6.d0*fi_sq - 15.d0*fi + 10.d0 )
                P_phi(ivar,i,j,k) = pphiloc

                Ppr_phi(ivar) = 30.d0*fi_sq * (fi_sq - 2.d0 * fi + 1.d0)

                Pph(1) = Pph(1) + pphiloc
                Pph(nphc) = Pph(nphc) + pphiloc
                Gph(nphc) = Gph(nphc) + fi_sq*(1.0-fi)**2

                DO jvar = 1, vartot
                    if(jvar.ne.ivar) Gph(1) = Gph(1) + fi_sq*phi(jvar,i,j,k)**2
                END DO

            end do

        end do

        Pph(1) = 1.d0 - Pph(1)

! =====
! Calculate the interface concentrations of the components
! The matrix dimension: zdim x zdim
!                       zdim = numph x numel

        do nelc = 1,numel

            i0 = zdim - numel + nelc
            Zr(i0) = Cons(nelc, i,j,k)

            do nphc = 1, numph

                j0 = numph * (nelc - 1) + nphc
                Zintf(i0, j0) = Pph(nphc)

            end do
        end do

        err_flag = 2
        Zinv = Zintf
        call invert (Zinv, zdim, err_flag)

! ===== initialization before update
        Xintf = 0.d0

        do nphc = 1, numph
        do nelc = 1, numel

            ! coordinates
            i0 = numph * (nelc - 1) + nphc
        
            do j0 = 1, zdim

                Xintf(nelc, nphc) = Xintf(nelc, nphc) + Zinv(i0,j0) * Zr(j0)

            end do

        end do
        end do

! ===== interface concentrations are updated

! =====
!     Calculate thermodynamic driving force and their derivatives wrt C

! =====
! Gibbs free energy
        GibbsF = 0.d0

        do nphc = 1, numph

            do nelc = 1, numel

                i1 = nelc + numel + 1

                GibbsF(nphc) = GibbsF(nphc) + Gparams(nelc,nphc) * ( Xintf(nelc, nphc) - Gparams(i1,nphc) )**2

            end do

            GibbsF(nphc) = GibbsF(nphc) + Gparams(nelc,nphc)

        end do

! ===== Gibbs free energy has been updated

! =====
! Calculate the chemical potentials for components to be used
! in the Cahn-Hilliard equation

        do nelc = 1, numel

            i1 = nelc + numel + 1
            MuC(nelc) = 2.0 * Gparams(nelc,1) * ( Xintf(nelc,1) - Gparams(i1,1) )

        end do

! =====
        ! at this moment the first value is saved
        ! remaining summation will do in the following loop
        fchloc = Pph(1) * GibbsF(1)
        fwloc = Wellph(1) * Gph(1)

! =====
! calculate quantities required to solve the solute diffusion equation
! per [Kim, Kim and Suzuki, Pys. Rev. E (1999)]

!       variants
        ivn = 0
        do nphc = 2, numph
            ivn = ivn + Nvar(nphc - 1)
            fvn = ivn + Nvar(nphc) -1

            do ivar = ivn, fvn

! Calculate the energy due to energy well

                term = ( grad_coeff(1,1,ivar)*grad_x(ivar,i,j,k)**2 &
                    + grad_coeff(2,2,ivar)*grad_y(ivar,i,j,k)**2 &
                    + grad_coeff(3,3,ivar)*grad_z(ivar,i,j,k)**2 &
                    + (grad_coeff(1,2,ivar)+grad_coeff(2,1,ivar))*grad_x(ivar,i,j,k)*grad_y(ivar,i,j,k) &
                    + (grad_coeff(1,3,ivar)+grad_coeff(3,1,ivar))*grad_x(ivar,i,j,k)*grad_z(ivar,i,j,k) &
                    + (grad_coeff(2,3,ivar)+grad_coeff(3,2,ivar))*grad_y(ivar,i,j,k)*grad_z(ivar,i,j,k))

                f_int = f_int + 0.5*term

! Calculate derivative of free energy wrt to phi for use in G-L equation
!
! Calculate derivative of system energy with respect to phi
! It is the sum of the bulk energy and the gradient energy

                fi = phi(ivar,i,j,k)
                termch = GibbsF(nphc) - GibbsF(1)

                do nelc = 1, numel

                    xdiffph = ( Xintf(nelc,1) - Xintf(nelc,nphc) )
                    termch = termch + MuC(nelc) * xdiffph

                    tprc = Ppr_phi(ivar) * xdiffph
                    Fconx(nelc,i,j,k) = Fconx(nelc,i,j,k) + tprc * grad_x(ivar,i,j,k)
                    Fcony(nelc,i,j,k) = Fcony(nelc,i,j,k) + tprc * grad_y(ivar,i,j,k)
                    Fconz(nelc,i,j,k) = Fconz(nelc,i,j,k) + tprc * grad_z(ivar,i,j,k)

                end do

                well_grad = 2.0 * Wellph(nphc) * fi * (1.0-fi) * (1.0-2.d0*fi)

                do jvar = 1, vartot
                    if(jvar.ne.ivar) well_grad = well_grad + 4.0 * Wellph(1) * fi * phi(jvar,i,j,k)**2
                end do

                f_grad(ivar,i,j,k) = termch * Ppr_phi(ivar) + well_grad

            end do

            fchloc = fchloc + Pph(nphc) * GibbsF(nphc)
            fwloc = fwloc + Wellph(nphc) * Gph(nphc)

        end do

! =====

        fch = fch + fchloc
! Non-dimensionalize energy wells
        fw = fw + fwloc

! =====

! =====
!      Prevent sum of phi's at any i,j,k exceeding 1.0 when three or four variants
!!     exist together at that i,j,k

        do jvar = 1, vartot

            term11 = 0.d0
            n_tilda = 1

            if(phi(jvar,i,j,k).gt.0.0) then
                do ivar = 1, vartot
                    if(phi(ivar,i,j,k).gt.0.0) then
                        if(ivar.ne.jvar) then
                            term11 =  term11 + f_grad(ivar,i,j,k) - k_del_sq(ivar,i,j,k)
                            n_tilda = n_tilda + 1
                        end if
                    end if
                end do
            end if

            f_diff(jvar,i,j,k) = f_grad(jvar,i,j,k) - term11/n_tilda
        end do

      end do
      end do
      end do
! =========================================================================

! =========================================================================
! elastic energy calculation
      
      if(step.ge.elastep0) then
        el_en = 0.0
        el_grad = 0.0

!        if(el_visit.eq.0) then
!            step_1 = elastep0
!        else
!            step_1 = step
!        end if

        time(2) = MPI_Wtime()

! =====
! Calculate eleastic energy due to the presence of alpha variants
        call eldis(phi, el_en, el_grad, u1, u2, u3, &
                    grad_1_u1, grad_1_u2, grad_1_u3, &
                    grad_2_u1, grad_2_u2, grad_2_u3, &
                    grad_3_u1, grad_3_u2, grad_3_u3)

        time(2) = MPI_Wtime() - time(2)
        eldis_time = eldis_time + time(2)

      else

        el_grad = 0.d0
        el_en = 0.d0

      end if

      f_diff = f_diff + el_grad
    
! =========================================================================


      call MPI_Reduce(f_int, fint_tot, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(fw, fw_tot, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(fch, fch_tot, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(solid_frac, solid_tot, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      solid_tot = solid_tot / float(Nx*Ny*Nz)


      el_sum = sum(el_en)

      call MPI_Reduce(el_sum, elen_tot, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      if ( myid .eq. 0 )  then

        e_tot = fint_tot + fch_tot + elen_tot + fw_tot

        write(98,10) step, fint_tot, fw_tot, fch_tot, &
                elen_tot, e_tot, glob_min, glob_max, phi_max, solid_tot
      end if
      call flush(98)
      
      do nelc = 1, numel

        call f_trans(Cons(nelc,:,:,:), dft_Cons(nelc,:,:,:) )

        call f_trans(Fconx(nelc,:,:,:), dft_Fconx(nelc,:,:,:) )
        call f_trans(Fcony(nelc,:,:,:), dft_Fcony(nelc,:,:,:) )
        call f_trans(Fconz(nelc,:,:,:), dft_Fconz(nelc,:,:,:) )

      end do


! =========================================================================
! ====== update fields

      do k = fst(3), fen(3)
      do j = fst(2), fen(2)
      do i = fst(1), fen(1)

        k_sq = kf_sq(1,i,j,k)+kf_sq(2,i,j,k)+kf_sq(3,i,j,k)

        do nelc = 1, numel

            term1 = kf(1,i,j,k) * dft_Fconx(nelc,i,j,k) &
                    + kf(2,i,j,k) * dft_Fcony(nelc,i,j,k) &
                    + kf(3,i,j,k) * dft_Fconz(nelc,i,j,k)

            term2 = 1.d0 - t_step * Cmobs(nelc) * k_sq
            dft_Cons(nelc, i,j,k) = ( dft_Cons(nelc, i,j,k) + t_step * Cmobs(nelc) * term1 )/term2

        end do

      end do
      end do
      end do

      do nelc = 1, numel

        call inv_trans(dft_Cons(nelc,:,:,:), Cons(nelc,:,:,:) )

      end do


! =========================================================================
! Solve time-dependent G-L equation

      DO ivar = 1, vartot

        dummy(:,:,:) = phi(ivar,:,:,:)
        call f_trans(dummy(:,:,:), dft_dummy(:,:,:) )

        dummy2(:,:,:) = f_diff(ivar,:,:,:)
        call f_trans(dummy2(:,:,:), dft_dummy2(:,:,:) )

        dummy3(:,:,:) = phi_noise(ivar,:,:,:)
        call f_trans(dummy3(:,:,:), dft_dummy3(:,:,:) )

        do k = fst(3), fen(3)
        do j = fst(2), fen(2)
        do i = fst(1), fen(1)

            k1_sq = kf_sq(1,i,j,k)
            k2_sq = kf_sq(2,i,j,k)
            k3_sq = kf_sq(3,i,j,k)
            k1k2 = kf(1,i,j,k)*kf(2,i,j,k)
            k2k3 = kf(2,i,j,k)*kf(3,i,j,k)
            k3k1 = kf(3,i,j,k)*kf(1,i,j,k)

            term1 = t_step * b_mob * dft_dummy2(i,j,k)

            term2 = 1.d0 - t_step * b_mob * (grad_coeff(1,1,ivar) * k1_sq &
                    + grad_coeff(2,2,ivar)*k2_sq + grad_coeff(3,3,ivar)*k3_sq &
                    + (grad_coeff(1,2,ivar)+grad_coeff(2,1,ivar)) * k1k2 &
                    + (grad_coeff(2,3,ivar)+grad_coeff(3,2,ivar)) * k2k3 &
                    + (grad_coeff(3,1,ivar)+grad_coeff(1,3,ivar)) * k3k1 )

            term3 = dft_dummy3(i,j,k)*t_step

            dft_dummy(i,j,k) = ( dft_dummy(i,j,k) - term1 + term3 ) / term2

        end do
        end do
        end do

!     perfrom inverse transformation to get the new phi values
!     This new phi which will be used in the Cahn-Hilliard equation
!     in the next step

        call inv_trans(dft_dummy(:,:,:), dummy(:,:,:) )

        phi(ivar,:,:,:) = dummy(:,:,:)

      END DO 
! =========================================================================

! =========================================================================
! ===== check values


      con_min = 1.0 
      con_max = 0.0
      max_phi = 0.0

      Cloc = 0.d0
      Caviter = 0.d0
      solid_frac = 0.0

      do k = ist(3), ien(3)
      do j = ist(2), ien(2)
      do i = ist(1), fen(1)

        ! about components
        do nelc = 1, numel

            if(Cons(nelc,i,j,k).lt.con_min) con_min = Cons(nelc,i,j,k)
            if(Cons(nelc,i,j,k).gt.con_max) con_max = Cons(nelc,i,j,k)

            Cloc(nelc) = Cloc(nelc) + Cons(nelc,i,j,k)

        end do

        ! about phases
        ! also update solid fraction
        do ivar = 1, vartot

            if(phi(ivar,i,j,k).lt.0.0) phi(ivar,i,j,k)=0.0
            if(phi(ivar,i,j,k).gt.1.0) phi(ivar,i,j,k)=1.0
            if(phi(ivar,i,j,k).gt.max_phi) max_phi= phi(ivar,i,j,k)
            if(phi(ivar,i,j,k).gt.0.0) solid_frac = solid_frac + phi(ivar,i,j,k)

        end do

      end do
      end do
      end do

! =====
! calculate average concentrations for mass conservation

      do nelc = 1, numel

        call MPI_Allreduce(Cloc(nelc), Caviter(nelc), 1, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, MPI_COMM_WORLD, ierr)

        Caviter(nelc) = Caviter(nelc)*ooNxyz

      end do

      call MPI_Allreduce(con_min, glob_min, 1, MPI_DOUBLE_PRECISION, &
           MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(con_max, glob_max, 1, MPI_DOUBLE_PRECISION, &
           MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(max_phi, phi_max, 1, MPI_DOUBLE_PRECISION, &
           MPI_MAX, MPI_COMM_WORLD, ierr)

      if(myid.eq.0) write(99,*) step, (Caviter(i), i = 1,numel)
      call flush(99)

! =====
! increase a time step
      step = step + 1 

! =====
!     Data output for post-processing and restart

      if(mod(step,ifreq).eq.0) then

        write(dir_name,'(a5,i7.7)') 'step_',step

        if(myid.eq.0) then
            call system('mkdir -p '//dir_name)
            open(9,file=dir_name//'/data_extract_mult.in')
            write(9,*) Nx, Ny, Nz, nprocs, vartot, numel
            write(9,*) dims
            close(9)
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        open (4, file=dir_name//'/data_save.'//file_num, &
            form='unformatted')
        rewind(4)
        write(4) ist, ien
        write(4) Cons, phi, el_en
        write(4) u1, u2, u3
        write(4) grad_1_u1, grad_1_u2, grad_1_u3
        write(4) grad_2_u1, grad_2_u2, grad_2_u3
        write(4) grad_3_u1, grad_3_u2, grad_3_u3
        call flush(4)

        close (4)

      end if

! align MPI processors
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

! back to the initial part of the loop
      if(step.lt.Nt) go to 4
! end main loop
! =========================================================================

! output formats ================================================
3     format(3i6,3e15.7)
9     format(8e15.7)
10    format(i8, 10e15.7)
11    format(3i8,2e15.7)

! =========================================================================
! ===== deallocate: free memory

      deallocate (Cons, dft_Cons)
      deallocate (phi, phi_noise, P_phi, Ppr_phi, dft_phi, dft_phi_noise)
      deallocate (f_grad, f_diff, el_grad)
      deallocate (u1, u2, u3)
      deallocate (grad_1_u1, grad_1_u2, grad_1_u3)
      deallocate (grad_2_u1, grad_2_u2, grad_2_u3)
      deallocate (grad_3_u1, grad_3_u2, grad_3_u3)
      deallocate (grad_x, grad_y, grad_z)
      deallocate (dft_grad_x, dft_grad_y, dft_grad_z)
      deallocate (Fconx, Fcony, Fconz)
      deallocate (dft_Fconx, dft_Fcony, dft_Fconz)
      deallocate (dft_f_grad, dft_f_diff, dft_k_del_sq)
      deallocate (el_en, grad_coeff)
      deallocate (kf, k_del_sq, kf_sq)
      deallocate (dummy, dummy2, dummy3, dft_dummy, dft_dummy2, dft_dummy3)
      deallocate (Nvar)
      deallocate (Zr, Zinv, Zintf, Xintf)
      deallocate (Cloc, Caviter, Cmobs, MuC, Cav, Cimax, Cigrad)
      deallocate (Pph, Gph, Ctype, CEla, SEla)
      deallocate (GibbsF, Gparams, EIntfph, SumE, Wellph)
      deallocate (trans_ppt, TStrph)

! =========================================================================
      call p3dfft_clean

      call MPI_Finalize(ierr)
! =========================================================================

! =========================================================================
! end program
      end
! =========================================================================



    
