      module ModInitials

      use ModFunc
      use ModGparams

      implicit none
! =========================================================================

! =========================================================================
      contains
! Contains: eltensor, InitElasticMat, InitTransGrad, k_space,
!           InitFields, NucAdd

! =========================================================================
! ===== eltensor
        subroutine eltensor(C, cryst_type)

        implicit none
        real *8 C(6,6), temp(6,6), sum1
        integer cryst_type
        integer i, j, m, n, p, q, s, t, ij, mn, pq, st

        ! Cubic crystals - C_11, C_12, C_44 are given

        if ( cryst_type .eq. 1 ) then
            C(1,3) = C(1,2)
            C(2,2) = C(1,1)
            C(2,3) = C(1,3)
            C(3,3) = C(1,1)
            C(5,5) = C(4,4)
            C(6,6) = C(4,4)
        end if

        ! Tetragonal crystala - C_11, C_12, C_13, C_33, C_44 and C_66 are
        ! given

        if ( cryst_type.eq. 2) then
            C(2,2) = C(1,1)
            C(5,5) = C(4,4)
            C(2,3) = C(1,3)
        end if

        ! Hexagonal crystals - C_11, C_12, C_13, C_33, C_44 are given

        if ( cryst_type .eq. 3 ) then
            C(2,2) = C(1,1)
            C(2,3) = C(1,3)
            C(5,5) = C(4,4)
            C(6,6) = 0.5d0*(C(1,1)-C(1,2))
        end if

        ! Trigonal crystals - C_11, C_12, C_13, C_33, C_44, C_15 are given

        if ( cryst_type .eq. 4 ) then
            C(2,2) = C(1,1)
            C(2,3) = C(1,3)
            C(2,5) = -C(1,5)
            C(4,6) = -C(1,5)
            C(5,5) = C(4,4)
            C(6,6) = 0.5d0*(C(1,1)-C(1,2))
        end if

        ! Fill in lower part of C matrix by symmetry

        do j=2,6
        do i=1,j-1
            C(j,i) = C(i,j)
        end do
        end do

        return
        end subroutine eltensor
! =========================================================================

! ===== InitElasticMat
        subroutine InitElasticMat()

        integer i, j
        real *8 Ctemp(6,6)

        do i = 1, numph
            Ctemp = CEla(:,:,i)

            if (Ctemp(1,1) .le. 0 ) then
                Ctemp = 0.d0
                CEla(:,:,i) = Ctemp
                SEla(:,:,i) = Ctemp

            else

                ! Complete the elasticity matrix
                call eltensor(Ctemp, Ctype(i))
            
                ! update calculated results
                CEla(:,:,i) = Ctemp(:,:)

                ! Invert the elasticty matrix for the precipitate phase
                err_flag = 20

                call invert(Ctemp, 6, err_flag)
            
                ! save results
                SEla(:,:,i) = Ctemp
            end if
        end do

        ! check values
        do j = 1,6
            do i = 1,6
                write(93,*) i, j, CEla(i,j,1), CEla(i,j,2)
            end do
        end do
        call flush(93)

        return
        end subroutine InitElasticMat
! =====
! =========================================================================


! ===== InitTransGrad
        subroutine InitTransGrad(nphc, i0)

! Calculate Transformation strain and gradient coefficients
! for all variatns

        integer nphc, i0, nv
        integer fvn, ivar, i, j

        real *8 TStr(3), EIntf(3)
        real *8 Q(3,3), Q_trans(3,3)
        real *8 A(3,3), A_pr(3,3)
        real *8 Temp(3,3), R(3,3)
        real *8 B(3,3), B_pr(3,3)
        real *8 SYM(3,3)

        logical :: file_exists

! =====

        nv = Nvar(nphc)
        EIntf = EIntfph(:,nphc)
        TStr = TStrph(:,nphc)

        ! final variant number
        fvn = i0 + nv - 1

        if (myid.eq.0) then

            write(*,*) "\n Phase number: ", nphc
            write(*,*) "Variant number from ", i0
            write(*,*) "               to   ", fvn

        end if

! =====
        if (nv.eq.1) then

            trans_ppt(1,i0) = TStr(1)
            trans_ppt(2,i0) = TStr(2)
            trans_ppt(3,i0) = TStr(3)

            grad_coeff(1,1,i0) = b_width * EIntf(1)
            grad_coeff(2,2,i0) = b_width * EIntf(2)
            grad_coeff(3,3,i0) = b_width * EIntf(3)

        end if

! =====

        if (nv.eq.3) then

            ivar = i0
            trans_ppt(1, ivar) = TStr(1)
            trans_ppt(2, ivar) = TStr(2)
            trans_ppt(3, ivar) = TStr(3)

            ivar = i0 + 1
            trans_ppt(1, ivar) = TStr(2)
            trans_ppt(2, ivar) = TStr(1)
            trans_ppt(3, ivar) = TStr(3)

            ivar = i0 + 2
            trans_ppt(1, ivar) = TStr(2)
            trans_ppt(2, ivar) = TStr(3)
            trans_ppt(3, ivar) = TStr(1)

            do ivar = i0, fvn

                grad_coeff(1,1,ivar) = b_width * EIntf(1)
                grad_coeff(2,2,ivar) = b_width * EIntf(2)
                grad_coeff(3,3,ivar) = b_width * EIntf(3)

            end do


        end if
! =====

        if (nv.eq.12) then

            INQUIRE(FILE="sym_var.inp", EXIST=file_exists)

            if (file_exists) then

                open(21, file='sym_var.inp', status='old')

            else

                initerr = 2
                if(myid.eq.0) write(*,*) "\n Check file: inp file does not exist. \n "
                return

            end if

! Rotation matrix to go from I:[1 0 0],[0 1 0],[0 0 1] to delta
! habit plane taken from Zhou paper

            Q(1,1) =  1.d0/dsqrt(2.0)
            Q(2,1) =  0.d0
            Q(3,1) = -1.d0/dsqrt(2.0)

            Q(1,2) =  1.d0/dsqrt(3.0)
            Q(2,2) =  1.d0/dsqrt(3.0)
            Q(3,2) =  1.d0/dsqrt(3.0)

            Q(1,3) =  1.d0/dsqrt(6.0)
            Q(2,3) = -2.d0/dsqrt(6.0)
            Q(3,3) =  1.d0/dsqrt(6.0)

            if(myid.eq.0) then
                write(*,*) 'Q'

                do i = 1,3
                    write(*,*) (Q(i,j), j=1,3)
                end do

            end if

            Q_trans = transpose(Q)
            if(myid.eq.0) then
                write(*,*) 'Q_trans'

                do i = 1,3
                    write(*,*) (q_trans(i,j), j=1,3)
                end do
            end if

! Create gradient coefficient matrix in local frame K
! =====
            A = 0.d0

            A(1,1) = b_width * EIntf(1)
            A(2,2) = b_width * EIntf(2)
            A(3,3) = b_width * EIntf(3)

            if(myid.eq.0) then
                write(*,*) 'A'

                do i = 1,3
                    write(*,*) (a(i,j), j=1,3)
                end do
            end if

! Rotate gradient coefficient matrix from local frame K to
! reference frame I using [Q][A][Q^T]
! =====
            Temp = matmul(Q, A)

            if(myid.eq.0) then
                write(*,*) 'Temp'

                do i = 1,3
                    write(*,*) (Temp(i,j), j=1,3)
                end do
            end if

! =====
            R = matmul(Temp, Q_trans)

            if(myid.eq.0) then
                write(*,*) 'R'

                do i = 1,3
                    write(*,*) (R(i,j), j=1,3)
                end do
            end if

! =====
            A = R

! =====
! Create transformation strain matrix in local frame J
            B = 0.d0

            B(1,1) = TStr(1)
            B(2,2) = TStr(2)
            B(3,3) = TStr(3)

! Rotate transformation strain matrix from local frame J to
! reference frame I using [Q][B][Q^T]
! =====
            Temp = matmul(Q,B)
            B = matmul(Temp,transpose(Q))

! Loop over variants
            do ivar = i0, fvn
      
! Read in the rotation matrix for symmetry operation of that variant
                do i=1,3
                    read(21,*) (SYM(i,j), j=1,3)
                end do

! Rotate gradient coefficient matrix from reference frame I to
! computational frame J using [S^T][A][S]
! =====
                Temp = matmul(transpose(SYM), A)

                A_pr = matmul(Temp, SYM)
                grad_coeff(:,:,ivar) = A_pr(:,:)

! Rotate transformation strain matrix from reference frame I to
! computational frame J using [S^T][B][S]
! =====
                Temp = matmul(transpose(SYM), B)
                B_pr = matmul(Temp, SYM)

                trans_ppt(1,ivar) = B_pr(1,1)
                trans_ppt(2,ivar) = B_pr(2,2)
                trans_ppt(3,ivar) = B_pr(3,3)

                trans_ppt(4,ivar) = B_pr(2,3)
                trans_ppt(5,ivar) = B_pr(1,3)
                trans_ppt(6,ivar) = B_pr(1,2)

            end do

            close(21)
            
        end if

        return
        end subroutine InitTransGrad

! =========================================================================

! ===== k_space
        subroutine k_space()

        integer i, j, k, ii, jj, kk
        real *8 pi
        double complex lamda1, lamda2, lamda3
        double complex unit, unit1, unit2, unit3

        pi = datan(1.d0)*4.0
        unit = (0.0,1.0)
        unit = 2.d0*pi*unit
        unit1 = unit/float(Nx)
        unit2 = unit/float(Ny)
        unit3 = unit/float(Nz)

!       !$acc parallel loop collapse (3) gang
        do i = fst(3), fen(3)
            do j = fst(2), fen(2)
                do k = fst(1), fen(1)

            ii = i-1

            if(ii.le.Nx/2) then
                lamda1 = float(ii) * unit1
            else
                lamda1 = float(ii-Nx) * unit1
            end if

                jj = j-1

                if(jj.le.Ny/2) then
                    lamda2 = float(jj) * unit2
                else
                    lamda2 = float(jj-Ny) * unit2
                end if

                    kk = k-1

                    if(kk.le.Nz/2) then
                        lamda3 = float(kk) * unit3
                    else
                        lamda3 = float(kk-Nz) * unit3
                    end if

                    kf(1,k,j,i) = lamda1
                    kf(2,k,j,i) = lamda2
                    kf(3,k,j,i) = lamda3

                    if (2*ii.eq.Nx) kf(1,k,j,i) = 0.d0
                    if (2*jj.eq.Ny) kf(2,k,j,i) = 0.d0
                    if (2*kk.eq.Nz) kf(3,k,j,i) = 0.d0

                    kf_sq(1,k,j,i) = lamda1*lamda1
                    kf_sq(2,k,j,i) = lamda2*lamda2
                    kf_sq(3,k,j,i) = lamda3*lamda3

                end do
            end do
        end do

        return

        end subroutine k_space

! =========================================================================

! ===== InitFields

        subroutine InitFields(nucstatus, constat, phi, Cons, phinoise)

        include 'mpif.h'

        ! status
        character nucstatus*10
        character constat*16

        integer r0ppt, rdist
        integer i, j, k, ii, jj, kk, ix, jy, kz, ivar
        integer mnx, mny, mnz
        integer pcount, nppt, pcount_tot
        integer randvar
        real *8 term1, term2, term3, sum_ijk
        real *8 cdifflim

        ! random number
        real *8 ran_2

        ! field informations
        real *8 phi(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 Cons(numel, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 phinoise(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))

        ! ===== Main Subroutine

        ! ===== set useful parameters
        mnx = Nx/2
        mny = Ny/2
        mnz = Nz/2

        cdifflim = 1.e-9

        r0ppt = radppt(1)**2 + radppt(2)**2 + radppt(3)**2

        ! ===== Initialization: for all fields
        phi = phi_init
        phinoise = 0.d0

        do i = 1, numel

            Cons(i,:,:,:) = Cav(i)

        end do
      
! ===== Single nucleus case
! a single nucleus is present inside the matrix
        if(nucstatus.eq.'single') then

            if (targv.eq.0) targv = 1
            if(myid.eq.0) write(*,*) "Nucleus information: single variant ", targv

            do k = -radppt(3), radppt(3)
            do j = -radppt(2), radppt(2)
            do i = -radppt(1), radppt(1)

                term1 = (float(i)/float(radppt(1)))**2
                term2 = (float(j)/float(radppt(2)))**2
                term3 = (float(k)/float(radppt(3)))**2
                sum_ijk = term1 + term2 + term3

                if(sum_ijk.le.1.0) then
                    ii = i + mnx
                    jj = j + mny
                    kk = k + mnz

                    if(ii.ge.ist(1).AND.ii.le.ien(1).AND. &
                    jj.ge.ist(2).AND.jj.le.ien(2).AND. &
                    kk.ge.ist(3).AND.kk.le.ien(3)) then

                        phi(targv,ii,jj,kk) = 1.0

                    end if
                end if

            end do
            end do
            end do


            ! ===== Temporal: update later
            if ((constat.eq.'cppt')) then
                if (myid .eq. 0) write(*,*) "Nucleus: different initial ppt size"

                do k  = ist(3), ien(3)
                do j  = ist(2), ien(2)
                do i  = ist(1), ien(1)

                    rdist = (i - mnx)**2 + (j - mnx)**2 + (k - mnx)**2

                    if(rdist.le.r0ppt) then

                        phi(targv,i,j,k) = 1.0

                    end if

                end do
                end do
                end do

            end if
            ! =====

        end if

! ===== end - case: single

! ===== Multiple nucleation case =====
! Initially three nuclei of different variants placed in the domain
! Subsequent nucleation is based on the composite nucleus model

        if(nucstatus.eq.'multiple') then
            nppt = 1

            if(targv.eq.0) then
                if(myid.eq.0) write(*,*) "Nucleus status: Multiple"
            else
                if(myid.eq.0) write(*,*) "Nucleus status: Multiple with a variant", targv
            endif

            ! set a center position of a precipitate
            ! r0ppt = radppt(1)**2 + radppt(2)**2 + radppt(3)**2

13          if(myid.eq.0) then

16              ix = ran_2(iseed)*Nx + 1.0
                jy = ran_2(iseed)*Ny + 1.0
                kz = ran_2(iseed)*Nz + 1.0

17              if(targv.eq.0) then
                    ! at this moment do not consider random nucleation of laves phase
                    randvar = var_lo + ran_2(iseed)*(var_hi-var_lo) + 1.0
                else
                    randvar = targv
                endif

            end if

            call MPI_Bcast(ix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(jy, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(kz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(randvar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

! ====== check whether precipitates overlap the others
            pcount = 0

            !=====
            ! r0ppt = (radppt(1)+addrad)**2 + (radppt(2)+addrad)**2 + (radppt(3)+addrad)**2
            ! r0ppt = (radppt(1))**2 + (radppt(2))**2 + (radppt(3))**2
            !=====

            do k  = ist(3), ien(3)
            do j  = ist(2), ien(2)
            do i  = ist(1), ien(1)

                rdist = (i - ix)**2 + (j - jy)**2 + (k-kz)**2

                do ivar = 1, vartot
                    if((rdist.le.r0ppt).AND.(phi(ivar,i,j,k).gt.0.0)) then
                        pcount = pcount +1
                    end if
                end do

            end do
            end do
            end do

            call MPI_Allreduce(pcount, pcount_tot, 1, MPI_INTEGER, &
                            MPI_SUM, MPI_COMM_WORLD, ierr)

            if(pcount_tot.gt.0) go to 16

! ===== change the precipitate

            if(pcount_tot.eq.0) then


                do k  = ist(3), ien(3)
                do j  = ist(2), ien(2)
                do i  = ist(1), ien(1)


                    rdist = (i - ix)**2 + (j - jy)**2 + (k-kz)**2

                    if(rdist.le.r0ppt) then

                      phi(randvar,i,j,k) = 1.0
!                     do ivar = 1, vartot
!                       phi(ivar,i,j,k) = 1.0/float(20.0)
!                     end do

                    end if

                end do
                end do
                end do

                nppt = nppt + 1

            end if

            ! align
            call MPI_Barrier(MPI_COMM_WORLD, ierr)

! =====
            if(nppt.le.nmaxppt) go to 13

        end if
! ===== end - case: multiple


! ===== random nucleation
        if(nucstatus.eq.'random') then
            if(myid.eq.0) write(*,*) "Nucleus status: Random"

            do k = ist(3), ien(3)
            do j = ist(2), ien(2)
            do i = ist(1), ien(1)

                randvar = ran_2(iseed)*vartot + 1.0
                phi(randvar,i,j,k)=0.1

            end do
            end do
            end do
        end if
! ===== end - case: random
! ===== Set up: Phase values done

! ===== update: concentration fields
        if(myid .eq. 0) write(*,*) "Cocentration: ", constat

        if (constat.eq.'cppt') then
            if(myid .eq. 0) then
                write (*,*) "Change: concentration inside the target variant"
                write (*,*) "Target: ", targv
                write (*,*) "Outside PPT: ", (Cav(j), j=1, numel)
                write (*,*) "Inside PPT : ", (Cimax(j), j=1, numel)
            end if

            do k  = ist(3), ien(3)
            do j  = ist(2), ien(2)
            do i  = ist(1), ien(1)

                ! note: instead of using vartot, any other value can be substituted for a purpose
                if (phi(targv,i,j,k) .gt. 0.9) then

                    do ii = 1, numel

                        Cons(ii,i,j,k) = Cimax(ii)

                    end do


                endif

            end do
            end do
            end do

        endif

        return
        end subroutine InitFields

! ===== Done initilization
! =========================================================================

        subroutine NucAdd(nucstatus, phi)

        include 'mpif.h'

        character nucstatus*10
        integer r0ppt, rdist
        integer i, j, k, ii, jj, kk, ix, jy, kz, ivar
        integer r1, r2, r3
        integer mnx, mny, mnz
        integer pcount, nppt, pcount_tot
        integer randvar
        integer rboxnx, rboxny, rboxnz

        ! random number
        real *8 ran_2

        ! field informations
        real *8 phi(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))

! =====

        ! set useful parameters
        mnx = Nx/2
        mny = Ny/2
        mnz = Nz/2
        r0ppt = radppt(1)**2 + radppt(2)**2 + radppt(3)**2

        ! box size for nucleation
        rboxnx = Nx - 2 * (radppt(1) + addrad(1) + 1)
        rboxny = Ny - 2 * (radppt(2) + addrad(2) + 1)
        rboxnz = Nz - 2 * (radppt(3) + addrad(3) + 1)

        ! ===== Multiple nucleation case
        ! Initially three nuclei of different variants placed in the domain
        ! Subsequent nucleation is based on the composite nucleus model

        write(*,*) "Tempo: check add variant subroutine "

        if(nucstatus.eq.'addmult') then
            nppt = 1

            if(targv.eq.0) then
                if(myid.eq.0) write(*,*) "Nucleus status: Multiple"
            else
                if(myid.eq.0) write(*,*) "Nucleus status: Multiple with a variant", targv
            endif

13          if(myid.eq.0) then

16              ix = ran_2(iseed) * rboxnx + 1.0 + radppt(1) + addrad(1)
                jy = ran_2(iseed) * rboxny + 1.0 + radppt(2) + addrad(2)
                kz = ran_2(iseed) * rboxnz + 1.0 + radppt(3) + addrad(3)

17              if(targv.eq.0) then

                    ! at this moment do not consider random nucleation of laves phase
                    randvar = ran_2(iseed)*(vartot-1) + 1.0

                else

                    randvar = targv

                endif


                r1 = ran_2(iseed) * addrad(1)
                r2 = ran_2(iseed) * addrad(2)
                r3 = ran_2(iseed) * addrad(3)

            end if

            call MPI_Bcast(ix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(jy, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(kz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            call MPI_Bcast(r1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(r2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(r3, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            call MPI_Bcast(randvar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            ! ====== check whether precipitates overlap the others
            pcount = 0

            !=====
            r0ppt = (radppt(1) + r1)**2 + (radppt(2) + r2)**2 + (radppt(3) + r3)**2
            !=====

            do k  = ist(3), ien(3)
            do j  = ist(2), ien(2)
            do i  = ist(1), ien(1)

                rdist = (i - ix)**2 + (j - jy)**2 + (k-kz)**2

                do ivar = 1, vartot
                    if((rdist.le.r0ppt).AND.(phi(ivar,i,j,k).gt.0.01)) then
                        pcount = pcount +1
                    end if
                end do

            end do
            end do
            end do

            call MPI_Allreduce(pcount, pcount_tot, 1, MPI_INTEGER, &
                            MPI_SUM, MPI_COMM_WORLD, ierr)

            if(pcount_tot.gt.0) go to 16

! ===== change the precipitate
            if(pcount_tot.eq.0) then

                do k  = ist(3), ien(3)
                do j  = ist(2), ien(2)
                do i  = ist(1), ien(1)


                    rdist = (i - ix)**2 + (j - jy)**2 + (k-kz)**2

                    if(rdist.le.r0ppt) then

                        phi(randvar,i,j,k) = 1.0
                        phi(vartot,i,j,k) = 0.0

                    end if

                end do
                end do
                end do

                nppt = nppt + 1

            end if

            ! align
            call MPI_Barrier(MPI_COMM_WORLD, ierr)

            if(nppt.le.nmaxppt) go to 13

        end if
! ===== end - case:

      return
      end subroutine NucAdd

! ===== Done add nucleations

! =========================================================================

      end module ModInitials
