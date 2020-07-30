      module ModElastic

      use ModGparams
      use ModFunc
      implicit none
! =========================================================================

      contains
! eldis, elener
! =========================================================================
        subroutine eldis(phi, el_en, el_grad, u1, u2, u3, &
                        grad_1_u1, grad_1_u2, grad_1_u3, &
                        grad_2_u1, grad_2_u2, grad_2_u3, &
                        grad_3_u1, grad_3_u2, grad_3_u3)

        include 'mpif.h'

        ! declare parameters
        integer ii, jj, kk, i, j, k, l, ij, kl
        integer ivn, fvn, ivar
        integer niter
        integer nphc

        real *8 C_mean(6,6)
        real *8 gradu(3,3), C(6,6), C_prime(6,6)
        real *8 Qinv(3, 3, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))

        ! phi
        real *8 pphi_sum, sum_0
        real *8 phi(vartot, ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))

        ! energy related
        real *8 el_en(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 el_grad(vartot,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))

        real *8 u1(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 u2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 u3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))

        real *8 u1_old(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 u2_old(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 u3_old(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))

        real *8 grad_1_u1(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_1_u2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_1_u3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_2_u1(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_2_u2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_2_u3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_3_u1(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_3_u2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_3_u3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))

        real *8 x_1(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 x_2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 x_3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 x_4(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 x_5(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 x_6(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))

        ! FFT - complex fields
        double complex dft_P_phi(vartot, fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))

        double complex dft_u1(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_u2(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_u3(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))

        double complex dft_grad_1_u1(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_grad_2_u1(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_grad_3_u1(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))

        double complex dft_grad_1_u2(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_grad_2_u2(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_grad_3_u2(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))

        double complex dft_grad_1_u3(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_grad_2_u3(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_grad_3_u3(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))

        double complex dft_x_1(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_x_2(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_x_3(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_x_4(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_x_5(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))
        double complex dft_x_6(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))

        real *8 dummy(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 Q(3,3), err_loc(6), err(6), err_old(6)
        real *8 cfac, newfac

        double complex R(3), kf_j, kpr
        double complex dft_dummy(fst(1):fen(1), fst(2):fen(2), fst(3):fen(3))

        ! =====
        logical converged

        ! end declaring parameters
        ! =====

        newfac = 0.5d0
        niter = 0
        cfac = 1.d0
        converged = .false.

        do ivar = 1, vartot

        ! trans_mat is always 0 in the current version: trans = trans_ppt
        ! trans(:,ivar) = trans_ppt(:,ivar) - trans_mat(:)

            call f_trans( P_phi(ivar,:,:,:), dft_P_phi(ivar,:,:,:) )

        end do

        ! =====
        do kk = fst(3), fen(3)
        do jj = fst(2), fen(2)
        do ii = fst(1), fen(1)

            Q = 0.d0
            R = (0.d0,0.d0)

            do i=1,3
            do j=1,3
                ij = i

                if ( i.ne.j) ij = 9 - i - j
                kf_j = kf(j,ii,jj,kk)

                sum_0 = 0.d0

                do k=1,3
                do l=1,3
                    kl = k
                    if ( k.ne.l ) kl = 9 - k - l

                    kpr = kf(k,ii,jj,kk) * kf_j

                    if ( k.eq.j ) kpr = kf_sq(k,ii,jj,kk)

                    ! Q(i,l) = Q(i,l) + C_alph(ij,kl,4)*kpr
                    ! need to update: currently using delta.
                    Q(i,l) = Q(i,l) + CEla(ij,kl,2) * kpr

                    ivn = 0
                    do nphc = 2, numph
                        ivn = ivn + Nvar(nphc - 1)
                        fvn = ivn + Nvar(nphc) -1

                        do ivar = ivn, fvn

                            sum_0 = sum_0 + dft_P_phi(ivar,ii,jj,kk) * trans_ppt(kl,ivar) * CEla(ij,kl,nphc)


                        end do
                    end do

                end do
                end do

                R(i) = R(i) + sum_0 * kf_j

            end do
            end do

            ! =====
            if ( ii.EQ.1 .AND. jj.EQ.1 .AND. kk.EQ.1 ) then
                if ( step .eq. elastep0 ) then

                    dft_u1(ii,jj,kk) = 0.d0
                    dft_u2(ii,jj,kk) = 0.d0
                    dft_u3(ii,jj,kk) = 0.d0

                    dft_grad_1_u1(ii,jj,kk) = (0.d0, 0.d0)
                    dft_grad_1_u2(ii,jj,kk) = (0.d0, 0.d0)
                    dft_grad_1_u3(ii,jj,kk) = (0.d0, 0.d0)

                    dft_grad_2_u1(ii,jj,kk) = (0.d0, 0.d0)
                    dft_grad_2_u2(ii,jj,kk) = (0.d0, 0.d0)
                    dft_grad_2_u3(ii,jj,kk) = (0.d0, 0.d0)

                    dft_grad_3_u1(ii,jj,kk) = (0.d0, 0.d0)
                    dft_grad_3_u2(ii,jj,kk) = (0.d0, 0.d0)
                    dft_grad_3_u3(ii,jj,kk) = (0.d0, 0.d0)

                end if
            else
                err_flag = 3

                call invert(Q, 3, err_flag)

                Qinv(:,:,ii,jj,kk) = Q(:,:)

                if ( step .eq. elastep0 ) then
                    dft_u1(ii,jj,kk) = Q(1,1)*R(1) + Q(1,2)*R(2) + Q(1,3)*R(3)
                    dft_u2(ii,jj,kk) = Q(2,1)*R(1) + Q(2,2)*R(2) + Q(2,3)*R(3)
                    dft_u3(ii,jj,kk) = Q(3,1)*R(1) + Q(3,2)*R(2) + Q(3,3)*R(3)

                    dft_grad_1_u1(ii,jj,kk) = kf(1,ii,jj,kk)*dft_u1(ii,jj,kk)
                    dft_grad_1_u2(ii,jj,kk) = kf(1,ii,jj,kk)*dft_u2(ii,jj,kk)
                    dft_grad_1_u3(ii,jj,kk) = kf(1,ii,jj,kk)*dft_u3(ii,jj,kk)

                    dft_grad_2_u1(ii,jj,kk) = kf(2,ii,jj,kk)*dft_u1(ii,jj,kk)
                    dft_grad_2_u2(ii,jj,kk) = kf(2,ii,jj,kk)*dft_u2(ii,jj,kk)
                    dft_grad_2_u3(ii,jj,kk) = kf(2,ii,jj,kk)*dft_u3(ii,jj,kk)

                    dft_grad_3_u1(ii,jj,kk) = kf(3,ii,jj,kk)*dft_u1(ii,jj,kk)
                    dft_grad_3_u2(ii,jj,kk) = kf(3,ii,jj,kk)*dft_u2(ii,jj,kk)
                    dft_grad_3_u3(ii,jj,kk) = kf(3,ii,jj,kk)*dft_u3(ii,jj,kk)

                end if
            end if

        end do
        end do
        end do

        ! =====
        if ( step .eq. elastep0 ) then

            call inv_trans(dft_u1, u1)
            call inv_trans(dft_u2, u2)
            call inv_trans(dft_u3, u3)

            call inv_trans(dft_grad_1_u1, grad_1_u1)
            call inv_trans(dft_grad_1_u2, grad_1_u2)
            call inv_trans(dft_grad_1_u3, grad_1_u3)

            call inv_trans(dft_grad_2_u1, grad_2_u1)
            call inv_trans(dft_grad_2_u2, grad_2_u2)
            call inv_trans(dft_grad_2_u3, grad_2_u3)

            call inv_trans(dft_grad_3_u1, grad_3_u1)
            call inv_trans(dft_grad_3_u2, grad_3_u2)
            call inv_trans(dft_grad_3_u3, grad_3_u3)

        end if
        ! =====

        ! =====

2       u1_old = u1
        u2_old = u2
        u3_old = u3

        niter = niter + 1

        ! =====
        do kk = ist(3), ien(3)
        do jj = ist(2), ien(2)
        do ii = ist(1), ien(1)

            pphi_sum = 0.d0
            C = 0.d0

            ivn = 0
            do nphc = 2, numph
                ivn = ivn + Nvar(nphc - 1)
                fvn = ivn + Nvar(nphc) -1

                do ivar = ivn, fvn

                    pphi_sum = pphi_sum + P_phi(ivar,ii,jj,kk)
                    C = C + P_phi(ivar,ii,jj,kk) * SEla(:,:,nphc)

                end do
            end do

            if(pphi_sum.lt.1.d0) C = C + (1.d0 - pphi_sum) * SEla(:,:,1)

            err_flag = 4

            call invert(C, 6, err_flag)

            ! need to update later: currently using delta.
            C_prime = CEla(:,:,2) - C(:,:)

            gradu(1,1) = grad_1_u1(ii,jj,kk)
            gradu(1,2) = grad_2_u1(ii,jj,kk)
            gradu(1,3) = grad_3_u1(ii,jj,kk)

            gradu(2,1) = grad_1_u2(ii,jj,kk)
            gradu(2,2) = grad_2_u2(ii,jj,kk)
            gradu(2,3) = grad_3_u2(ii,jj,kk)

            gradu(3,1) = grad_1_u3(ii,jj,kk)
            gradu(3,2) = grad_2_u3(ii,jj,kk)
            gradu(3,3) = grad_3_u3(ii,jj,kk)

            ! =====
            do i=1,3
            do j=i,3

                ij = i
                if ( i.ne.j) ij = 9 - i - j

                sum_0 = 0.d0

                do k=1,3
                do l=1,3

                    kl = k
                    if ( k.ne.l ) kl = 9 - k - l

                    sum_0 = sum_0 + cfac * C_prime(ij,kl) * ( gradu(k,l) + e_mean(kl) )

                    do ivar = 1, vartot

                        sum_0 = sum_0 + C(ij,kl)*P_phi(ivar,ii,jj,kk)*trans_ppt(kl,ivar)

                    end do

                end do
                end do

                if ( ij .eq. 1 ) x_1(ii,jj,kk) = sum_0
                if ( ij .eq. 2 ) x_2(ii,jj,kk) = sum_0
                if ( ij .eq. 3 ) x_3(ii,jj,kk) = sum_0
                if ( ij .eq. 4 ) x_4(ii,jj,kk) = sum_0
                if ( ij .eq. 5 ) x_5(ii,jj,kk) = sum_0
                if ( ij .eq. 6 ) x_6(ii,jj,kk) = sum_0

            end do
            end do

        end do
        end do
        end do

        call f_trans(x_1, dft_x_1)
        call f_trans(x_2, dft_x_2)
        call f_trans(x_3, dft_x_3)
        call f_trans(x_4, dft_x_4)
        call f_trans(x_5, dft_x_5)
        call f_trans(x_6, dft_x_6)

        ! u1 and u2 calculated above are the zeroth order approximations
        ! they will be used as inputs to successive higher order calculations
        ! until they converge

        do kk = fst(3), fen(3)
        do jj = fst(2), fen(2)
        do ii = fst(1), fen(1)

            R = (0.d0, 0.d0)

            R(1) = dft_x_1(ii,jj,kk) * kf(1,ii,jj,kk)
            R(1) = R(1) + dft_x_6(ii,jj,kk) * kf(2,ii,jj,kk)
            R(1) = R(1) + dft_x_5(ii,jj,kk) * kf(3,ii,jj,kk)

            R(2) = dft_x_6(ii,jj,kk) * kf(1,ii,jj,kk)
            R(2) = R(2) + dft_x_2(ii,jj,kk) * kf(2,ii,jj,kk)
            R(2) = R(2) + dft_x_4(ii,jj,kk) * kf(3,ii,jj,kk)

            R(3) = dft_x_5(ii,jj,kk)*kf(1,ii,jj,kk)
            R(3) = R(3) + dft_x_4(ii,jj,kk)*kf(2,ii,jj,kk)
            R(3) = R(3) + dft_x_3(ii,jj,kk)*kf(3,ii,jj,kk)

            if ( ii.EQ.1 .AND. jj.EQ.1 .AND. kk.EQ.1 ) then

                dft_u1(ii,jj,kk) = 0.d0
                dft_u2(ii,jj,kk) = 0.d0
                dft_u3(ii,jj,kk) = 0.d0

            else

                dft_u1(ii,jj,kk) = Qinv(1,1,ii,jj,kk)*R(1) &
                                + Qinv(1,2,ii,jj,kk)*R(2) &
                                + Qinv(1,3,ii,jj,kk)*R(3)

                dft_u2(ii,jj,kk) = Qinv(2,1,ii,jj,kk)*R(1) &
                                + Qinv(2,2,ii,jj,kk)*R(2) &
                                + Qinv(2,3,ii,jj,kk)*R(3)

                dft_u3(ii,jj,kk) = Qinv(3,1,ii,jj,kk)*R(1) &
                                + Qinv(3,2,ii,jj,kk)*R(2) &
                                + Qinv(3,3,ii,jj,kk)*R(3)

            end if

            dft_grad_1_u1(ii,jj,kk) = kf(1,ii,jj,kk) * dft_u1(ii,jj,kk)
            dft_grad_1_u2(ii,jj,kk) = kf(1,ii,jj,kk) * dft_u2(ii,jj,kk)
            dft_grad_1_u3(ii,jj,kk) = kf(1,ii,jj,kk) * dft_u3(ii,jj,kk)

            dft_grad_2_u1(ii,jj,kk) = kf(2,ii,jj,kk) * dft_u1(ii,jj,kk)
            dft_grad_2_u2(ii,jj,kk) = kf(2,ii,jj,kk) * dft_u2(ii,jj,kk)
            dft_grad_2_u3(ii,jj,kk) = kf(2,ii,jj,kk) * dft_u3(ii,jj,kk)

            dft_grad_3_u1(ii,jj,kk) = kf(3,ii,jj,kk) * dft_u1(ii,jj,kk)
            dft_grad_3_u2(ii,jj,kk) = kf(3,ii,jj,kk) * dft_u2(ii,jj,kk)
            dft_grad_3_u3(ii,jj,kk) = kf(3,ii,jj,kk) * dft_u3(ii,jj,kk)

        end do
        end do
        end do

        call inv_trans(dft_u1, u1)
        call inv_trans(dft_u2, u2)
        call inv_trans(dft_u3, u3)

        call inv_trans(dft_grad_1_u1, grad_1_u1)
        call inv_trans(dft_grad_1_u2, grad_1_u2)
        call inv_trans(dft_grad_1_u3, grad_1_u3)

        call inv_trans(dft_grad_2_u1, grad_2_u1)
        call inv_trans(dft_grad_2_u2, grad_2_u2)
        call inv_trans(dft_grad_2_u3, grad_2_u3)

        call inv_trans(dft_grad_3_u1, grad_3_u1)
        call inv_trans(dft_grad_3_u2, grad_3_u2)
        call inv_trans(dft_grad_3_u3, grad_3_u3)

        ! =====
        err_loc = 0.d0

        do k = ist(3), ien(3)
        do j = ist(2), ien(2)
        do i = ist(1), ien(1)

            err_loc(1) = err_loc(1) + (u1(i,j,k) - u1_old(i,j,k))**2
            err_loc(2) = err_loc(2) + u1_old(i,j,k)**2

            err_loc(3) = err_loc(3) + (u2(i,j,k) - u2_old(i,j,k))**2
            err_loc(4) = err_loc(4) + u2_old(i,j,k)**2

            err_loc(5) = err_loc(5) + (u3(i,j,k) - u3_old(i,j,k))**2
            err_loc(6) = err_loc(6) + u3_old(i,j,k)**2

        end do
        end do
        end do

        call MPI_Allreduce(err_loc(1), err(1), 6, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, MPI_COMM_WORLD, ierr)

        if(err(1)/err(2).lt.1.d-7 .AND. err(3)/err(4).lt.1.d-7 .AND. &
                    err(5)/err(6).lt.1.d-7 ) converged = .true.

        if ( niter .gt. 1 .AND. .not. converged ) then
            if ( err(1)/err(2) .gt. err_old(1)/err_old(2) .OR. &
                err(3)/err(4) .gt. err_old(3)/err_old(4) .OR. &
                err(5)/err(6) .gt. err_old(5)/err_old(6) ) then

                cfac = cfac * newfac
                newfac = newfac * 0.5d0

                if ( newfac .lt. 1.d-7 ) then
                    write(98,*) 'newfac is too small ... trouble converging'
                    write(*,*) converged, myid
                    STOP
                end if

            niter = 0
            converged = .false.

            go to 2

            end if
        end if

        err_old = err

        if ( niter .eq. 10 ) then

            if ( myid .eq. 0 ) write(98,*) 'Decreasing error over 10 iterations'
            converged = .true.

        end if

        if ( .not. converged .OR. niter.lt.4 ) go to 2

        if( cfac .lt. 1.d0-1.d-6 ) then

            cfac = cfac + 0.1d0
            if ( cfac .gt. 1.d0 ) cfac = 1.d0

            niter = 0
            converged = .false.

            go to 2

        end if

        call elener(phi, el_en, el_grad, grad_1_u1, grad_1_u2, grad_1_u3, &
                    grad_2_u1, grad_2_u2, grad_2_u3, &
                    grad_3_u1, grad_3_u2, grad_3_u3)

        return

        end subroutine eldis

! =========================================================================

        subroutine elener(phi, el_en, el_grad, grad_1_u1, grad_1_u2, grad_1_u3, &
                        grad_2_u1, grad_2_u2, grad_2_u3, &
                        grad_3_u1, grad_3_u2, grad_3_u3)

        ! declare parameters

        ! domain
        integer i, j, k, ii, jj
        integer ivn, fvn, ivar
        integer nphc

        real *8 pphi_sum, pphiloc
        real *8 trans_eps(6)
        real *8 transsum(vartot)

        ! phi
        real *8 phi(vartot,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))

        ! energy related
        real *8 grad_1_u1(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_1_u2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_1_u3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))

        real *8 grad_2_u1(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_2_u2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_2_u3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))

        real *8 grad_3_u1(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_3_u2(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 grad_3_u3(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))

        real *8 el_en(ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 el_grad(vartot,ist(1):ien(1), ist(2):ien(2), ist(3):ien(3))
        real *8 el_ppt, el_mat

        ! additional
        real *8 eps_ppt(6), eps_mat(6)
        real *8 C(6,6), eps(6), sig(6)

        ! end declaring parameters

        transsum = sum(trans_ppt, dim = 1)

        if(myid.eq.0) write(95,*) 'trans_ppt(:,var)'
        if(myid.eq.0) write(95,*) -1*trans_ppt(:,:)
        call flush(95)

        do k = ist(3), ien(3)
        do j = ist(2), ien(2)
        do i = ist(1), ien(1)

            !trans_eps= trans_mat
            trans_eps = 0.d0
            pphi_sum = 0.d0

            ! Precipitates: variants
            do ivar = 1, vartot

                fi = phi(ivar,i,j,k)
                fi_sq = fi * fi

                pphiloc = P_phi(ivar,i,j,k)
                pphi_sum = pphi_sum + pphiloc

                Ppr_phi(ivar) = 30.d0 * fi_sq * (fi_sq - 2.d0 * fi + 1.d0)

                trans_eps(1) = trans_eps(1) + pphiloc * trans_ppt(1,ivar)
                trans_eps(2) = trans_eps(2) + pphiloc * trans_ppt(2,ivar)
                trans_eps(3) = trans_eps(3) + pphiloc * trans_ppt(3,ivar)
                trans_eps(4) = trans_eps(4) + pphiloc * trans_ppt(4,ivar)
                trans_eps(5) = trans_eps(5) + pphiloc * trans_ppt(5,ivar)
                trans_eps(6) = trans_eps(6) + pphiloc * trans_ppt(6,ivar)


            end do

            ! Compute elastic strain tensor

            eps(1) = e_mean(1) + grad_1_u1(i,j,k) - trans_eps(1)
            eps(2) = e_mean(2) + grad_2_u2(i,j,k) - trans_eps(2)
            eps(3) = e_mean(3) + grad_3_u3(i,j,k) - trans_eps(3)

            eps(4) = 2.d0 * e_mean(4) + grad_2_u3(i,j,k) + grad_3_u2(i,j,k) &
                    - 2.d0 * trans_eps(4)
            eps(5) = 2.d0 * e_mean(5) + grad_1_u3(i,j,k) + grad_3_u1(i,j,k) &
                    - 2.d0 * trans_eps(5)
            eps(6) = 2.d0 * e_mean(6) + grad_2_u1(i,j,k) + grad_1_u2(i,j,k) &
                    - 2.d0 * trans_eps(6)

            ! =====
            C = 0.0

            ivn = 0
            do nphc = 2, numph
                ivn = ivn + Nvar(nphc - 1)
                fvn = ivn + Nvar(nphc) -1

                do ivar = ivn, fvn

                    C = C + P_phi(ivar,i,j,k) * SEla(:,:,nphc)

                end do
            end do

            C = C  + (1.d0 - pphi_sum)*SEla(:,:,1)

            err_flag = 6

            call invert(C, 6, err_flag)

            ! =====
            ! this loop and the below loop can be combined: tried but it does not work
            !  Compute stress from elastic strain

            do ii = 1,6
                sig(ii) = 0.d0

                sig(ii) = sig(ii) + C(ii,1) * eps(1)
                sig(ii) = sig(ii) + C(ii,2) * eps(2)
                sig(ii) = sig(ii) + C(ii,3) * eps(3)
                sig(ii) = sig(ii) + C(ii,4) * eps(4)
                sig(ii) = sig(ii) + C(ii,5) * eps(5)
                sig(ii) = sig(ii) + C(ii,6) * eps(6)

            end do

            ! =====
            ! Compute elastic strain energy - matrix

            el_mat = 0.d0

            do ii = 1,6
                eps_mat(ii) = 0.d0

                eps_mat(ii) = eps_mat(ii) + SEla(ii,1,1) * sig(1)
                eps_mat(ii) = eps_mat(ii) + SEla(ii,2,1) * sig(2)
                eps_mat(ii) = eps_mat(ii) + SEla(ii,3,1) * sig(3)
                eps_mat(ii) = eps_mat(ii) + SEla(ii,4,1) * sig(4)
                eps_mat(ii) = eps_mat(ii) + SEla(ii,5,1) * sig(5)
                eps_mat(ii) = eps_mat(ii) + SEla(ii,6,1) * sig(6)

                el_mat = el_mat + sig(ii) * eps_mat(ii)

            end do

            el_mat = 0.5d0 * el_mat

            el_en(i,j,k) = (1.d0 - pphi_sum) * el_mat

            ! =====
            ! Compute elastic strain energy - precipitates

            ivn = 0

            do nphc = 2, numph
                ivn = ivn + Nvar(nphc - 1)
                fvn = ivn + Nvar(nphc) -1

                do ivar = ivn, fvn

                if (transsum(ivar) .gt. 1.e-9) then
                    el_ppt = 0.d0

                    do ii = 1,6

                        eps_ppt(ii) = 0.d0

                        eps_ppt(ii) = eps_ppt(ii) + SEla(ii,1,nphc) * sig(1)
                        eps_ppt(ii) = eps_ppt(ii) + SEla(ii,2,nphc) * sig(2)
                        eps_ppt(ii) = eps_ppt(ii) + SEla(ii,3,nphc) * sig(3)
                        eps_ppt(ii) = eps_ppt(ii) + SEla(ii,4,nphc) * sig(4)
                        eps_ppt(ii) = eps_ppt(ii) + SEla(ii,5,nphc) * sig(5)
                        eps_ppt(ii) = eps_ppt(ii) + SEla(ii,6,nphc) * sig(6)

                        el_ppt = el_ppt + sig(ii)*eps_ppt(ii)

                    end do

                    el_ppt = 0.5d0 * el_ppt

                    el_en(i,j,k) = el_en(i,j,k) + el_ppt * P_phi(ivar,i,j,k)

                    el_grad(ivar,i,j,k) = ( sig(1)*trans_ppt(1,ivar) + sig(2) * trans_ppt(2,ivar) + sig(3) * trans_ppt(3,ivar) &
                                    + 2.d0 * (sig(4) * trans_ppt(4,ivar) + sig(5)*trans_ppt(5,ivar) + sig(6)*trans_ppt(6,ivar)) )

                    el_grad(ivar,i,j,k) = Ppr_phi(ivar) * (el_mat - el_ppt) - el_grad(ivar,i,j,k) * Ppr_phi(ivar)

                end if
                end do
            end do


        end do
        end do
        end do

        ! =====
        return
        end subroutine elener

! =========================================================================

      end module ModElastic
