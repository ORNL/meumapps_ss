      module ModFunc

      use ModGparams
      implicit none
! =========================================================================

! =========================================================================
      contains
! f_trans, inv_trans, invert, and langevin
! ===== ran2: not working - update later

! =========================================================================
! ===== f_trans

        subroutine f_trans(Ain, Aout)

        use p3dfft

        real *8 Ain( ist(1):ien(1), ist(2):ien(2), ist(3):ien(3) )
        double complex Aout( fst(1):fen(1), fst(2):fen(2), fst(3):fen(3) )

        call p3dfft_ftran_r2c(Ain, Aout, 'fft')
        Aout = Aout * ooNxyz

        return
        end subroutine f_trans

! =========================================================================
! ===== inv_trans
        subroutine inv_trans(Ain, Aout)

        use p3dfft

        double complex Ain( fst(1):fen(1), fst(2):fen(2), fst(3):fen(3) )
        real *8 Aout( ist(1):ien(1), ist(2):ien(2), ist(3):ien(3) )

        call p3dfft_btran_c2r(Ain, Aout, 'tff')

        return
        end subroutine inv_trans

! =========================================================================
! ===== invert
        subroutine invert(c,n,flag)

        integer n, job, ipvt(n), i, j, flag
        real *8 c(n,n), c_prime(n,n), z(n), rcond, det

        ! routines for inverting the C matrix using linpack library

        c_prime = c
        call dgeco(c,n,n,ipvt,rcond,z)
        job = 11

        if(rcond.gt.1.d-16) then
            call dgedi(c,n,n,ipvt,det,z,job)
        else
            write(*,*) 'rcond=', rcond, 'flag=', flag
            write(*,*) 'Poor cond. no; halting program'
            write(*,*) 'det=', det

            do i = 1, 4
                write(*,3) (c_prime(i,j), j=1,n)
            end do

            stop
        end if

3       format(4e15.7)

        return
        end subroutine invert


! =========================================================================
! ===== langevin: noise
        subroutine langevin(iseed,sd,harvest)

        integer i, iseed
        real *8 v1, v2, rsq, ran_2, sd, harvest

        harvest = 0.0

1       v1 = ran_2(iseed)
        v2 = ran_2(iseed)

        v1 = 2.d0*v1 - 1.d0
        v2 = 2.d0*v2 - 1.d0

        rsq = v1*v1 + v2*v2

        if(rsq.gt.0.d0.AND.rsq.lt.1.d0) then
            rsq = sqrt(-2.d0*log(rsq)/rsq)
            harvest = v2*rsq*sd
        else
            go to 1
        end if

        return
        end subroutine langevin
! =========================================================================


! =========================================================================
      end module ModFunc
