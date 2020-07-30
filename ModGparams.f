      module ModGparams

      implicit none
! =========================================================================
! simulation parameters
      integer Nx, Ny, Nz
      integer step, Nt, noise_step
      real *8 del_N, del_N_sq, ooNxyz, disc_ener
      real *8 b_width

      integer initerr, err_flag
      integer el_visit, nrun, elastep0

! for mpi / p3dfft / complex fields
      integer ndim, dims(2)
      integer ist(3), ien(3), isize(3)
      integer fst(3), fen(3), fsize(3)
      integer nprocs, myid, ierr

! for nucleation / precipitates / variants
      integer iseed, targv, nmaxppt
      integer radppt(3), addrad(3)
      integer vartot, var_lo, var_hi

! for phi related
      integer numph
      real *8 phi_init
      real *8 fi, fi_sq

! for elements:concentrations
      integer zdim
      integer numel
      real *8 diff_coeff, xdiffph

! for randum values
      real *8 sd, noise

! for energy related
      integer n_tilda
      real *8 fw, fwloc, fw_tot
      real *8 fch, fchloc, fch_tot, solid_frac, solid_tot
      real *8 b_mob
      real *8 well_grad
      real *8 f_int, fint_tot, intfantiph
      real *8 el_sum, elen_tot, e_tot
! for strain/stress matrices
      real *8 e_mean(6)

! (almost) not using:
      real *8 m_star, l_star

! additional
      character nucstat*10, file_num*4, dir_name*12, buffchar*16
      character constat*16
      logical file_exists

! Arrays =================================
! for setting
      integer, dimension(:), allocatable :: Nvar
      real *8, dimension(:), allocatable :: Zr
      real *8, dimension(:,:), allocatable :: Zinv, Zintf
      real *8, dimension(:,:), allocatable :: Xintf

! for phase
      real *8, dimension(:), allocatable :: Pph, Gph
      real *8, dimension(:), allocatable :: Ppr_phi
      real *8, dimension(:,:,:,:), allocatable :: P_phi

! for concentrations
      real *8, dimension(:), allocatable :: Cav, Cimax, Cigrad
      real *8, dimension(:), allocatable :: Cmobs, MuC


! for FFT related / complex fields
      double complex, dimension(:,:,:,:), allocatable :: kf, kf_sq

! for energy related
      real *8, dimension(:), allocatable :: Wellph, SumE

      real *8, dimension(:,:), allocatable :: Gparams
      real *8, dimension(:,:), allocatable :: EIntfph
      real *8, dimension(:,:), allocatable :: TStrph

      real *8, dimension(:,:,:), allocatable :: grad_coeff
! for stress / strain matrices
      integer, dimension(:), allocatable :: Ctype
      real *8, dimension(:,:,:), allocatable :: CEla
      real *8, dimension(:,:,:), allocatable :: SEla

      real *8, dimension(:,:), allocatable :: trans_ppt

! =========================================================================
      contains

! =========================================================================
      end module ModGparams

