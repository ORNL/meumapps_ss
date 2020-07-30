      program data_ext_modular
      implicit none
      integer Nx, Ny, Nz, nprocs, dims(2), myid, ist(3), ien(3)
      integer i, j, k, ii, kk, var, ivar, step, imod, jmod, kmod, var_max, nel
      real *8, dimension(:,:,:,:), allocatable :: cons, phi
      real *8, dimension(:,:,:), allocatable :: elen
      real *8, dimension(:), allocatable :: phi_count
      real *8 term, phi_sq, phi_sum, phi_tot, phi_max, phi_temp
      character file_num*4

      open(11, file='data_extract_mult.in', status='old')
      read(11,*) Nx, Ny, Nz, nprocs, var, nel
      read(11,*) dims(1), dims(2)
      close(11)
      write(95,*) Nx, Ny, Nz, nprocs, var, dims(1), dims(2)
      call flush(95)
      open(1, file='shape_con', status='unknown')
      open(3, file='shape_phi', status='unknown')
      open(4, file='phi_var', status='unknown')

      allocate(phi_count(var))

      phi_tot = 0.0
      phi_count = 0.0

      write(1,*) 'zone', ' ', 'i=', Nx, ' ', 'j=', Ny, ' ', 'k=', Nz, 'f=', 'point'
      write(3,*) 'zone', ' ', 'i=', Nx, ' ', 'j=', Ny, ' ', 'k=', Nz, 'f=', 'point'

      myid = 0
      do kk=1,dims(2)
         do j=1,dims(1)
            myid = myid + 1
            write(file_num,'(i4.4)') myid
            open(12,file='data_save.'//file_num, status='old', &
                form='unformatted')
            read(12) ist, ien
            if ( j.eq. 1 ) allocate( cons(nel, Nx, Ny, ist(3):ien(3)) )
            if ( j.eq. 1 ) allocate( phi(var,Nx, Ny, ist(3):ien(3)) )
            if ( j.eq. 1 ) allocate( elen(Nx, Ny, ist(3):ien(3)) )
            read(12) cons(1:nel,1:Nx,ist(2):ien(2),ist(3):ien(3)),  &
                     phi(1:var,1:Nx,ist(2):ien(2),ist(3):ien(3)), &
                     elen(1:Nx,ist(2):ien(2),ist(3):ien(3))
            close(12)
         end do

         do k = ist(3), ien(3)
         do j = 1, Ny
         do i = 1, Nx
           write(1,3) i, j, k, (cons(ii,i,j,k), ii=1,nel), elen(i,j,k)
           if(k.eq.Nz/2.AND.j.eq.Ny/2) write(22,4) i,(cons(ii,i,j,k),ii=1,nel)

         end do
         end do
         end do
         
         do k = ist(3), ien(3)
         do j = 1, Ny
         do i = 1, Nx

         phi_sq = 0.0
         phi_sum = 0.0
         phi_max = 0.85
         var_max = 0

         do ivar = 1, var
            phi_temp = phi(ivar,i,j,k)
            if(phi_temp.gt.phi_max) then
!              phi_max = phi_temp
               var_max = ivar
            end if
!           if(phi_temp.gt.0.3) var_max = ivar
            phi_sq = phi_sq+ phi_temp**2
            phi_sum = phi_sum + phi_temp
            phi_count(ivar) = phi_count(ivar) + phi_temp
         end do 
!        write(*,*) 'var_max=', var_max
         phi_tot = phi_tot + phi_sum 
!        end if

         write(3,5) i, j, k, var_max, phi_sq, phi_sum

7        continue
         
         end do
         end do
         end do

         deallocate (cons)
         deallocate (phi)
         deallocate (elen)

      end do

      write(*,*) phi_count, phi_tot
      phi_count = phi_count/phi_tot
      write(4,6) (phi_count(ivar),ivar=1, var)
1     format("zone i=",i3," j=",i3," k=",i3," f=point")
3     format(3i6,4e15.7)
4     format(i6,e15.7,e15.7)
5     format(4i6,18e15.7)
6     format(12f8.4)
      end
