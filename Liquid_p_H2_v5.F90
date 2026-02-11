program Liquid_p_H2_v5
   use iso_fortran_env
   use omp_lib
   implicit none

   integer,  parameter :: dp    = real64
   integer,  parameter :: n     = 108
   real(dp), parameter :: pi    = 3.14159265358979323846264338327950288419_dp
   
   real(dp), parameter :: kB_SI = 1.380649e-23_dp
   real(dp), parameter :: amu_to_kg = 1.66054e-27_dp 
   real(dp), parameter :: mass_H2_amu = 2.016_dp   
   real(dp), parameter :: mass_H2_kg = mass_H2_amu * amu_to_kg  
   real(dp), parameter :: Ang_to_m = 1.0e-10_dp  
   real(dp), parameter :: ps_to_s = 1.0e-12_dp     
   
   ! Derived conversion factor for F[K/Å] -> a[Å/ps^2]
   ! a = F · k_B / m, with unit conversions:
   ! F[N] = F[K/Å] · k_B[J/K] / (10^-10 m/Å) = F[K/Å] · k_B · 10^10
   ! a[m/s^2] = F[N] / m[kg] = F[K/Å] · k_B × 10^10 / m 
   ! a[Å/ps^2] = a[m/s^2] · 10^10 Å/m / (10^12 ps/s)^2 = a[m/s^2] · 10^-14
   ! Therefore: a[Å/ps^2] = F[K/Å] · (k_B / m) · 10^-4
   real(dp), parameter :: force_to_accel = kB_SI / mass_H2_kg * 1.0e-4_dp
   ! This equals 0.4124 Å^2/(ps^2·K)
   
   real(dp), parameter :: temp  = 14.0_dp      
   real(dp), parameter :: tau   = 0.0025_dp        
   real(dp), parameter :: dt    = 0.0005_dp          ! Time step in ps (0.5 fs)
   integer,  parameter :: equil_steps = 200000        ! Equilibration steps (10 ps)
   integer,  parameter :: steps = 20000               ! Production steps (4 ps)
   integer,  parameter :: ntau  = steps/2            ! Max time lag for VACF
   integer,  parameter :: ntraj = 1000  
   
   real(dp), parameter :: lattice_param = 5.54_dp
   real(dp), parameter :: box_length = 3.0_dp * lattice_param
   
   real(dp), parameter :: sg_alpha  = 1.713_dp    
   real(dp), parameter :: sg_beta   = 2.9615_dp      ! In 1/Angstrom
   real(dp), parameter :: sg_gamma  = 0.03546_dp   
   real(dp), parameter :: sg_rc     = 4.403_dp       ! Cutoff for damping function in Angstrom
   real(dp), parameter :: sg_C6     = 84178.9_dp     ! In K*Angstrom^6
   real(dp), parameter :: sg_C8     = 417858.4_dp  
   real(dp), parameter :: sg_C9     = 147037.3_dp   
   real(dp), parameter :: sg_C10    = 2617496.9_dp 
   real(dp), parameter :: cutoff    = 10.31_dp       ! Cutoff for intermolecular interaction

   real(dp)    :: r(3,n), v(3,n), f(3,n)
   real(dp)    :: vacf(ntau), vacf_sum(ntau), time_lag(ntau)
   real(dp),    allocatable :: v_series(:, :, :)     ! Velocity time series (steps, 3, n)
   integer     :: i, j, k, traj, coord
   real(dp)    :: D_integral, D_coefficient


   call random_seed()

   vacf = 0.0_dp
   vacf_sum = 0.0_dp

   call omp_set_num_threads(5)
   print *, 'Using', omp_get_max_threads(), 'OpenMP threads'

   do i = 1, ntau
      time_lag(i) = (i-1) * dt
   end do
   
   !$omp parallel do private(r, v, f, vacf, v_series, i, j) &
   !$omp& reduction(+:vacf_sum) schedule(dynamic)
   do traj = 1, ntraj
      !$omp critical
      print *, 'Trajectory', traj, '/', ntraj
      !$omp end critical

      allocate(v_series(steps, 3, n))

      call init_fcc(r, v, lattice_param, temp)
      
      do i = 1, equil_steps
         call step_vv(r, v, f, dt, box_length)
         call thermostat(v, dt, tau, temp)
      end do

      do i = 1, steps
         call step_vv(r, v, f, dt, box_length)
         v_series(i, :, :) = v(:, :)
      end do
     
      call compute_vacf(v_series, steps, n, vacf, ntau)
      vacf_sum = vacf_sum + vacf

      deallocate(v_series)
   end do
   !$omp end parallel do

   vacf = vacf_sum / real(ntraj, dp)
   call write_two_col('liquid_pH2_N256_VACF.dat', time_lag, vacf, ntau)
   
   ! D = (1/3) * integral(VACF)
   D_integral = 0.5_dp * vacf(1) * dt
   do i = 2, ntau-1
      D_integral = D_integral + vacf(i) * dt
   end do
   D_integral = D_integral + 0.5_dp * vacf(ntau) * dt
   
   D_coefficient = D_integral / 3.0_dp
   
   print *, 'VACF(0) =', vacf(1), 'Angstrom²/ps²'
   print *, 'D =', D_coefficient, 'Angstrom²/ps'


contains

   subroutine init_fcc(r, v, a, T)
      real(dp), intent(out) :: r(:,:), v(:,:)
      real(dp), intent(in)  :: a, T
      integer :: ix, iy, iz, idx
      real(dp) :: vcm(3)
      real(dp) :: v_scale
      
      real(dp), parameter :: basis(3,4) = reshape([ &
         0.0_dp, 0.0_dp, 0.0_dp, & 
         0.5_dp, 0.5_dp, 0.0_dp, &  
         0.5_dp, 0.0_dp, 0.5_dp, & 
         0.0_dp, 0.5_dp, 0.5_dp  & 
      ], [3, 4])
      
      idx = 0
      do iz = 0, 2
         do iy = 0, 2
            do ix = 0, 2
               do j = 1, 4
                  idx = idx + 1
                  r(1, idx) = (real(ix,dp) + basis(1,j)) * a
                  r(2, idx) = (real(iy,dp) + basis(2,j)) * a
                  r(3, idx) = (real(iz,dp) + basis(3,j)) * a
               end do
            end do
         end do
      end do
      
      v_scale = sqrt(force_to_accel * T)
      
      do i = 1, n
         do coord = 1, 3
            v(coord, i) = v_scale * randn()
         end do
      end do
      
      vcm = 0.0_dp
      do i = 1, n
         vcm(:) = vcm(:) + v(:, i)
      end do
      vcm = vcm / real(n, dp)
      
      do i = 1, n
         v(:, i) = v(:, i) - vcm(:)
      end do
      
   end subroutine init_fcc

   pure subroutine force_silvera_goldman(r, f, L)
      real(dp), intent(in)  :: r(:,:)
      real(dp), intent(out) :: f(:,:)
      real(dp), intent(in)  :: L
      integer  :: i, j
      real(dp) :: rij(3), rij_mag, rij_mag2
      real(dp) :: fij(3), fmag
      real(dp) :: V_rep, V_att, dV_rep, dV_att
      real(dp) :: r2, r4, r6, r8, r9, r10
      real(dp) :: fc, dfc
      
      f = 0.0_dp
      
      do i = 1, n
         do j = i+1, n  
            rij = r(:, j) - r(:, i)
            rij(1) = rij(1) - L * nint(rij(1) / L)
            rij(2) = rij(2) - L * nint(rij(2) / L)
            rij(3) = rij(3) - L * nint(rij(3) / L)
            
            rij_mag2 = rij(1)**2 + rij(2)**2 + rij(3)**2
            rij_mag = sqrt(rij_mag2)
            
            if (rij_mag < cutoff) then
               if (rij_mag <= sg_rc) then
                  fc = exp(-(sg_rc/rij_mag - 1.0_dp)**2)
                  dfc = fc * 2.0_dp * (sg_rc/rij_mag - 1.0_dp) * sg_rc / rij_mag2
               else
                  fc = 1.0_dp
                  dfc = 0.0_dp
               end if
               
               V_rep = exp(sg_alpha - sg_beta*rij_mag - sg_gamma*rij_mag2)
               dV_rep = -V_rep * (sg_beta + 2.0_dp*sg_gamma*rij_mag)
               
               r2 = rij_mag2
               r4 = r2 * r2
               r6 = r4 * r2
               r8 = r6 * r2
               r9 = r8 * rij_mag
               r10 = r8 * r2
               
               V_att = sg_C6/r6 + sg_C8/r8 - sg_C9/r9 + sg_C10/r10
               dV_att = -6.0_dp*sg_C6/r6 - 8.0_dp*sg_C8/r8 + 9.0_dp*sg_C9/r9 - 10.0_dp*sg_C10/r10
               dV_att = dV_att / rij_mag

               fmag = -(dV_rep - dfc*V_att - fc*dV_att)

               fij = (fmag / rij_mag) * rij
               
               f(:, i) = f(:, i) + fij
               f(:, j) = f(:, j) - fij
            end if
         end do
      end do
      
   end subroutine force_silvera_goldman

   pure subroutine step_vv(r, v, f, dt, L)
      ! Units: r[Å], v[Å/ps], f[K/Å], dt[ps]
      real(dp), intent(inout) :: r(:,:), v(:,:), f(:,:)
      real(dp), intent(in)    :: dt, L
      real(dp) :: halfdt
      integer  :: i, coord

      halfdt = 0.5_dp * dt
      call force_silvera_goldman(r, f, L)   
      v = v + halfdt * f * force_to_accel
      r = r + dt * v
      do i = 1, n
         do coord = 1, 3
            r(coord, i) = r(coord, i) - L * floor(r(coord, i) / L)
         end do
      end do
      call force_silvera_goldman(r, f, L)
      v = v + halfdt * f * force_to_accel
   end subroutine step_vv

   subroutine thermostat(v, dt, tau, T)
      real(dp), intent(inout) :: v(:,:)
      real(dp), intent(in)    :: dt, tau, T
      integer  :: dof, kshape
      real(dp) :: K, Kbar, c, s, r1, sum_r2, alpha2, alpha, factor
      real(dp) :: v_squared_sum
      dof  = 3 * n
      v_squared_sum = sum(v*v)
      K    = 0.5_dp * v_squared_sum / force_to_accel
      Kbar = 0.5_dp * real(dof,dp) * T
      c = exp(-dt/tau)
      s = 1.0_dp - c
      r1     = randn()
      kshape = (dof-2)/2
      sum_r2 = rand_gamma(kshape) + randn()**2
      factor = (Kbar / (real(dof,dp)*max(K, tiny(1.0_dp))))
      alpha2 = c                                            &
            + factor * s * (r1*r1 + sum_r2)               &
            + 2.0_dp * exp(-0.5_dp*dt/tau)                &
               * sqrt( factor * s ) * r1
      alpha2 = max(alpha2, 0.0_dp)
      alpha  = sqrt(alpha2)
      v      = alpha * v
   end subroutine thermostat

   real(dp) function rand_gamma(k)
      integer, intent(in) :: k
      real(dp) :: va(k)
      integer  :: i
      call random_number(va)
      do i = 1, k
         va(i) = max(va(i), 1.0e-12_dp)
      end do
      va = -log(va)
      rand_gamma = 2.0_dp * sum(va)
   end function rand_gamma

   real(dp) function randn()
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1    = max(u1, 1.0e-12_dp)
      randn = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*pi*u2)
   end function randn

   subroutine compute_vacf(v_series, nsteps, natoms, vacf_out, ntau_out)
      real(dp), intent(in)  :: v_series(:, :, :)  ! (steps, 3, natoms)
      integer,  intent(in)  :: nsteps, natoms, ntau_out
      real(dp), intent(out) :: vacf_out(:)
      integer :: i, j, k, coord, n_origins
      real(dp) :: dot_product_val

      vacf_out = 0.0_dp

      do i = 1, ntau_out
         n_origins = nsteps - (i-1)
         do j = 1, n_origins
            do k = 1, natoms
               dot_product_val = 0.0_dp
               do coord = 1, 3
                  dot_product_val = dot_product_val + v_series(j, coord, k) * v_series(j+(i-1), coord, k)
               end do
               vacf_out(i) = vacf_out(i) + dot_product_val
            end do
         end do
         vacf_out(i) = vacf_out(i) / (real(n_origins, dp) * real(natoms, dp))
      end do
   end subroutine compute_vacf

   subroutine write_two_col(fname, x, y, n)
      character(*), intent(in) :: fname
      real(dp),     intent(in) :: x(:), y(:)
      integer,      intent(in) :: n
      integer :: u,i
      open(newunit=u, file=fname, status='replace', action='write')
         write(u,'(A)') '# Time[ps]    VACF[Angstrom²/ps²]'
         do i=1,n
               write(u,*) x(i), y(i)
         end do
      close(u)
      print *, 'saved to:  ', fname
   end subroutine write_two_col

end program Liquid_p_H2_v5



