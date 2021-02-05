# BC: simplified version removing unused functions
f90_template = """module solve_ida

   implicit none

   integer :: neq = {neq}
   integer :: iout(50)
   real*8 :: rout(50)
   real*8 :: y0({neq}), yp0({neq})
   real*8 :: diff({neq}), mas({neq}, {neq})
   real*8 :: jac({neq}, {neq})
   real*8 :: rates({nrates})
   real*8 :: dypdr({neq}, {nrates})
   integer :: dvacdy({nvac}, {neq})

end module solve_ida

subroutine initialize(neqin, y0in, rtol, atol, ipar, rpar, id_vec)

   use solve_ida, only: neq, iout, rout, y0, yp0, mas, diff, dypdr, dvacdy

   implicit none

   integer, intent(in) :: neqin, ipar(*)
   real*8, intent(in) :: y0in(neqin), rtol, atol(*)
   real*8, intent(in) :: rpar(*)
   real*8, intent(in) :: id_vec(neqin)
   real*8 :: constr_vec(neqin)
   real*8 :: t0, yptmp(neqin)
   integer :: nthreads, iatol, ier
   integer :: i
   integer :: meth, itmeth
   integer :: myid

   dypdr = 0
{dypdrcalc}

   dvacdy = 0
{dvacdycalc}

   iatol = 2
   constr_vec = 1.d0

   y0 = y0in
   yp0 = 0
   yptmp = 0
   diff = id_vec
   mas = 0
   t0 = 0
   meth = 2  ! 1 = Adams (nonstiff), 2 = BDF (stiff)
   itmeth = 2  ! 1 = functional iteration, 2 = Newton iteration

   do i = 1, neq
      mas(i, i) = id_vec(i)
   enddo

!   ! Calculate yp
   call fidaresfun(0.d0, y0, yptmp, yp0, ipar, rpar, ier)

   ! initialize Sundials
   call fnvinits(2, neq, ier)
   ! allocate memory
   call fidamalloc(t0, y0, yp0, iatol, rtol, atol, iout, rout, ipar, rpar, ier)
   ! set maximum number of steps (default = 500)
   call fidasetiin('MAX_NSTEPS', 50000, ier)
   ! set algebraic variables
   call fidasetvin('ID_VEC', id_vec, ier)
   ! set constraints (all yi >= 0.)
   call fidasetvin('CONSTR_VEC', constr_vec, ier)

  call FSUNDenseMatInit(2, neq, neq, ier)
  call FSUNLAPACKDENSEINIT(2, ier)
  call FIDALSINIT(ier)

end subroutine initialize

subroutine find_steady_state(neqin, nrates, dt, maxiter, epsilon, t1, u1, du1, r1)

   use solve_ida, only: y0, yp0, iout, rout, rates, dypdr

   implicit none

   integer, intent(in) :: neqin, nrates, maxiter
   real*8, intent(in) :: dt, epsilon

   real*8 :: rpar(1)
   integer :: ipar(1)

   real*8, intent(out) :: t1, u1(neqin), du1(neqin), r1(nrates)

   real*8 :: tout, epsilon2
   real*8 :: dutmp(neqin), du0(neqin)
   integer :: itask, ier
   integer :: i

   logical :: converged = .FALSE.

   epsilon2 = epsilon**2
   i = 0
   itask = 1
   tout = 0.0d0
   u1 = y0
   du1 = yp0
   t1 = 0.d0
   du0 = 0.d0

   call fidacalcic(1, dt, ier)

   do while (.not. converged)
      if (tout - t1 < dt * 0.01) then
         tout = tout + dt
      end if

      call fidasolve(tout, t1, u1, du1, itask, ier)

      i = i + 1

      call fidaresfun(tout, u1, du0, dutmp, ipar, rpar, ier)

      if (maxval(dutmp**2) < epsilon2) then
         converged = .TRUE.
      end if
      if (i >= maxiter) then
         print *, "ODE NOT CONVERGED!"
         exit
      end if
   end do
   
   call ratecalc({neq}, u1)
   r1 = rates

end subroutine find_steady_state

subroutine solve(neqin, nrates, nt, tfinal, t1, u1, du1, r1)

   use solve_ida, only: y0, yp0, iout, rout, rates

   implicit none

   integer, intent(in) :: neqin, nt, nrates
   real*8, intent(in) :: tfinal

   real*8 :: rpar(1)
   integer :: ipar(1)

   real*8, intent(out) :: t1(nt)
   real*8, intent(out) :: u1(neqin, nt), du1(neqin, nt)
   real*8, intent(out) :: r1(nrates, nt)

   real*8 :: dt, tout
   integer :: itask, ier
   integer :: i

   itask = 1
   dt = tfinal / (nt - 1)
   tout = 0.0d0
   u1 = 0
   du1 = 0
   u1(:, 1) = y0
   du1(:, 1) = yp0
   t1(1) = 0.d0
   call ratecalc({neq}, u1(:, 1))
   r1(:, 1) = rates

   do i = 2, nt
      tout = tout + dt
      do while (tout - t1(i) > dt * 0.01)
         call fidasolve(tout, t1(i), u1(:, i), du1(:, i), itask, ier)
      end do
      r1(:, i) = rates
   end do

end subroutine solve

subroutine finalize

   implicit none

   call fidafree

end subroutine finalize

subroutine fidaresfun(tres, yin, ypin, res, ipar, rpar, reserr)

   use solve_ida, only: neq, diff, dypdr, rates

   implicit none

   integer, intent(in) :: ipar(*)
   integer, intent(out) :: reserr
   real*8, intent(in) :: tres, rpar(*)
   real*8, intent(in) :: yin(neq), ypin(neq)
   real*8, intent(out) :: res(neq)
   real*8 :: y(neq)

   integer :: i

   reserr = 0

   y = yin
   res = 0


   do i = 1, neq
      if (y(i) < -1d-10)  then
!         y(i) = 0.d0
         reserr = 1
      endif
   enddo

   call ratecalc({neq}, y)

   res = matmul(dypdr, rates) - diff * ypin
   
end subroutine fidaresfun

subroutine ratecalc(neqin, yin)

   use solve_ida, only: rates

   implicit none

   integer, intent(in) :: neqin
   real*8, intent(in) :: yin(neqin)
   real*8 :: y(neqin)
   real*8 :: vac({nvac})

   integer :: i

   y = yin


   vac = 0
{vaccalc}

   do i = 1, {nvac}
      if (vac(i) < -1d-10) then
         vac(i) = 0.d0
      endif
   enddo

   rates = 0
{ratecalc}

end subroutine ratecalc

   """


pyf_template = """!    -*- f90 -*-
! Note: the context of this file is case sensitive.

python module {modname} ! in
    interface  ! in :{modname}
        module solve_ida ! in :{modname}:{modname}.f90
            integer dimension(50) :: iout
            real*8 dimension(50) :: rout
            real*8 dimension({neq}) :: y0
            real*8 dimension({neq}) :: yp0
            real*8 dimension({neq}) :: diff
            real*8 dimension({neq},{neq}) :: mas
            real*8 dimension({neq},{neq}) :: jac
            real*8 dimension({nrates}) :: rates
            real*8 dimension({neq},{nrates}) :: dypdr
            integer dimension({nvac},{neq}) :: dvacdy
            integer, optional :: neq={neq}
        end module solve_ida
        subroutine initialize(neqin,y0in,rtol,atol,ipar,rpar,id_vec) ! in :{modname}:{modname}.f90
            use solve_ida, only: neq,iout,rout,y0,yp0,mas,diff,dypdr,dvacdy
            integer, optional,intent(in),check(len(y0in)>=neqin),depend(y0in) :: neqin=len(y0in)
            real*8 dimension(neqin),intent(in) :: y0in
            real*8 intent(in) :: rtol
            real*8 dimension(*),intent(in) :: atol
            integer dimension(*),intent(in) :: ipar
            real*8 dimension(*),intent(in) :: rpar
            real*8 dimension(neqin),intent(in),depend(neqin) :: id_vec
        end subroutine initialize
        subroutine find_steady_state(neqin,nrates,dt,maxiter,epsilon,t1,u1,du1,r1) ! in :{modname}:{modname}.f90
            use solve_ida, only: y0,yp0,iout,rout,rates,dypdr
            integer intent(in) :: neqin
            integer intent(in) :: nrates
            real*8 intent(in) :: dt
            integer intent(in) :: maxiter
            real*8 intent(in) :: epsilon
            real*8 intent(out) :: t1
            real*8 intent(out),dimension(neqin),depend(neqin) :: u1
            real*8 intent(out),dimension(neqin),depend(neqin) :: du1
            real*8 intent(out),dimension(nrates),depend(nrates) :: r1
        end subroutine find_steady_state
        subroutine solve(neqin,nrates,nt,tfinal,t1,u1,du1,r1) ! in :{modname}:{modname}.f90
            use solve_ida, only: y0,yp0,iout,rout,rates
            integer intent(in) :: neqin
            integer intent(in) :: nrates
            integer intent(in) :: nt
            real*8 intent(in) :: tfinal
            real*8 intent(out),dimension(nt),depend(nt) :: t1
            real*8 intent(out),dimension(neqin,nt),depend(neqin,nt) :: u1
            real*8 intent(out),dimension(neqin,nt),depend(neqin,nt) :: du1
            real*8 intent(out),dimension(nrates,nt),depend(nrates,nt) :: r1
        end subroutine solve
        subroutine finalize ! in :{modname}:{modname}.f90
        end subroutine finalize
    end interface
end python module {modname}

! This file was auto-generated with f2py (version:2).
! See http://cens.ioc.ee/projects/f2py2e/"""
