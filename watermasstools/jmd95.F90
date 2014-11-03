!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module eos
     implicit none
     private
! !PUBLIC MEMBER FUNCTIONS:

   public :: state,        &
             state_surface 
!EOP
!BOC

    integer, parameter ::               &
       char_len       = 256                    ,&
       char_len_long  = 512                    ,&
       log_kind       = kind(.true.)           ,&
       int_kind       = kind(1)                ,&
       i4             = selected_int_kind(6)   ,&
       i8             = selected_int_kind(13)  ,&
       r4             = selected_real_kind(6)  ,&
       r8             = selected_real_kind(13)

    real (r8), parameter :: &
       c0     =    0.0_r8   ,&
       c1     =    1.0_r8   ,&
       c2     =    2.0_r8   ,&
       c3     =    3.0_r8   ,&
       c4     =    4.0_r8   ,&
       c5     =    5.0_r8   ,&
       c8     =    8.0_r8   ,&
       c10    =   10.0_r8   ,&
       c16    =   16.0_r8   ,&
       c1000  = 1000.0_r8   ,&
       c10000 =10000.0_r8   ,&
       c1p5   =    1.5_r8   ,&
       p33    = c1/c3       ,&
       p5     = 0.500_r8    ,&
       p25    = 0.250_r8    ,&
       p125   = 0.125_r8    ,&
       p001   = 0.001_r8    ,&
       eps    = 1.0e-10_r8  ,&
       eps2   = 1.0e-20_r8  ,&
       bignum = 1.0e+30_r8


!-----------------------------------------------------------------------
!
!  valid ranges and pressure as function of depth
!
!-----------------------------------------------------------------------

   real (r8), parameter :: & 
      tmin = -5.0_r8, &
      tmax = 50.0_r8, &! valid temperature range for level k
      smin = 0.0_r8, &
      smax = 50.0_r8 ! valid salinity    range for level k

!-----------------------------------------------------------------------
!
!  choices for eos type and valid range checks
!
!-----------------------------------------------------------------------

   integer (i4), parameter :: &
      state_type_jmcd       = 1,    &! integer ids for state choice
      state_type_mwjf       = 2,    &
      state_type_polynomial = 3,    &
      state_type_linear     = 4

   integer (i4), parameter ::    &
      state_itype           = 1 ! input state type chosen

!-----------------------------------------------------------------------
!
!  UNESCO EOS constants and JMcD bulk modulus constants
!
!-----------------------------------------------------------------------

   !*** for density of fresh water (standard UNESCO)

   real (r8), parameter ::              &
      unt0 =   999.842594_r8,           &
      unt1 =  6.793952e-2_r8,           &
      unt2 = -9.095290e-3_r8,           &
      unt3 =  1.001685e-4_r8,           &
      unt4 = -1.120083e-6_r8,           &
      unt5 =  6.536332e-9_r8

   !*** for dependence of surface density on salinity (UNESCO)

   real (r8), parameter ::              &
      uns1t0 =  0.824493_r8 ,           &
      uns1t1 = -4.0899e-3_r8,           &
      uns1t2 =  7.6438e-5_r8,           &
      uns1t3 = -8.2467e-7_r8,           &
      uns1t4 =  5.3875e-9_r8,           &
      unsqt0 = -5.72466e-3_r8,          &
      unsqt1 =  1.0227e-4_r8,           &
      unsqt2 = -1.6546e-6_r8,           &
      uns2t0 =  4.8314e-4_r8

   !*** from Table A1 of Jackett and McDougall

  real (r8), parameter ::              &
      bup0s0t0 =  1.965933e+4_r8,       &
      bup0s0t1 =  1.444304e+2_r8,       &
      bup0s0t2 = -1.706103_r8   ,       &
      bup0s0t3 =  9.648704e-3_r8,       &
      bup0s0t4 = -4.190253e-5_r8

   real (r8), parameter ::              &
      bup0s1t0 =  5.284855e+1_r8,       &
      bup0s1t1 = -3.101089e-1_r8,       &
      bup0s1t2 =  6.283263e-3_r8,       &
      bup0s1t3 = -5.084188e-5_r8

   real (r8), parameter ::              &
      bup0sqt0 =  3.886640e-1_r8,       &
      bup0sqt1 =  9.085835e-3_r8,       &
      bup0sqt2 = -4.619924e-4_r8

   real (r8), parameter ::              &
      bup1s0t0 =  3.186519_r8   ,       &
      bup1s0t1 =  2.212276e-2_r8,       &
      bup1s0t2 = -2.984642e-4_r8,       &
      bup1s0t3 =  1.956415e-6_r8 

   real (r8), parameter ::              &
      bup1s1t0 =  6.704388e-3_r8,       &
      bup1s1t1 = -1.847318e-4_r8,       &
      bup1s1t2 =  2.059331e-7_r8,       &
      bup1sqt0 =  1.480266e-4_r8 

   real (r8), parameter ::              &
      bup2s0t0 =  2.102898e-4_r8,       &
      bup2s0t1 = -1.202016e-5_r8,       &
      bup2s0t2 =  1.394680e-7_r8,       &
      bup2s1t0 = -2.040237e-6_r8,       &
      bup2s1t1 =  6.128773e-8_r8,       &
      bup2s1t2 =  6.207323e-10_r8



!BOP
! !MODULE: state_mod
!
! !DESCRIPTION:
!  This module contains routines necessary for computing the density 
!  from model temperature and salinity using an equation of state.
!
!  The module supports four forms of EOS:
!  \begin{enumerate}
!     \item The UNESCO equation of state computed using the 
!           potential-temperature-based bulk modulus from Jackett and 
!           McDougall, JTECH, Vol.12, pp 381-389, April, 1995.
!     \item The faster and more accurate alternative to the UNESCO eos
!           of McDougall, Wright, Jackett and Feistel (hereafter 
!           MWJF, 2001 submission to JTECH).
!     \item a polynomial fit to the full UNESCO EOS
!     \item a simple linear EOS based on constant expansion coeffs
!  \end{enumerate}
!
! !REVISION HISTORY:
!  SVN:$Id: state_mod.F90 26114 2010-12-17 20:29:34Z njn01 $

! !USES:

!!!! Can't use any of the other modules

!    use kinds_mod
!    use blocks
!    use distribution
!    use domain
!    use constants
!    use grid
!    use io
!    use broadcast
!    use time_management
!    use exit_mod
!    use shr_sys_mod

!#ifdef CCSMCOUPLED
!   !*** ccsm
!   use shr_vmath_mod
!#endif




!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: state
! !INTERFACE:

! subroutine state(k, kk, TEMPK, SALTK, this_block, &
!                         RHOOUT, RHOFULL, DRHODT, DRHODS)
 subroutine state(PRESREF, TEMPK, SALTK,   &
                        nx_block, ny_block,      &
                        RHOOUT, DRHODT, DRHODS)

! !DESCRIPTION:
!  Returns the density of water at level k from equation of state
!  $\rho = \rho(d,\theta,S)$ where $d$ is depth, $\theta$ is
!  potential temperature, and $S$ is salinity. the density can be
!  returned as a perturbation (RHOOUT) or as the full density
!  (RHOFULL). Note that only the polynomial EOS choice will return
!  a perturbation density; in other cases the full density is returned
!  regardless of which argument is requested.
!
!  This routine also computes derivatives of density with respect
!  to temperature and salinity at level k from equation of state
!  if requested (ie the optional arguments are present).
!
!  If $k = kk$ are equal the density for level k is returned.
!  If $k \neq kk$ the density returned is that for a parcel
!  adiabatically displaced from level k to level kk.
!
! !REVISION HISTORY:
!  same as module

   implicit none
!   private
!   save

! !INPUT PARAMETERS:

!   INTEGER*8, intent(in) :: &
!      k,                    &! depth level index
!      kk                     ! level to which water is adiabatically 
                            ! displaced

   INTEGER*8, intent(in) :: &
       nx_block,               &
       ny_block

   REAL*8, intent(in) :: PRESREF
!f2py intent(in) PRESREF 
      !PRES,                    &! pressure 
                     ! pressure to which water is adiabatically 
                            ! displaced

   REAL*8, dimension(nx_block,ny_block), intent(in) :: & 
      TEMPK,             &! temperature at level k
      SALTK               ! salinity    at level k
!f2py intent(in) TEMPK
!f2py intent(in) SALTK

!   type (block), intent(in) :: &
!      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

!   REAL*8, dimension(nx_block,ny_block), optional, intent(out) :: & 
   REAL*8, dimension(nx_block,ny_block), intent(out) :: & 
      RHOOUT,  &! perturbation density of water
      DRHODT,  &! derivative of density with respect to temperature
      DRHODS    ! derivative of density with respect to salinity
!f2py intent(out) RHOOUT
!f2py intent(out) DRHODT
!f2py intent(out) DRHODT


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------



   INTEGER*8 :: &
      ib,ie,jb,je,       &! extent of physical domain
      bid,               &! local block index
      out_of_range        ! counter for out-of-range T,S values

   REAL*8, dimension(nx_block,ny_block) :: &
      TQ,SQ,             &! adjusted T,S
      BULK_MOD,          &! Bulk modulus
      RHO_S,             &! density at the surface
      DRDT0,             &! d(density)/d(temperature), for surface
      DRDS0,             &! d(density)/d(salinity   ), for surface
      DKDT,              &! d(bulk modulus)/d(pot. temp.)
      DKDS,              &! d(bulk modulus)/d(salinity  )
      SQR,DENOMK,        &! work arrays
      WORK1, WORK2, WORK3, WORK4, T2

   REAL*8 :: p, p2 ! temporary pressure scalars

!-----------------------------------------------------------------------
!
!  first check for valid range if requested
!
!-----------------------------------------------------------------------

!   bid = this_block%local_id

   !select case (state_range_iopt)
   !case (state_range_ignore)

      !*** prevent problems with garbage on land points or ghost cells
      
      TQ = min(TEMPK, c1000)
      TQ = max(TQ,   -c1000)
      SQ = min(SALTK, c1000)
      SQ = max(SALTK, c0)


      TQ = min(TEMPK,tmax)
      TQ = max(TQ,tmin)
      SQ = min(SALTK,smax)
      SQ = max(SQ,smin)

!   end select

!-----------------------------------------------------------------------
!
!  now compute density or expansion coefficients
!
!-----------------------------------------------------------------------

   !select case (state_itype)

!-----------------------------------------------------------------------
!
!  Jackett and McDougall EOS
!
!-----------------------------------------------------------------------

      ! formulas need pressure in bars, not decibars
      p = 0.10 * PRESREF
      p2  = p*p

      !SQ  = c1000*SQ
      SQR = sqrt(SQ)
      T2  = TQ*TQ


      !***
      !*** first calculate surface (p=0) values from UNESCO eqns.
      !***

      WORK1 = uns1t0 + uns1t1*TQ + & 
             (uns1t2 + uns1t3*TQ + uns1t4*T2)*T2
      WORK2 = SQR*(unsqt0 + unsqt1*TQ + unsqt2*T2)

      RHO_S = unt1*TQ + (unt2 + unt3*TQ + (unt4 + unt5*TQ)*T2)*T2 &
                      + (uns2t0*SQ + WORK1 + WORK2)*SQ


      !***
      !*** now calculate bulk modulus at pressure p from 
      !*** Jackett and McDougall formula
      !***

      WORK3 = bup0s1t0 + bup0s1t1*TQ +                    &
             (bup0s1t2 + bup0s1t3*TQ)*T2 +                &
              p *(bup1s1t0 + bup1s1t1*TQ + bup1s1t2*T2) + &
              p2*(bup2s1t0 + bup2s1t1*TQ + bup2s1t2*T2)
      WORK4 = SQR*(bup0sqt0 + bup0sqt1*TQ + bup0sqt2*T2 + &
                   bup1sqt0*p)

      BULK_MOD  = bup0s0t0 + bup0s0t1*TQ +                    &
                  (bup0s0t2 + bup0s0t3*TQ + bup0s0t4*T2)*T2 + &
                  p *(bup1s0t0 + bup1s0t1*TQ +                &
                     (bup1s0t2 + bup1s0t3*TQ)*T2) +           &
                  p2*(bup2s0t0 + bup2s0t1*TQ + bup2s0t2*T2) + &
                  SQ*(WORK3 + WORK4)

      DENOMK = c1/(BULK_MOD - p)

      !***
      !*** now calculate required fields
      !***

      RHOOUT = ((unt0 + RHO_S)*BULK_MOD*DENOMK)

      !if (present(DRHODT)) then
      DRDT0 =  unt1 + c2*unt2*TQ +                      &
               (c3*unt3 + c4*unt4*TQ + c5*unt5*T2)*T2 + &
               (uns1t1 + c2*uns1t2*TQ +                 &
                (c3*uns1t3 + c4*uns1t4*TQ)*T2 +         &
                (unsqt1 + c2*unsqt2*TQ)*SQR )*SQ

      DKDT  = bup0s0t1 + c2*bup0s0t2*TQ +                       &
              (c3*bup0s0t3 + c4*bup0s0t4*TQ)*T2 +               &
              p *(bup1s0t1 + c2*bup1s0t2*TQ + c3*bup1s0t3*T2) + &
              p2*(bup2s0t1 + c2*bup2s0t2*TQ) +                  &
              SQ*(bup0s1t1 + c2*bup0s1t2*TQ + c3*bup0s1t3*T2 +  &
                  p  *(bup1s1t1 + c2*bup1s1t2*TQ) +             &
                  p2 *(bup2s1t1 + c2*bup2s1t2*TQ) +             &
                  SQR*(bup0sqt1 + c2*bup0sqt2*TQ))

      DRHODT = (DENOMK*(DRDT0*BULK_MOD -   &
                        p*(unt0+RHO_S)*DKDT*DENOMK))

      DRDS0  = c2*uns2t0*SQ + WORK1 + c1p5*WORK2
      DKDS = WORK3 + c1p5*WORK4

      DRHODS = DENOMK*(DRDS0*BULK_MOD -                    &
                      p*(unt0+RHO_S)*DKDS*DENOMK)

 end subroutine state

 subroutine state_surface(TEMPK, SALTK,         &
                        nx_block, ny_block,     &
                        RHOOUT, DRHODT, DRHODS)
                        
  implicit none

   INTEGER*8, intent(in) ::    &
       nx_block,               &
       ny_block

   REAL*8, dimension(nx_block,ny_block), intent(in) :: & 
      TEMPK,             & ! temperature at level k
      SALTK                ! salinity    at level k  
!f2py intent(in) TEMPK
!f2py intent(in) SALTK

   REAL*8, dimension(nx_block,ny_block), intent(out) :: & 
      RHOOUT,  &! perturbation density of water
      DRHODT,  &! derivative of density with respect to temperature
      DRHODS    ! derivative of density with respect to salinity
!f2py intent(out) RHOOUT
!f2py intent(out) DRHODT
!f2py intent(out) DRHODS

   REAL*8, dimension(nx_block,ny_block) :: &
      TQ,SQ,             &! adjusted T,S
      RHO_S,             &! density at the surface
      SQR,               &! work arrays
      WORK1, WORK2, T2

    REAL*8 :: p, p2 ! temporary pressure scalars

    !*** prevent problems with garbage on land points or ghost cells

    TQ = min(TEMPK,tmax)
    TQ = max(TQ,tmin)
    SQ = min(SALTK,smax)
    SQ = max(SQ,smin)

    p = c0
    p2  = p*p

    SQR = sqrt(SQ)
    T2  = TQ*TQ

    !***
    !*** calculate surface (p=0) values from UNESCO eqns.
    !***

    WORK1 = uns1t0 + uns1t1*TQ + & 
           (uns1t2 + uns1t3*TQ + uns1t4*T2)*T2
    WORK2 = SQR*(unsqt0 + unsqt1*TQ + unsqt2*T2)

    RHO_S = unt1*TQ + (unt2 + unt3*TQ + (unt4 + unt5*TQ)*T2)*T2 &
                    + (uns2t0*SQ + WORK1 + WORK2)*SQ

    RHOOUT = unt0 + RHO_S

    DRHODT =  unt1 + c2*unt2*TQ +                      &
             (c3*unt3 + c4*unt4*TQ + c5*unt5*T2)*T2 + &
             (uns1t1 + c2*uns1t2*TQ +                 &
              (c3*uns1t3 + c4*uns1t4*TQ)*T2 +         &
              (unsqt1 + c2*unsqt2*TQ)*SQR )*SQ

    DRHODS  = c2*uns2t0*SQ + WORK1 + c1p5*WORK2

 end subroutine state_surface
 
 end module eos
