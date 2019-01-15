#include "fabm_driver.h"


module cserra_npztemp

use fabm_types

   implicit none 

!  default: all is private.
   private

   integer, parameter :: nb_sp = 20


! !PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_cserra_npztemp
   
!     Variable identifiers
      type (type_state_variable_id)        :: id_n, id_g1(nb_sp), id_g2(nb_sp), id_b, id_pom, id_doc 
      type (type_dependency_id)            :: id_par, id_temp
      type (type_global_dependency_id)     :: id_yearday
      !type (type_diagnostic_variable_id)   :: id_inv_L(nb_sp), id_inv_B(nb_sp)
      type (type_diagnostic_variable_id)   :: id_g1_tot, id_g2_tot, id_NPP
      !type (type_diagnostic_variable_id)   :: id_VL1(nb_sp), id_VL2(nb_sp), id_VB1(nb_sp), id_VB2(nb_sp)
      type (type_diagnostic_variable_id)   :: id_VL1_tot, id_VL2_tot, id_VB1_tot, id_VB2_tot
      !type (type_diagnostic_variable_id)   :: id_JL1(nb_sp), id_JL2(nb_sp)
      type (type_diagnostic_variable_id)   :: id_JL1_tot, id_JN1_tot, id_JF1_tot, id_JR1_tot, id_Jleak1_tot, id_JNleak1_tot, id_mu1_tot, id_loss_g1
      type (type_diagnostic_variable_id)   :: id_JL2_tot, id_JN2_tot, id_JF2_tot, id_JR2_tot, id_Jleak2_tot, id_JNleak2_tot, id_mu2_tot, id_loss_g2, id_death2
      !type (type_diagnostic_variable_id)   :: id_J_L2, id_J_N2, id_JF2, id_J_R2, id_mu2
      type (type_diagnostic_variable_id)   :: id_J_Rb, id_mub, id_Jdoc, id_Jpom, id_JNb
      type (type_diagnostic_variable_id)   :: id_mean_invL1, id_mean_invL2, id_mean_invB1, id_mean_invB2, id_loss_b
      type (type_diagnostic_variable_id)   :: id_A_N1(nb_sp),id_A_N2(nb_sp), id_A_L1(nb_sp), id_A_L2(nb_sp), id_A_F1(nb_sp), id_A_F2(nb_sp), id_M1(nb_sp), id_M2(nb_sp)
      type (type_diagnostic_variable_id)   :: id_Llvl1, id_Nlvl1, id_Flvl1, id_Llvl2, id_Nlvl2, id_Flvl2, id_Cgains1, id_Ngains1, id_Cgains2, id_Ngains2
      type (type_diagnostic_variable_id)   :: id_JCleak1(nb_sp), id_JCleak2(nb_sp), id_JNleak1(nb_sp), id_JNleak2(nb_sp)
      type (type_horizontal_diagnostic_variable_id) :: id_export_pom, id_rate_n

!     Model parameters
    real(rk) :: n,g1(nb_sp),g2(nb_sp),b,pom,doc
    real(rk) :: Vst1, Vst2, c_L, a_L, c_N, c_F, k, m, mu_mort, C_cn, loss, lossP2
    real(rk) :: V_L1(nb_sp), V_B1(nb_sp), Vtot1, A_L1, A_N1, A_F1, M1, J_L1(nb_sp), J_N1(nb_sp), J_F1(nb_sp)
    real(rk) :: J_R1(nb_sp), Jtot1, mu1(nb_sp), Jleak1(nb_sp), JNleak1(nb_sp)
    real(rk) :: V_L2(nb_sp), V_B2(nb_sp), Vtot2, A_L2, A_N2, A_F2, M2, J_L2(nb_sp), J_N2(nb_sp), J_F2(nb_sp)
    real(rk) :: J_R2(nb_sp), Jtot2, mu2(nb_sp), Jleak2(nb_sp), JNleak2(nb_sp), death2
    real(rk) :: c_DOC, Vstb, c_POM, Vtotb, JNleakB, lossB, JleakB
    real(rk) :: A_doc, A_Nb, A_POM, Mb, J_DOCb, J_POMb, J_Nb, J_Rb, J_Rb_d, Jtotb, mub, mu_POM, pom_frac   
    real(rk) :: JL1_tot, JN1_tot, JF1_tot, JR1_tot, Jleak1_tot, JNleak1_tot
    real(rk) :: JL2_tot, JN2_tot, JF2_tot, JR2_tot, Jleak2_tot, JNleak2_tot
    real(rk) :: g1_tot, g2_tot, NPP
    real(rk) :: J_im, J_em, f, im
    real(rk) :: max_inv, delta_inv
    real(rk) :: kc, J_LL1
    real(rk) :: mean_invL1, mean_invL2, mean_invB1, mean_invB2
    real(rk) :: invL1_each(nb_sp), invL2_each(nb_sp), invB1_each(nb_sp), invB2_each(nb_sp)  
    real(rk) :: Q10_N, Q10_m, Q10_r, Q10_mumort, Tref, TrefN
    real(rk) :: A_N1_d(nb_sp), A_L1_d(nb_sp), A_F1_d(nb_sp), A_N2_d(nb_sp), A_L2_d(nb_sp), A_F2_d(nb_sp), M1_d(nb_sp), M2_d(nb_sp)
    integer :: nb_sp

   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: get_light_extinction
   end type


!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the model namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_cserra_npztemp), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   integer             :: i,j
   character(LEN=nb_sp)    :: index

!EOP


real(rk), parameter :: Vst1 = 1E-5_rk 
real(rk), parameter :: Vst2 = 1E-2_rk
real(rk), parameter :: Vstb = 1E-8_rk
real(rk), parameter :: a_L = 0.43_rk !10.0_rk
real(rk), parameter :: c_L = 0.01_rk
real(rk), parameter :: c_N = 3.75E-5_rk !0.000162_rk
real(rk), parameter :: c_F = 0.01_rk
real(rk), parameter :: k = 0.05_rk !0.5_rk
real(rk), parameter :: m = 4.0_rk !2.3_rk
real(rk), parameter :: c_DOC = 8.8E-6_rk !2E-4_rk !7.6E-11_rk
real(rk), parameter :: c_POM = 80_rk !1E-11_rk
real(rk), parameter :: mu_mort = 0.005_rk !0.01_rk
real(rk), parameter :: C_cn = 5.68_rk
real(rk), parameter :: loss = 0.02_rk !0.005_rk
real(rk), parameter :: lossB = 0.02_rk !0.005_rk
real(rk), parameter :: lossP2= 0.02_rk !0.005_rk
real(rk), parameter :: Vtotb=2E-8_rk
real(rk), parameter :: f=0.05_rk
real(rk), parameter :: im=1E-7_rk
real(rk), parameter :: kc=8.8E-7_rk
real(rk), parameter :: max_inv=0.9_rk
real(rk), parameter :: Q10_N=1.5_rk
real(rk), parameter :: Q10_m=2.0_rk
real(rk), parameter :: Q10_r=2.0_rk
real(rk), parameter :: Q10_mumort=2.0_rk
real(rk), parameter :: Tref=18.0_rk
real(rk), parameter :: TrefN=18.0_rk
real(rk), parameter :: mu_POM=4.0_rk
real(rk), parameter :: pom_frac=0.3_rk

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day and are converted here to values per second.

    call self%get_parameter(self%Vst1,   'Vst1', 'mugC cell^-1 ',  'Structural biomass og G1', default=Vst1)
    call self%get_parameter(self%Vst2,   'Vst2', 'mugC cell^-1 ',  'Structural biomass og G2', default=Vst2)
    call self%get_parameter(self%Vstb,   'Vstb', 'mugC cell^-1 ',  'Structural biomass og B', default=Vstb)
    call self%get_parameter(self%a_L,    'a_L', 'mugC d^{-1} (W m^{-2})^{-1} (mugC^{2/3})^{-1}','Affinity per investment in photosynthesis', default=a_L, scale_factor=d_per_s) 
    call self%get_parameter(self%c_L,    'c_L',  'mugC d^(-1) (W m^{-2})^{-1} (mugC^{2/3})^{-1}',  'Maximum light affinity', default=c_L, scale_factor=d_per_s) 
    call self%get_parameter(self%c_N,    'c_N',  'L d^{-1} (mugC^{1/3})^{-1}',  'Maximum nutrient affinity', default=c_N, scale_factor=d_per_s)
    call self%get_parameter(self%c_F,    'c_F',  'L d^{-1} (mugC)^{-1}',  'Maximum food affinity', default=c_F, scale_factor=d_per_s)
    call self%get_parameter(self%k,      'k',    'd^{-1}',  'Metabolic cost', default=k, scale_factor=d_per_s)
    call self%get_parameter(self%m,      'm',    'd^{-1}',  'Biosynthetic rate', default=m, scale_factor=d_per_s)
    call self%get_parameter(self%c_DOC,  'c_DOC',      'L d^{-1} (mugC)^{-1}',  'Maximum DOC affinity', default=c_DOC, scale_factor=d_per_s)
    call self%get_parameter(self%c_POM,   'c_POM',      'L d^{-1} (mugC)^{-1}',  'Maximum POM affinity', default=c_POM)
    call self%get_parameter(self%mu_mort, 'mu_mort',    'mugC^{1/4} d^{-1}',  'Loss rate to predators',  default=mu_mort, scale_factor=d_per_s)
    call self%get_parameter(self%C_cn,    'C_cn',       'mugC mugN^{-1}',  'C:N ratio', default=C_cn)
    call self%get_parameter(self%loss,    'loss',       'd^{-1}',  'Loss rate G1', default=loss, scale_factor=d_per_s)
    call self%get_parameter(self%lossB,   'lossB',      'd^{-1}',  'Loss rate B',  default=lossB, scale_factor=d_per_s)
    call self%get_parameter(self%lossP2,  'lossP2',     'd^{-1}',  'Loss rate G2', default=lossP2, scale_factor=d_per_s)
    call self%get_parameter(self%Vtotb,   'Vtotb',      'mugC cell^-1',  'Total structural biomass of bacteria', default=VtotB)
    call self%get_parameter(self%f,       'f',     'd^{-1}',  'im/emmigration rate', default=f, scale_factor=d_per_s)
    call self%get_parameter(self%im,      'im',      'mugC cells l^-1',  'background concentration', default=im)
    call self%get_parameter(self%kc,      'kc',      'm² mugC⁻1',  'attenuation coefficient from plankton', default=kc)
    call self%get_parameter(self%max_inv, 'max_inv',      '-',  'values of max inv for light', default=max_inv)
    call self%get_parameter(self%Q10_N, 'Q10_N',      '-',  'Q10 of N uptake', default=Q10_N)
    call self%get_parameter(self%Q10_m, 'Q10_m',      '-',  'Q10 of biosynthesis', default=Q10_m)
    call self%get_parameter(self%Q10_r, 'Q10_r',      '-',  'Q10 of respiration', default=Q10_r)
    call self%get_parameter(self%Q10_mumort, 'Q10_mumort',      '-',  'Q10 of mortality HTL', default=Q10_mumort)
    call self%get_parameter(self%Tref, 'Tref',      'degrees celsius',  'refference temperature for Q10 params.', default=Tref)
    call self%get_parameter(self%TrefN, 'TrefN',      'degrees celsius',  'refference temperature for Q10 N doc params.', default=TrefN)
    call self%get_parameter(self%mu_POM, 'mu_POM',      'days⁻1',  'maximum pom hydrolysis rate', default=mu_POM, scale_factor=d_per_s)
    call self%get_parameter(self%pom_frac, 'pom_frac',      '-',  'fraction of dead material going to DOC', default=pom_frac)

!Register state variables
   call self%register_state_variable(self%id_n,'N','mugN /l','nutrient concentration',     1.0_rk,minimum=0.0_rk)
   call self%register_state_variable(self%id_b,'B','mugC cells /l','B concentration',      0.02_rk,minimum=0.0_rk)
   call self%register_state_variable(self%id_doc,'DOC','mugC /l','DOC concentration',      0.02_rk,minimum=0.0_rk)
   call self%register_state_variable(self%id_pom,'POM','mugC /l','POM concentration',      1E-5_rk,minimum=0.0_rk,vertical_movement=-1.0_rk/86400.0_rk)

do i=1,nb_sp

 write (index,'(i0)') i

 !call self%register_state_variable(self%id_g1(i),'G1'//trim(index),'mug C /l','G1 biomass'//trim(index),1E-5_rk,minimum=0.0_rk,vertical_movement=-1_rk/86400_rk)
 call self%register_state_variable(self%id_g2(i),'G2'//trim(index),'mug C /l','G2 biomass'//trim(index),2.0_rk,minimum=0.0_rk,vertical_movement=-0.0_rk/86400.0_rk)
 call self%register_state_variable(self%id_g1(i),'G1'//trim(index),'mug C /l','G1 biomass'//trim(index),5.0_rk,minimum=0.0_rk,vertical_movement=-0.0_rk/86400.0_rk)
      !call self%register_diagnostic_variable(self%id_VL1(i),'VL1'//trim(index),'-','Photosynthetic biomass of G1 sp. '//trim(index))
      !call self%register_diagnostic_variable(self%id_VL2(i),'VL2'//trim(index),'-','Photosynthetic biomass of G2 sp. '//trim(index))
      !call self%register_diagnostic_variable(self%id_VB1(i),'VB1'//trim(index),'-','Biosynthetic biomass of G1 sp. '//trim(index))
      !call self%register_diagnostic_variable(self%id_VB2(i),'VB2'//trim(index),'-','Biosynthetic biomass of G2 sp. '//trim(index))

      !call self%register_diagnostic_variable(self%id_JL1(i),'J_L1'//trim(index),'-','light up for species '//trim(index))
      !call self%register_diagnostic_variable(self%id_JL2(i),'J_L2'//trim(index),'-','light up species '//trim(index))
call self%register_diagnostic_variable(self%id_A_N1(i),'A_N1_'//trim(index),'-','Affinity for nutrients G1_'//trim(index))
call self%register_diagnostic_variable(self%id_A_N2(i),'A_N2_'//trim(index),'-','Affinity for nutrients G2_'//trim(index))
call self%register_diagnostic_variable(self%id_A_L1(i),'A_L1_'//trim(index),'-','Affinity for light G1_'//trim(index))
call self%register_diagnostic_variable(self%id_A_L2(i),'A_L2_'//trim(index),'-','Affinity for light G2_'//trim(index))
call self%register_diagnostic_variable(self%id_A_F1(i),'A_F1_'//trim(index),'-','Affinity for food G1_'//trim(index))
call self%register_diagnostic_variable(self%id_A_F2(i),'A_F2_'//trim(index),'-','Affinity for food G2_'//trim(index))
call self%register_diagnostic_variable(self%id_M1(i),'M1_'//trim(index),'-','max gortwh rate G1_'//trim(index))
call self%register_diagnostic_variable(self%id_M2(i),'M2_'//trim(index),'-','max gortwh rate G2_'//trim(index))

call self%register_diagnostic_variable(self%id_JCleak1(i),'JCleak1_'//trim(index),'days⁻1','JCleak G1_'//trim(index))
call self%register_diagnostic_variable(self%id_JCleak2(i),'JCleak2_'//trim(index),'days⁻1','JCleak G2_'//trim(index))
call self%register_diagnostic_variable(self%id_JNleak1(i),'JNleak1_'//trim(index),'days⁻1','JNleak G1_'//trim(index))
call self%register_diagnostic_variable(self%id_JNleak2(i),'JNleak2_'//trim(index),'days⁻1','JNleak G2_'//trim(index))

end do


!Register diagnostics
call self%register_diagnostic_variable(self%id_g1_tot,'G1_tot','mugC cells l⁻1','Total biomass of G1')
call self%register_diagnostic_variable(self%id_g2_tot,'G2_tot','mugC cells l⁻1','Total biomass of G2')
call self%register_diagnostic_variable(self%id_NPP,'NPP','mugC l⁻1 days⁻1','total carbon fixed by all')
call self%register_diagnostic_variable(self%id_J_Rb,'J_Rb','mugC days⁻1','total respiration from bacteria')
call self%register_diagnostic_variable(self%id_mub,'mub','days⁻1','growth rate bacteria')
call self%register_diagnostic_variable(self%id_Jdoc,'Jdoc','mugC days⁻1','DOC uptake from bacteria')
call self%register_diagnostic_variable(self%id_Jpom,'Jpom','mugC days⁻1','POM degradation')
call self%register_diagnostic_variable(self%id_JNb,'JNb','mugC days⁻1','Nuptake from bacteria')


call self%register_diagnostic_variable(self%id_mean_invL1,'mean_invL1','-','Mean photo investment for G1')
call self%register_diagnostic_variable(self%id_mean_invL2,'mean_invL2','-','Mean photo investment for G2')
call self%register_diagnostic_variable(self%id_mean_invB1,'mean_invB1','-','Mean synth investment for G1')
call self%register_diagnostic_variable(self%id_mean_invB2,'mean_invB2','-','Mean synth investment for G2')

call self%register_diagnostic_variable(self%id_JL1_tot,'JL1_tot','mugC l⁻1 days⁻1','total carbon fixed by G1')
call self%register_diagnostic_variable(self%id_JN1_tot,'JN1_tot','mugN l⁻1 days⁻1','total nutrients uptake rate by G1')
call self%register_diagnostic_variable(self%id_JF1_tot,'JF1_tot','mugC l⁻1 days⁻1','total food uptake rate by G1')
call self%register_diagnostic_variable(self%id_JR1_tot,'JR1_tot','mugC l⁻1 days⁻1','total respiration by G1')
call self%register_diagnostic_variable(self%id_Jleak1_tot,'Jleak1_tot','mugC l⁻1 days⁻1','total C leak rate by G1')
call self%register_diagnostic_variable(self%id_JNleak1_tot,'JNleak1_tot','mugC l⁻1 days⁻1','total N leak rate by G1')
call self%register_diagnostic_variable(self%id_mu1_tot,'mu1_tot','mugC l⁻1 days⁻1','total growth G1')

call self%register_diagnostic_variable(self%id_JL2_tot,'JL2_tot','mugC l⁻1 days⁻1','total carbon fixed by G2')
call self%register_diagnostic_variable(self%id_JN2_tot,'JN2_tot','mugN l⁻1 days⁻1','total nutrients uptake rate by G2')
call self%register_diagnostic_variable(self%id_JF2_tot,'JF2_tot','mugC l⁻1 days⁻1','total food uptake rate by G2')
call self%register_diagnostic_variable(self%id_JR2_tot,'JR2_tot','mugC l⁻1 days⁻1','total respiration by G2')
call self%register_diagnostic_variable(self%id_Jleak2_tot,'Jleak2_tot','mugC l⁻1 days⁻1','total C leak rate by G2')
call self%register_diagnostic_variable(self%id_JNleak2_tot,'JNleak2_tot','mugC l⁻1 days⁻1','total N leak rate by G2')
call self%register_diagnostic_variable(self%id_mu2_tot,'mu2_tot','mugC l⁻1 days⁻1','total growth G2')

call self%register_diagnostic_variable(self%id_VL1_tot,'VL1_tot','mugC l⁻1 ','Total photosynthetic biomass of G1')
call self%register_diagnostic_variable(self%id_VL2_tot,'VL2_tot','mugC l⁻1 ','Total photosynthetic biomass of G2')
call self%register_diagnostic_variable(self%id_VB1_tot,'VB1_tot','mugC l⁻1 ','Total Biosynthetic biomass of G1')
call self%register_diagnostic_variable(self%id_VB2_tot,'VB2_tot','mugC l⁻1 ','Total Biosynthetic biomass of G2')

call self%register_diagnostic_variable(self%id_loss_g1,'loss_g1','mug l⁻1 d⁻1 ','background mortality G1')
call self%register_diagnostic_variable(self%id_loss_g2,'loss_g2','mug l⁻1 d⁻1 ','background mortality G2')
call self%register_diagnostic_variable(self%id_loss_b,'loss_b','mug l⁻1 d⁻1 ','background mortality B')
call self%register_diagnostic_variable(self%id_death2,'death2','mug l⁻1 d⁻1 ','loss to HTL')

call self%register_diagnostic_variable(self%id_Llvl1,'Llvl1','-','Llvl1')
call self%register_diagnostic_variable(self%id_Nlvl1,'Nlvl1','-','Nlvl1')
call self%register_diagnostic_variable(self%id_Flvl1,'Flvl1','-','Flvl1')
call self%register_diagnostic_variable(self%id_Llvl2,'Llvl2','-','Llvl2')
call self%register_diagnostic_variable(self%id_Nlvl2,'Nlvl2','-','Nlvl2')
call self%register_diagnostic_variable(self%id_Flvl2,'Flvl2','-','Flvl2')

call self%register_diagnostic_variable(self%id_Cgains1,'Cgains1','-','Cgains1')
call self%register_diagnostic_variable(self%id_Ngains1,'Ngains1','-','Ngains1')
call self%register_diagnostic_variable(self%id_Cgains2,'Cgains2','-','Cgains2')
call self%register_diagnostic_variable(self%id_Ngains2,'Ngains2','-','Ngains2')



call self%register_horizontal_diagnostic_variable(self%id_export_pom,'export','-','export')
call self%register_horizontal_diagnostic_variable(self%id_rate_n,'rate_n','-','bottom exchange of n')

   !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_export_pom,rate_pom)

!Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_yearday, standard_variables%number_of_days_since_start_of_the_year) 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 



end subroutine initialize
!EOC

   
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!EOC
 !IROUTINE: Right hand sides of mixosize model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_cserra_npztemp), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
    integer i,j
    real(rk) :: n,g1(nb_sp),g2(nb_sp),b,pom,doc
    real(rk) :: dndt,dg1dt(nb_sp), dg2dt(nb_sp), dbdt,ddocdt,dpomdt, dg1dtmin(nb_sp), dg2dtmin(nb_sp)
	real(rk) :: light, temp
	real(rk) :: Vst1, Vst2, c_L, a_L, c_N, c_F, k, m, mu_mort, C_cn, loss, lossP2
	real(rk) :: V_L1(nb_sp), V_B1(nb_sp), Vtot1, A_L1, A_N1, A_F1, M1, J_L1(nb_sp), J_N1(nb_sp), J_F1(nb_sp), J_R1(nb_sp)
        real(rk) :: Jtot1, mu1(nb_sp), Jleak1(nb_sp), JNleak1(nb_sp)
	real(rk) :: V_L2(nb_sp), V_B2(nb_sp), Vtot2, A_L2, A_N2, A_F2, M2, J_L2(nb_sp), J_N2(nb_sp), J_F2(nb_sp), J_R2(nb_sp) 
        real(rk) :: Jtot2, mu2(nb_sp), Jleak2(nb_sp), JNleak2(nb_sp), death2
	real(rk) :: c_DOC, Vstb, c_POM, Vtotb, JNleakB, lossB, JleakB, mu_POM
	real(rk) :: A_doc, A_Nb, A_POM, Mb, J_DOCb, J_POMb, J_Nb, J_Rb, J_Rb_d, Jtotb, mub, lossb_tot, pom_frac   
	real(rk) :: max_inv, delta_inv
    real(rk) :: J_im, J_em, f, im
    real(rk) :: JL1_tot, JN1_tot, JF1_tot, JR1_tot, Jleak1_tot, JNleak1_tot, mu1_tot, loss1_tot, death2_tot, Llvl1, Nlvl1, Flvl1, Cgains1, Ngains1
    real(rk) :: JL2_tot, JN2_tot, JF2_tot, JR2_tot, Jleak2_tot, JNleak2_tot, mu2_tot, loss2_tot, Llvl2, Nlvl2, Flvl2, Cgains2, Ngains2
    real(rk) :: g1_tot, g2_tot, NPP
    real(rk) :: mean_invL1, mean_invL2, mean_invB1, mean_invB2
    real(rk) :: invL1_each(nb_sp), invL2_each(nb_sp), invB1_each(nb_sp), invB2_each(nb_sp)
    real(rk) :: inv_L(nb_sp), inv_B(nb_sp)
    real(rk) :: A_N1_d(nb_sp), A_L1_d(nb_sp), A_F1_d(nb_sp), A_N2_d(nb_sp), A_L2_d(nb_sp), A_F2_d(nb_sp), M1_d(nb_sp), M2_d(nb_sp)
    real(rk) :: Q10_N, Q10_m, Q10_r, Q10_mumort, Tref, TrefN
    real(rk) :: c_Nt, c_DOCt, mt, kt, mu_mortt, c_POMt, mu_POMt
    real(rk) :: VL1_tot, VL2_tot, VB1_tot, VB2_tot
    real(rk) :: fdiff, time, imN

   real(rk), parameter        :: secs_pr_day = 86400.0_rk
   real(rk), parameter        :: Pi = 3.1415927_rk



! Enter spatial loops (if any)
   _LOOP_BEGIN_
   
! Retrieve current (local) state variable values.
    _GET_(self%id_n,n)  ! nutrient

  do i=1,nb_sp
    _GET_(self%id_g1(i),g1(i))
    _GET_(self%id_g2(i),g2(i))
  end do
    
    _GET_(self%id_b,b)
    _GET_(self%id_pom,pom)
    _GET_(self%id_doc,doc)

! Retrieve current environmental conditions.
   _GET_(self%id_par,light)             ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)
   _GET_GLOBAL_(self%id_yearday,time)  
 
!initializing the total rates
JL1_tot= 0
JN1_tot= 0
JF1_tot= 0
JR1_tot= 0
Jleak1_tot= 0
JNleak1_tot=0
loss1_tot=0
JL2_tot= 0
JN2_tot= 0
JF2_tot= 0
JR2_tot= 0
Jleak2_tot= 0
JNleak2_tot= 0
loss2_tot=0
VL1_tot=0
VL2_tot=0
VB1_tot=0
VB2_tot=0
mu1_tot=0
mu2_tot=0
death2_tot=0

 Cgains1=0
 Ngains1=0
 Cgains2=0
 Ngains2=0

Llvl1=0
Nlvl1=0
Flvl1=0
Llvl2=0
Nlvl2=0
Flvl2=0

delta_inv=self%max_inv/nb_sp

inv_L(1)=0.01_rk 
inv_B(1)=1.0_rk-inv_L(1)
do j=2,nb_sp
 inv_L(j)=inv_L(j-1)+delta_inv
 inv_B(j)=1.0_rk-inv_L(j)
end do


temp=temp+3.0_rk

 c_Nt=self%c_N*self%Q10_N**((temp-self%TrefN)/10.0_rk)
 c_DOCt=self%c_DOC*self%Q10_N**((temp-self%TrefN)/10.0_rk)
 mt=self%m*self%Q10_m**((temp-self%Tref)/10.0_rk)
 kt=self%k*self%Q10_r**((temp-self%Tref)/10.0_rk)
 mu_mortt=self%mu_mort*self%Q10_mumort**((temp-self%Tref)/10.0_rk)
 c_POMt=self%c_POM*self%Q10_N**((temp-self%Tref)/10.0_rk)
 mu_POMt=self%mu_POM*self%Q10_r**((temp-self%Tref)/10.0_rk)

do i=1,nb_sp

!small plankton ----------------------------------------------------------------------------------
    V_L1(i)=inv_L(i)*self%Vst1 !photosynthetic biomass
    V_B1(i)=inv_B(i)*self%Vst1 !biosynthetic biomass
    Vtot1=self%Vst1+V_L1(i)+V_B1(i) !total biomass
    A_L1=(self%c_L*self%a_L*inv_L(i)*self%Vst1*self%Vst1**(2.0_rk/3.0_rk))&
    	/((self%a_L*inv_L(i)*self%Vst1)+self%c_L*self%Vst1**(2.0_rk/3.0_rk)) !affinity for light
    A_N1=c_Nt*self%Vst1**(1.0_rk/3.0_rk)
    A_F1=self%c_F*self%Vst1
    M1=mt*inv_B(i)*self%Vst1
    J_L1(i)=M1*(A_L1*light)/((A_L1*light)+M1) !uptake of light
    J_N1(i)=M1*(self%C_cn*A_N1*n)/((self%C_cn*A_N1*n)+M1) !uptake of nutrients
    J_F1(i)=M1*(A_F1*b)/((A_F1*b)+M1)
    J_R1(i)=self%Vst1*kt*(1.0_rk+inv_L(i)+inv_B(i))
    Jtot1=min(J_L1(i)+J_F1(i)-J_R1(i), J_N1(i)+J_F1(i))
    mu1(i)=Jtot1/(self%Vst1*(1.0_rk+inv_L(i)+inv_B(i)))
    !Jleak1(i)=max(0.0_rk,0.3_rk*(J_L1(i)-J_R1(i)-J_N1(i)))
    Jleak1(i)=max(0.0_rk,(J_L1(i)-J_R1(i)-J_N1(i)))
    JNleak1(i)=max(0.0_rk,-(J_L1(i)-J_R1(i)-J_N1(i)))

JL1_tot=(JL1_tot + J_L1(i)*g1(i)/Vtot1)
JN1_tot= (JN1_tot + J_N1(i)*g1(i)/Vtot1)
JF1_tot= (JF1_tot + J_F1(i)*g1(i)/Vtot1)
JR1_tot=(JR1_tot + J_R1(i)*g1(i)/Vtot1)
Jleak1_tot= (Jleak1_tot + Jleak1(i)*g1(i)/Vtot1)
JNleak1_tot=(JNleak1_tot + JNleak1(i)*g1(i)/Vtot1)
mu1_tot=(mu1_tot + mu1(i)*g1(i))
loss1_tot=loss1_tot+self%loss*g1(i)**2.0_rk

 Cgains1=Cgains1+((J_L1(i)+J_F1(i)-J_R1(i))/((self%Vst1*(1.0_rk+inv_L(i)+inv_B(i)))))*g1(i)
 Ngains1=Ngains1+((J_N1(i)+J_F1(i))/(self%Vst1*(1.0_rk+inv_L(i)+inv_B(i))))*g1(i)

Llvl1=Llvl1+((A_L1*light)/((A_L1*light)+M1))*G1(i)
Nlvl1=Nlvl1+((self%C_cn*A_N1*n)/((self%C_cn*A_N1*n)+M1))*G1(i)
Flvl1=Flvl1+((A_F1*b)/((A_F1*b)+M1))*G1(i)


!large plankton ------------------------------------------------------------------------------------
    V_L2(i)=inv_L(i)*self%Vst2 !photosynthetic biomass
    V_B2(i)=inv_B(i)*self%Vst2 !biosynthetic biomass
    Vtot2=self%Vst2+V_L2(i)+V_B2(i) !total biomass
    A_L2=(self%c_L*self%a_L*inv_L(i)*self%Vst2*self%Vst2**(2.0_rk/3.0_rk))&
    	/((self%a_L*inv_L(i)*self%Vst2)+self%c_L*self%Vst2**(2.0_rk/3.0_rk)) !affinity for light
    A_N2=c_Nt*self%Vst2**(1.0_rk/3.0_rk)
    A_F2=self%c_F*self%Vst2
    M2=mt*inv_B(i)*self%Vst2
    J_L2(i)=M2*(A_L2*light)/((A_L2*light)+M2) !uptake of light
    J_N2(i)=M2*(self%C_cn*A_N2*n)/((self%C_cn*A_N2*n)+M2) !uptake of nutrients
    J_F2(i)=M2*(A_F2*sum(g1))/((A_F2*sum(g1))+M2) !hacemos la suma de G1 porque les da igual sus investments
    J_R2(i)=self%Vst2*kt*(1.0_rk+inv_L(i)+inv_B(i))
    Jtot2=min(J_L2(i)+J_F2(i)-J_R2(i), J_N2(i)+J_F2(i))
    mu2(i)=Jtot2/(self%Vst2*(1.0_rk+inv_L(i)+inv_B(i)))
    death2=mu_mortt*(self%Vst2**(-1.0_rk/4.0_rk))
    !Jleak2(i)=max(0.0_rk,0.3_rk*(J_L2(i)-J_R2(i)-J_N2(i)))
    Jleak2(i)=max(0.0_rk,(J_L2(i)-J_R2(i)-J_N2(i)))
    JNleak2(i)=max(0.0_rk,-(J_L2(i)-J_R2(i)-J_N2(i)))

JL2_tot=(JL2_tot + J_L2(i)*g2(i)/Vtot2)
JN2_tot= (JN2_tot + J_N2(i)*g2(i)/Vtot2)
JF2_tot= (JF2_tot + J_F2(i)*g2(i)/Vtot2)
JR2_tot=(JR2_tot + J_R2(i)*g2(i)/Vtot2)
Jleak2_tot= (Jleak2_tot + Jleak2(i)*g2(i)/Vtot2)
JNleak2_tot=(JNleak2_tot + JNleak2(i)*g2(i)/Vtot2)
mu2_tot=(mu2_tot + mu2(i)*g2(i))
loss2_tot=loss2_tot+self%loss*g2(i)**2.0_rk
death2_tot=(death2_tot + death2*g2(i))

Llvl2=Llvl2+((A_L2*light)/((A_L2*light)+M2))*G2(i)
Nlvl2=Nlvl2+((self%C_cn*A_N2*n)/((self%C_cn*A_N2*n)+M2))*G2(i)
Flvl2=Flvl2+((A_F2*b)/((A_F2*b)+M2))*G2(i)

 Cgains2=Cgains2+((J_L2(i)+J_F2(i)-J_R2(i))/((self%Vst2*(1.0_rk+inv_L(i)+inv_B(i)))))*g2(i)
 Ngains2=Ngains2+((J_N2(i)+J_F2(i))/(self%Vst2*(1.0_rk+inv_L(i)+inv_B(i))))*g2(i)

VL1_tot=VL1_tot+V_L1(i)*g1(i)/Vtot1
VL2_tot=VL2_tot+V_L2(i)*g2(i)/Vtot2
VB1_tot=VB1_tot+V_B1(i)*g1(i)/Vtot1
VB2_tot=VB2_tot+V_B2(i)*g2(i)/Vtot2

A_N1_d(i)=A_N1
A_L1_d(i)=A_L1
A_F1_d(i)=A_F1
A_N2_d(i)=A_N2
A_L2_d(i)=A_L2
A_F2_d(i)=A_F2
M1_d(i)=M1
M2_d(i)=M2

_SET_DIAGNOSTIC_(self%id_A_N1(i),A_N1_d(i))
_SET_DIAGNOSTIC_(self%id_A_L1(i),A_L1_d(i))
_SET_DIAGNOSTIC_(self%id_A_F1(i),A_F1_d(i))
_SET_DIAGNOSTIC_(self%id_A_N2(i),A_N2_d(i))
_SET_DIAGNOSTIC_(self%id_A_L2(i),A_L2_d(i))
_SET_DIAGNOSTIC_(self%id_A_F2(i),A_F2_d(i))
_SET_DIAGNOSTIC_(self%id_M1(i),M1_d(i))
_SET_DIAGNOSTIC_(self%id_M2(i),M2_d(i))

_SET_DIAGNOSTIC_(self%id_JCleak1(i),Jleak1(i))
_SET_DIAGNOSTIC_(self%id_JCleak2(i),Jleak2(i))
_SET_DIAGNOSTIC_(self%id_JNleak1(i),JNleak1(i))
_SET_DIAGNOSTIC_(self%id_JNleak2(i),JNleak2(i))

end do

!total biomass of generalists

g1_tot=sum(g1)
g2_tot=sum(g2)

_SET_DIAGNOSTIC_(self%id_g1_tot,g1_tot)
_SET_DIAGNOSTIC_(self%id_g2_tot,g2_tot)

_SET_DIAGNOSTIC_(self%id_VL1_tot,VL1_tot)
_SET_DIAGNOSTIC_(self%id_VL2_tot,VL2_tot)
_SET_DIAGNOSTIC_(self%id_VB1_tot,VB1_tot)
_SET_DIAGNOSTIC_(self%id_VB2_tot,VB2_tot)

NPP=JL1_tot+JL2_tot-JR1_tot-JR2_tot

_SET_DIAGNOSTIC_(self%id_NPP,NPP)

_SET_DIAGNOSTIC_(self%id_JL1_tot,JL1_tot)
_SET_DIAGNOSTIC_(self%id_JN1_tot,JN1_tot)
_SET_DIAGNOSTIC_(self%id_JF1_tot,JF1_tot)
_SET_DIAGNOSTIC_(self%id_JR1_tot,JR1_tot)
_SET_DIAGNOSTIC_(self%id_Jleak1_tot,Jleak1_tot)
_SET_DIAGNOSTIC_(self%id_JNleak1_tot,JNleak1_tot)
_SET_DIAGNOSTIC_(self%id_mu1_tot,mu1_tot)

_SET_DIAGNOSTIC_(self%id_JL2_tot,JL2_tot)
_SET_DIAGNOSTIC_(self%id_JN2_tot,JN2_tot)
_SET_DIAGNOSTIC_(self%id_JF2_tot,JF2_tot)
_SET_DIAGNOSTIC_(self%id_JR2_tot,JR2_tot)
_SET_DIAGNOSTIC_(self%id_Jleak2_tot,Jleak2_tot)
_SET_DIAGNOSTIC_(self%id_JNleak2_tot,JNleak2_tot)
_SET_DIAGNOSTIC_(self%id_mu2_tot,mu2_tot)

_SET_DIAGNOSTIC_(self%id_loss_g1,loss1_tot)
_SET_DIAGNOSTIC_(self%id_loss_g2,loss2_tot)
_SET_DIAGNOSTIC_(self%id_death2,death2_tot)

_SET_DIAGNOSTIC_(self%id_Llvl1,Llvl1)
_SET_DIAGNOSTIC_(self%id_Nlvl1,Nlvl1)
_SET_DIAGNOSTIC_(self%id_Flvl1,Flvl1)
_SET_DIAGNOSTIC_(self%id_Llvl2,Llvl2)
_SET_DIAGNOSTIC_(self%id_Nlvl2,Nlvl2)
_SET_DIAGNOSTIC_(self%id_Flvl2,Flvl2)

_SET_DIAGNOSTIC_(self%id_Cgains1,Cgains1)
_SET_DIAGNOSTIC_(self%id_Ngains1,Ngains1)
_SET_DIAGNOSTIC_(self%id_Cgains2,Cgains2)
_SET_DIAGNOSTIC_(self%id_Ngains2,Ngains2)



!mean investments
do i=1,nb_sp
invL1_each(i)=inv_L(i)*g1(i)
invL2_each(i)=inv_L(i)*g2(i)
invB1_each(i)=inv_B(i)*g1(i)
invB2_each(i)=inv_B(i)*g2(i)
end do

mean_invL1=sum(invL1_each)/sum(g1)
mean_invL2=sum(invL2_each)/sum(g2)
mean_invB1=sum(invB1_each)/sum(g1)
mean_invB2=sum(invB2_each)/sum(g2)

 _SET_DIAGNOSTIC_(self%id_mean_invL1,mean_invL1)
 _SET_DIAGNOSTIC_(self%id_mean_invL2,mean_invL2)
 _SET_DIAGNOSTIC_(self%id_mean_invB1,mean_invB1)
 _SET_DIAGNOSTIC_(self%id_mean_invB2,mean_invB2)


!Bacteria ------------------------------------------------------------------------------------------

    A_Nb=c_Nt*self%Vstb**(1.0_rk/3.0_rk)
    Mb=mt*self%Vstb*0.5_rk
    J_DOCb=(Mb*c_doct*doc)/(c_doct*doc+Mb) !uptake of DOC
    !J_POMb=(mu_POMt*pom)/(c_POMt+pom)
    J_POMb=(mu_POMt*(mu_POMt/c_POMt)*pom)/((mu_POMt/c_POMt)*pom+mu_POMt)
    !J_POMb=mu_POMt*pom
    J_Nb=Mb*(self%C_cn*A_Nb*n)/((self%C_cn*A_Nb*n)+Mb) !uptake of DIN
    J_Rb=self%Vtotb*kt
    Jtotb= min(J_DOCb - J_Rb , J_Nb)
    mub=Jtotb/self%Vtotb
    JleakB=max(0.0_rk,J_DOCb-J_Rb-J_Nb)
    JNleakB=max(0.0_rk,-(J_DOCb-J_Rb-J_Nb))


_SET_DIAGNOSTIC_(self%id_J_Rb,J_Rb)
_SET_DIAGNOSTIC_(self%id_mub,mub)
_SET_DIAGNOSTIC_(self%id_Jdoc,J_DOCb)
_SET_DIAGNOSTIC_(self%id_Jpom,J_POMb)
_SET_DIAGNOSTIC_(self%id_JNb,J_Nb)



!Immigration--------------------------------------------------------------------------------------------------
     
!if (light>1.0_rk) then
    !fdiff=0.00045_rk/86400.0_rk
!else 
    !fdiff=0.05_rk/86400.0_rk
!end if

    !fdiff= (0.5002 + 0.4998 *(cos(time*2*pi/365)))/86400.0_rk    
    !J_im=self%f*self%im/nb_sp !i/emmigration flows
    !J_em=self%f



!  Use size based model
!  Calculate change of everything

do i=1,nb_sp
   dg1dt(i)=mu1(i)*(g1(i)) - self%loss*(g1(i)**2.0_rk) - (g1(i)/sum(g1))*JF2_tot !+J_im -J_em*g1(i)
   dg1dtmin(i)=0.0_rk
   !dg1dtmin(i)=mu1(i)*g1(i)

   dg2dt(i)=mu2(i)*(g2(i)) - death2*(g2(i)) - self%lossP2*(g2(i)**2.0_rk) !+J_im -J_em*g2(i)
   dg2dtmin(i)=0.0_rk
   !dg2dtmin(i)=mu2(i)*g2(i)

end do


    dbdt=mub*(b) - self%lossB*(b**2.0_rk) - JF1_tot !+ J_im - J_em*b

    dndt=(-JN1_tot - JN2_tot - J_Nb*b/self%Vtotb &
         + JNleak1_tot + JNleakB*b/self%Vtotb  + JNleak2_tot + J_POMb*b + self%pom_frac*loss1_tot &
         + self%pom_frac*loss2_tot)/(self%C_cn) !+ J_Rb*b/self%Vtotb

    dpomdt=(1.0_rk-self%pom_frac)*loss1_tot + (1.0_rk-self%pom_frac)*loss2_tot - J_POMb*b

    ddocdt=Jleak1_tot + Jleak2_tot - J_DOCb*b/self%Vtotb + J_POMb*b + JleakB*b/self%Vtotb &
    + self%pom_frac*loss1_tot + self%pom_frac*loss2_tot


lossb_tot=self%lossB*(b**2.0_rk)
_SET_DIAGNOSTIC_(self%id_loss_b,lossb_tot)

!  Set temporal derivatives,  subtracted by the secs_pr_day


do i=1,nb_sp
!if(g1(i).lt.self%im)then
!_SET_ODE_(self%id_g1(i),dg1dtmin(i))   
!else
_SET_ODE_(self%id_g1(i),dg1dt(i))   
!end if

!if(g2(i).lt.self%im)then
!_SET_ODE_(self%id_g2(i),dg2dtmin(i))   
!else
_SET_ODE_(self%id_g2(i),dg2dt(i))
!end if

end do

!if(b.lt.self%im)then
!_SET_ODE_(self%id_b,dg1dtmin(1))  
!else
_SET_ODE_(self%id_b,dbdt)
!end if

!_SET_ODE_(self%id_b,dbdt)
_SET_ODE_(self%id_n,dndt)
_SET_ODE_(self%id_pom,dpomdt)
_SET_ODE_(self%id_doc,ddocdt)


   ! Leave spatial loops (if any)
   _LOOP_END_


   end subroutine do
!EOC

subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_cserra_npztemp),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_

   real(rk) :: n, pom, g1(nb_sp), g2(nb_sp)
   real(rk) :: rate_pom, rate_n, rate_g1(nb_sp), rate_g2(nb_sp)
   real(rk), parameter :: n_bottom=140.0_rk
   real(rk), parameter :: vel=1_rk/86400.0_rk !era 0.5
   integer  :: i
  
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_n,n)
   _GET_(self%id_pom,pom)

!do i=1,nb_sp

    !_GET_(self%id_g1(i),g1(i))
    !_GET_(self%id_g2(i),g2(i))

     !   rate_g1(i)=g1(i)*(-0.0_rk/86400.0_rk)
      !  rate_g2(i)=g2(i)*(-0.0_rk/86400.0_rk)

  ! _SET_BOTTOM_EXCHANGE_(self%id_g1(i),rate_g1(i))
  ! _SET_BOTTOM_EXCHANGE_(self%id_g2(i),rate_g2(i))

!end do

       rate_pom=pom*(-10.0_rk/86400.0_rk)
       rate_n=vel*(n_bottom-n)
       !rate_n=vel*n_bottom

   _SET_BOTTOM_EXCHANGE_(self%id_n,rate_n)
   _SET_BOTTOM_EXCHANGE_(self%id_pom,rate_pom)

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_export_pom,rate_pom)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rate_n,rate_n)
   
   _HORIZONTAL_LOOP_END_
end subroutine do_bottom

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_cserra_npztemp), intent(in)     :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !LOCAL VARIABLES:
   real(rk) :: g1(nb_sp),g2(nb_sp)
   real(rk) :: kc, im
   integer  :: i
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.

do i=1,nb_sp
    _GET_(self%id_g1(i),g1(i))
    _GET_(self%id_g2(i),g2(i))
end do

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kc*(self%im+sum(g1)+sum(g2)))

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC

!-----------------------------------------------------------------------




! !-----------------------------------------------------------------------
! 
    end module cserra_npztemp
! 
! !-----------------------------------------------------------------------
! ! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
! !-----------------------------------------------------------------------
