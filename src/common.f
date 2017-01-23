c/************************************************
c      common.f for PBTREESPH ver.p12
c  21 July, 2000    produced by D.KAWATA
c ************************************************/
c70*******************************************************************
      include 'global.f'

c /* for Particle */
c /* Mass */      
      common /PDD/m_p(0:MN-1),ms0_p(0:MN-1),eps_p(0:MN-1)
c /* Virtual Position Xn+1,Xn */      
      common /PDD/x_p(0:MN-1),y_p(0:MN-1),z_p(0:MN-1)
c /* Correct position */      
      common /PDD/xc_p(0:MN-1),yc_p(0:MN-1),zc_p(0:MN-1)
c /* Velosity Vn-1/2,Vn+1/2 */
      common /PDD/vx_p(0:MN-1),vy_p(0:MN-1),vz_p(0:MN-1)
c /* Vn+1 */      
      common /PDD/vnx_p(0:MN-1),vny_p(0:MN-1),vnz_p(0:MN-1)
c /* Virtual Vn+1 */      
      common /PDD/vvnx_p(0:MN-1),vvny_p(0:MN-1),vvnz_p(0:MN-1)
c /* Internal Energy Un+1,Un */       
      common /PDD/u_p(0:MN-1),pu_p(0:MN-1)
c /* Smoothing Length hn+1,hn */      
      common /PDD/h_p(0:MN-1)
c /* Density rhon+1,rhon */       
      common /PDD/rho_p(0:MN-1)
c /* Presure Pn+1,Pn */      
      common /PDD/p_p(0:MN-1)
c /* Acceleration */       
      common /PDD/ax_p(0:MN-1),ay_p(0:MN-1),az_p(0:MN-1)
c /* Sound Velocity */      
      common /PDD/cs_p(0:MN-1)
c /* div V */
      common /PDD/div_v_p(0:MN-1)
c /* |rot V| */      
      common /PDD/arot_v_p(0:MN-1)
c /* Max Viscosity Coefficient */      
      common maxmyu_p(0:MN-1)
c /* dv/dt */      
      common /PDD/dvx_p(0:MN-1),dvy_p(0:MN-1),dvz_p(0:MN-1)
c /* du/dt_n,du/dt_n+1 */      
      common /PDD/pdu_p(0:MN-1),ndu_p(0:MN-1)
c /* Number of neigbor within 2xh */      
      common numnb_p(0:MN-1)
c /*** For Star ***/
c /* energy of supernova */      
      common /PDD/Gsn_p(0:MN-1)
c /* For Heavy Elements unit Solar Mass */
      common /PYDD/mzHe_p(0:MN-1),mzC_p(0:MN-1),mzN_p(0:MN-1),
     &  mzO_p(0:MN-1),mzNe_p(0:MN-1),mzMg_p(0:MN-1),
     &  mzSi_p(0:MN-1),mzFe_p(0:MN-1),mzZ_p(0:MN-1)
c /* For Total Ejected mass by SNe unit Solar Mass */
c      common /PYDD/tnsn_p(0:MN-1),tmej_p(0:MN-1),tmzHe_p(0:MN-1),
c     & tmzC_p(0:MN-1),tmzN_p(0:MN-1),tmzO_p(0:MN-1),tmzNe_p(0:MN-1),
c     & tmzMg_p(0:MN-1),tmzSi_p(0:MN-1),tmzFe_p(0:MN-1),
c     & tmzZ_p(0:MN-1)
c /*** For Individual time step ***/
c /* Individual time bin */      
      common /PTDD/dt_p(0:MN-1)
c /* system time bin */      
      common /PTDD/lt_p(0:MN-1)
c /* Virtual Time bin */      
      common /PTDD/vdt_p(0:MN-1)
c /* for cooling */
      common /PTDD/ram_p(0:MN-1)
c /* star formtion time */
      common /PTDD/ts_p(0:MN-1)
c /* flag for cold phase 0: cold */
      common flagc_p(0:MN-1)
      
c /***** For Dark Matter *****/
c /* mass */      
      common /DMDD/m_dm(0:MNDM-1),eps_dm(0:MNDM-1)
c /* Virtual Position */      
      common /DMDD/x_dm(0:MNDM-1),y_dm(0:MNDM-1),z_dm(0:MNDM-1)
c /* Correct Position */      
      common /DMDD/xc_dm(0:MNDM-1),yc_dm(0:MNDM-1),zc_dm(0:MNDM-1)
c /* Velocity V(n+1/2) */
      common /DMDD/vx_dm(0:MNDM-1),vy_dm(0:MNDM-1),vz_dm(0:MNDM-1)
c /* Virtual Velocity V(n+1/2) */      
      common /DMDD/vnx_dm(0:MNDM-1),vny_dm(0:MNDM-1),vnz_dm(0:MNDM-1)
c /* Acceration */      
      common /DMDD/dvx_dm(0:MNDM-1),dvy_dm(0:MNDM-1),dvz_dm(0:MNDM-1)
c /* for Individual time steep */
c /* Individual time bin */      
      common /DMDD/dt_dm(0:MNDM-1)
c /* time in system time step */      
      common /DMDD/lt_dm(0:MNDM-1)
c /* virtual time bin */      
      common /DMDD/vdt_dm(0:MNDM-1)


c /***** For Individual Time Step *****/
c /* list of active particles for gas and star */      
      common list_ap(0:MN-1)
c /* list of active particles for DM */      
      common list_adm(0:MNDM-1)


c /***** for Work ****/
      common list(0:MNDM-1),node(0:MNDM-1)
c /* for work */
      common /WD/tx(0:MNDM-1),ty(0:MNDM-1),tz(0:MNDM-1)
      common /WD/tcx(0:7),tcy(0:7),tcz(0:7)
 
c /***** Definition of Structure *****/
c /* TimeStep DTn DTn-1*/      
      common /TMD/TM_dt
c /* Hubble time */      
      common /TMD/TM_t0
c /* Total Time */      
      common /TMD/TM_tot
c /* System time bin */      
      common /TMD/TM_sdt
c /* Minimum local time bin */      
      common /TMD/TM_ldtmin
c /* local Time in system */      
      common /TMD/TM_lt
      common /TMD/TM_k,TM_kn
      
c /***** For OFPara *****/
c /* This structure have Number of output file. */
c /* Make file for display by step */
      common /OFD/OF_count,OF_scount
      common /OFD/OF_odt,OF_todt
      
c /***** For System Information *****/
c /*** 0: live halo, otherwise: fixed halo ***/
      common SI_nh(MNG),SI_nd(MNG)
c /* number of high resolution DM ***/
      common SI_ndm1
c /* 0 = include cooling, 1 = without */      
      common SI_flagc
c /* 0 = include starformation 1 = without */      
      common SI_flags
c /* 0 = include hydrodynamics otherwise = no */      
      common SI_flaghyd
c /* flag for making output file */
      common SI_flagout
c /* Cosmological Value */
c      common /SYSD/SI_h0,SI_omg0,SI_lam0,SI_a
c/* table for scale factor and time */
c      common /TABLED/at_ttb(0:NATTABLE),ta_ttb(0:NATTABLE)
c /* for memory of feedback history */
      common /SYSD/SI_esnIa,SI_esnII,SI_snzg,SI_snzs
      common /SYSD/SI_esng,SI_esns

c /***** For Yields lookup table *****/
      common /YTB/z_ytb(0:NYTZ+1),t_ytb(0:NYTT+1),
     & mej_ytb(0:NYTT+1,0:NYTZ+1),nsn_ytb(0:NYTT+1,0:NYTZ+1),
     & mzHe_ytb(0:NYTT+1,0:NYTZ+1),mzC_ytb(0:NYTT+1,0:NYTZ+1),
     & mzN_ytb(0:NYTT+1,0:NYTZ+1),mzO_ytb(0:NYTT+1,0:NYTZ+1),
     & mzNe_ytb(0:NYTT+1,0:NYTZ+1),mzMg_ytb(0:NYTT+1,0:NYTZ+1),
     & mzSi_ytb(0:NYTT+1,0:NYTZ+1), mzFe_ytb(0:NYTT+1,0:NYTZ+1),
     & mzZ_ytb(0:NYTT+1,0:NYTZ+1)
c /***** For Cooling Lookup Table *****/
      common /CTB/z_ctb(0:NCTZ-1),cr_ctb(0:NCTT,0:NCTZ-1)


c /**** for Random Numbrer ***/
      common idum      
