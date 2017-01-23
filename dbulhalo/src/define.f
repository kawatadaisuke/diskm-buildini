c /***************************************************
c   define.f  for bulhalo
c  22 Aug. 2006  written by D.Kawata
c ****************************************************/
c70*******************************************************************

      implicit none
      include 'mpif.h'

      double precision M_PI
      parameter (M_PI=3.141592654d0)
c *** MNS for sample disk / NCPU
c *** MN for _rzg
      integer MNR,MN,MNS
      parameter (MNR=120000,MN=1000000,MNS=40000000)
c *** max number of DM halo particles
      integer MNDM
      parameter (MNDM=100)
c *** number of bins for radial profile ***
      integer NDR
      parameter (NDR=10000)
c *** max number of 2D grid size
      integer MN2D
      parameter (MN2D=1001)
c *** minimum radius ***
      double precision FRINI0,FRINI0D,FRINI0DZ
c      parameter (FRINI0=1.0e-7)
      parameter (FRINI0=1.0e-6)
      parameter (FRINI0D=1.0e-6)
      parameter (FRINI0DZ=1.0e-6)


c *** maximum number of components ***
      integer MNC
      parameter (MNC=10)
c *****   For Tree   *****
c * Max Number of node in tree *
      integer MAXNODE,MNNODE
      parameter (MAXNODE=MNS*5)
      parameter (MNNODE=MAXNODE)

c /*** Normalization Unit ***/
      double precision MUSM,LUPC,LUKPC,LUCM,VUKMS,G,TMUGYR,DU
     &  ,ERU,TUK
c /* mass unit ( solar mass ) */      
      parameter (MUSM=1.0e12)
c /* (solarmass)/(gram) */      
c      parameter (MUG=1.99e45)
c /* length unit / pc ( 100kpc/pc ) */      
      parameter (LUPC=1.0e5)
c /* length unit(100kpc)/kpc */      
      parameter (LUKPC=1.0e2)
c /* (100 kpc)/(cm) */      
      parameter (LUCM=3.086e23)
c /* km/s */      
      parameter (VUKMS=207.4d0)
c /* Gravitatrional Constant */
      parameter (G=1.0d0)
c /* time ( Gyr ) */
      parameter (TMUGYR=0.471d0)
c * Density *      
      parameter (DU=6.77e-26)
c * cooling (energy) rate (erg s^-1 cm^-3) *      
      parameter (ERU=1.96e-27)
c * tempature ( K ) *      
      parameter (TUK=1.0e4)

      double precision KCGS,MYU,MP,K_MU,GCGS,SMG,KPCCM,HUB0,TPRHO
     &  ,XZSOL
c /* average molecular weight */      
      parameter (MYU=0.6d0)
c /* proton mass */      
      parameter (MP=1.67265e-24)
c /* Boltzmann constarnt */    
      parameter (KCGS=1.381e-16)
c /* G(cgs) */      
      parameter (GCGS=6.672e-8)
c /* k/m unit /( cm^2 s^-2 K^-1) */       
      parameter (K_MU=4.3e10)
c * (KCGS/(MYU*MP))/K_MU)))*TUK *
      parameter (TPRHO=((KCGS/(MYU*MP))/K_MU))
c /* (solarmass)/(gram) */      
      parameter (SMG=1.99e33)
c /* (100 kpc)/(cm) */      
      parameter (KPCCM=3.086e21)
c /* Hubble constant 9.78x10^9xh^-1 yr h = 1.0*/
      parameter (HUB0=100.0d0/(VUKMS*10.0d0))
c * Solar metallicity *
      parameter (XZSOL=0.019d0)
c *** specific heat ***
      double precision GAM
      parameter (GAM=5.0d0/3.0d0)
c *** smoothing parameter ***
      double precision ETAH
      parameter (ETAH=2.4d0)

c *** For Arithmetic Value ***
      double precision THIRD,INF,MININF,V34
      parameter (THIRD=1.0d0/3.0d0)
      parameter (V34=3.0d0/4.0d0)
      parameter (INF=1.0d6)
      parameter (MININF=1.0d-10)
      integer MAXFILENAME,MAXSTEP
      parameter (MAXFILENAME=25)
      parameter (MAXSTEP=1000000)
c *** Definition Parameters ***
c *** for Tree ***
      double precision MGROOT,THETA
c * Margin for Size of root tree *      
      parameter (MGROOT=1.001d0)
c * Torelance Parameter *
      parameter (THETA=0.8d0)
c *** for kernel tables ***
      integer NKTAB
      parameter (NKTAB=10000)
c *** for eos
      integer MNEOS
      parameter (MNEOS=20)
c *** density limit ***
      double precision RHO_LIM
c      parameter (RHO_LIM=1.0e-3)
c testing at ver.23
c      parameter (RHO_LIM=1.0e-5)
c for exp2d
c      parameter (RHO_LIM=1.0e-20)
      parameter (RHO_LIM=1.0e-15)

c *** for particles ***
c *** DM ***
      double precision x_dm,y_dm,z_dm,vx_dm,vy_dm,vz_dm,m_dm
      common x_dm(0:MNDM-1),y_dm(0:MNDM-1),z_dm(0:MNDM-1)
     &  ,vx_dm(0:MNDM-1),vy_dm(0:MNDM-1),vz_dm(0:MNDM-1)
     &  ,m_dm(0:MNDM-1)
c *** baryon particles    
c *** sample particles
      double precision x_p,y_p,z_p,vx_p,vy_p,vz_p,m_p,h_p
     & ,dvr_p
      common x_p(0:MNS-1),y_p(0:MNS-1),z_p(0:MNS-1),vx_p(0:MNS-1)
     &  ,m_p(0:MNS-1),h_p(0:MNS-1),dvr_p(0:MNS-1)

c *** radial mass profile ***
      double precision mrdmh,rdmh,rhodmh,dlrdmh,lridmh,lrodmh
      common mrdmh(0:MNR-1),rdmh(0:MNR-1),rhodmh(0:MNR-1)
     &  ,dlrdmh,lridmh,lrodmh

c for halo hot gas radial mass profile
      double precision mrhg,rmhg,rhorhg,temhg,uhg
      common mrhg(0:MNR-1),rmhg(0:MNR-1),rhorhg(0:MNR-1),temhg(0:MNR)
     &  ,uhg(0:MNR-1)

c *** exp2d gas surface density profile
      integer ndrexp2d
      common ndrexp2d
      double precision mrexp2d,rhorexp2d,rexp2d
      common mrexp2d(0:MNR),rhorexp2d(0:MNR),rexp2d(0:MNR)
      double precision riexp2d,roexp2d,drexp2d
      common riexp2d,roexp2d,drexp2d

c *** kernel look up table ****
      double precision ds_tb
      common /KTB/ds_tb
      double precision s_tb,w_tb,dwds_s_tb,dwdsc_tb
      common /KTB/s_tb(0:NKTAB),w_tb(0:NKTAB),dwds_s_tb(0:NKTAB)
     &  ,dwdsc_tb(0:NKTAB)
      double precision dphidr_r_tb,dphidh_tb
      common /KTB/dphidr_r_tb(0:NKTAB),dphidh_tb(0:NKTAB)
      double precision d2phidr2_tb,d3phidr3_tb   
      common /KTB/d2phidr2_tb(0:NKTAB),d3phidr3_tb(0:NKTAB)

      integer nmet_eos,nnh_eos
      common nmet_eos,nnh_eos
      double precision lnh_eos,ltm_eos,myu_eos,lmet_eos
      common lnh_eos(0:MNEOS-1),ltm_eos(0:MNEOS-1,0:MNEOS-1)
     &  ,myu_eos(0:MNEOS-1,0:MNEOS-1),lmet_eos(0:MNEOS-1)

c ***** for Tree used also in ddecb.F ****
c * level info *
c      common level_tr,llist_tr(0:MNB-1)
c * Number of nodes *
      integer num_tr
      common /PDI/num_tr
c * number of contained particle *      
       integer np_tr
      common /PDI/np_tr(0:MAXNODE-1)
c * name of Particle *      
      integer pn_tr
      common /PDI/pn_tr(0:MAXNODE-1)
c * length of side *       
      double precision l_tr
      common /TREED/l_tr(0:MAXNODE-1)
c * Coordinate of center *      
      double precision cx_tr,cy_tr,cz_tr
      common /TREED/cx_tr(0:MAXNODE-1),cy_tr(0:MAXNODE-1),
     & cz_tr(0:MAXNODE-1)
c * first child node *      
      integer daughter_tr
      common daughter_tr(0:MAXNODE-1)
c * next node *      
      integer next_tr
      common next_tr(0:MAXNODE-1)
c * center of mass *      
      double precision cmx_tr,cmy_tr,cmz_tr
      common /TREED/cmx_tr(0:MAXNODE-1),cmy_tr(0:MAXNODE-1),
     &  cmz_tr(0:MAXNODE-1)
c * total of mass *      
      double precision mass_tr
      common /TREED/mass_tr(0:MAXNODE-1)
c * distance between cm and center, maximum softening *      
      double precision delta_tr,hm_tr
      common /TREED/delta_tr(0:MAXNODE-1),hm_tr(0:MAXNODE-1)
c * for Multipole Expansion *
      double precision mx_tr,my_tr,mz_tr
      common /TREED/mx_tr(0:MAXNODE-1),my_tr(0:MAXNODE-1)
     &  ,mz_tr(0:MAXNODE-1)
      double precision mxx_tr,myy_tr,mzz_tr
      common /TREED/mxx_tr(0:MAXNODE-1),myy_tr(0:MAXNODE-1),
     &  mzz_tr(0:MAXNODE-1)
      double precision mxy_tr,myz_tr,mzx_tr
      common /TREED/mxy_tr(0:MAXNODE-1),myz_tr(0:MAXNODE-1),
     &  mzx_tr(0:MAXNODE-1)

c *** rz-grid ***
      integer id_rzg
      common id_rzg(0:MN2D-1,0:MN2D-1)
      double precision lri_rzg,lro_rzg,dlr_rzg,lzi_rzg,lzo_rzg
     &  ,dlz_rzg
      common lri_rzg,lro_rzg,dlr_rzg,lzi_rzg,lzo_rzg,dlz_rzg
      double precision x_rzg,y_rzg,z_rzg,rho_rzg,tm_rzg,p_rzg
     & ,myu_rzg
      common x_rzg(0:MN-1),y_rzg(0:MN-1),z_rzg(0:MN-1)
     & ,rho_rzg(0:MN-1),tm_rzg(0:MN-1),p_rzg(0:MN-1),myu_rzg(0:MN-1)
      double precision vsigr_rzg,vc2_rzg
     & ,vsigph_rzg,vph_rzg,vsigz_rzg
      common vsigr_rzg(0:MN-1),vc2_rzg(0:MN-1)
     & ,vsigph_rzg(0:MN-1),vph_rzg(0:MN-1),vsigz_rzg(0:MN-1)
      double precision fx_rzg,fy_rzg,fz_rzg
      common fx_rzg(0:MN-1),fy_rzg(0:MN-1),fz_rzg(0:MN-1)
      double precision lmet_rzg,mz_rzg,fh_rzg
      common lmet_rzg(0:MN-1),mz_rzg(0:MN-1),fh_rzg(0:MN-1)

c *** for cooling tables ***
      integer MNT_CRTB,MNMET_CRTB,MNNH_CRTB,MNZ_CRTB
      parameter (MNT_CRTB=100)
      parameter (MNMET_CRTB=10)
      parameter (MNNH_CRTB=10)
      parameter (MNZ_CRTB=100)

c ***** for metal cooling and heating ***
      integer nz_crtb,nnh_crtb,nmet_crtb,nt_crtb
      double precision z_crtb,nh_crtb,met_crtb,t_crtb
      double precision cr_crtb,hr_crtb,myu_crtb
      common /CRTBI/nz_crtb,nnh_crtb,nmet_crtb,nt_crtb
      common /CRTBD/z_crtb(0:MNT_CRTB-1),nh_crtb(0:MNNH_CRTB)
     & ,met_crtb(0:MNMET_CRTB),t_crtb(0:MNT_CRTB)
      common /CRTBD/cr_crtb(0:MNT_CRTB,0:MNMET_CRTB
     &  ,0:MNNH_CRTB,0:MNZ_CRTB)
     &  ,hr_crtb(0:MNT_CRTB,0:MNMET_CRTB
     &  ,0:MNNH_CRTB,0:MNZ_CRTB)
     &  ,myu_crtb(0:MNT_CRTB,0:MNMET_CRTB
     &  ,0:MNNH_CRTB,0:MNZ_CRTB)

c ***** for Work ****
c * Non-active particle list *  
      integer nalist,talist
      common nalist(0:MNS-1)
c * temporary active particle list *      
      common talist(0:MNS-1)
c 
      integer list,node
      common list(0:MNS-1),node(0:MNS-1)
c * build tree *
c * Particle No. of not-finished Particle *      
      integer pn_nfp
      common pn_nfp(0:MNS-1)
c * Node  No. included this particle *      
      integer nd_nfp
      common nd_nfp(0:MNS-1)
c * Which subcell is this particle in? *      
      integer c_nfp
      common c_nfp(0:MNS-1)
c * label for not-finished particle *      
      integer label_nfp
      common label_nfp(0:MNS-1)
c * number of particle in subcell *      
      integer np_sc
      common np_sc(0:MNNODE-1)
c * parent node *      
      integer pare_tr
      common pare_tr(0:MNNODE-1)
c * node No. of subcell *      
      integer nd_sc,pare_sc,c_sc
      common nd_sc(0:MNNODE-1)
      common pare_sc(0:MNNODE-1)
      common c_sc(0:MNNODE-1)
c *for parent cell *
c * Is parent have a child ? no 0 yes pre-node *      
      integer flag_pc
      common flag_pc(0:MNNODE-1)
c * for work *
      double precision tx,ty,tz
      common /WD/tx(0:MNS-1),ty(0:MNS-1),tz(0:MNS-1)
      double precision tcx,tcy,tcz
      common /WD/tcx(0:7),tcy(0:7),tcz(0:7)

      integer NCPU
      parameter (NCPU=64)

c ***** for MPI *****
      integer nprocs,myrank
      integer jsta,jend            
      integer istatus,istatusr,ierr
      integer idisp,jjlen
      integer ireqs,ireqr
      common nprocs,myrank
      common jsta,jend                        
c      common istatus(MPI_STATUS_SIZE),istatusr(MPI_STATUS_SIZE),ierr
      common idisp(0:NCPU-1),jjlen(0:NCPU-1)
      common ireqs(0:NCPU-1),ireqr(0:NCPU-1)
