c /****************************************
c   define.f for TREESPH for ellfor ver.p12
c 21 July, 2000    produced by D.KAWATA
c ****************************************/
      implicit none

      double precision M_PI
      parameter (M_PI=3.141592654d0)
c /*** Definition Max Number ***/
c /*** Number of Galaxies ***/
      double precision MNG
      parameter (MNG=2)
c /*** for gas particles ***/
c /* Max Number of Particle */      
c      integer MN      
c      parameter (MN=3100000)
c      integer MNS
c      parameter (MNS=3100000)
c/**** for neigbor list ***/
c /* Neigbor < NBMG x smoothing length */
      double precision NBMG      
      parameter (NBMG=2.0d0)
c /* Hernquist & Katz 89 (2.17) , Ns */
c /* Characteristic Number of Neighbor */      
      double precision  NNBS
      parameter (NNBS=40.0d0)
      double precision NNBSMIN
      parameter (NNBSMIN=10.0d0)
      
c /*****   For Tree   *****/
c /* Max Number of node in tree */
      integer MAXNODE
      parameter (MAXNODE=10)

c /*** for DM particles ***/
c /* Max Number of Total DM Particle */
c      integer MNDM
c      parameter (MNDM=6100000)
c /*** For Tree of high lesolution DM particles***/
c /* Max Number of node in tree */      
      integer MAXNODEDM
      parameter (MAXNODEDM=10)
c /*** For Tree of low resolution DM particles ***/
c /* Max Number of node in tree */      
      integer MAXNODELDM
      parameter (MAXNODELDM=1)

c /*** Definition Parameters ***/
c /*** for Tree ***/
      double precision MGROOT,THETA
c /* Margin for Size of root tree */      
      parameter (MGROOT=1.001d0)
c /* Torelance Parameter */
      parameter (THETA=0.8d0)

c/*** for at Table ***/
      integer NATTABLE
      parameter (NATTABLE=100000)

c /*****   For system value   *****/
      integer IMAX
c /* Itteration limit */
      parameter (IMAX=10000)
      double precision FH,TLIMIT,CLIMIT,EPSC,DU_UPU
c /* h > FH*EPS */      
      parameter (FH=1.0d0)
c /*  cooling limit (10^4K) */      
      parameter (TLIMIT=1.0d0)
      parameter (CLIMIT=1.0d0)
c /* cooling error */      
      parameter (EPSC=1.0d-4)
c /* dlog U in updateu */
      parameter (DU_UPU=0.1d0)

c /*** For star formation ***/
      double precision MFSF,DTH,CSEFF
      double precision FVSN,SNEU
c /* mg_p < m_p x MFSF -> star particle */      
      parameter (MFSF=0.05d0)
c /* in our unit 6.77e-26 g/cm^3 ref. NW93 */      
c /*** 7e-26 ***/
c      parameter (DTH=1.034d0)
c /*** 2e-25 about nh = 0.1 ***/
      parameter (DTH=0.2954d0)
c /* star formation efficiency coefficient */
      parameter (CSEFF=0.5d0)
c /* supernova energy by solar mass */      
c      parameter (SNU=1.0e51)
c /* supernova energy fraction for kinematic energy */
       parameter (FVSN=0.0d0)
c /* supernova energy (unit 8.56e59 erg) */
c       parameter (SNEU=0.1d0*1.168224e-9)
       parameter (SNEU=10.0d0*1.168224e-9)
       
c /*** Solar Abundances from WW95 ***/
c /*** H1 ***/
      double precision XHSOL
      parameter (XHSOL=0.706d0)
c /*** Fe56 ***/
      double precision XFESOL
      parameter (XFESOL=0.00117d0)
      
c /*** For Artificial Viscosity ***/
      double precision V_ETA,V_ALPHA,V_BETA,SHAREV
      parameter (V_ETA=0.1d0)
      parameter (V_ALPHA=1.0d0)
      parameter (V_BETA=2.0d0)
      parameter (SHAREV=1.0d0)

c /*****   for Time Step   *****/
      double precision MGTL,MGTU,CCFL,CDYN,CGRAV
c /* Mergin Time Interval lower  */      
      parameter (MGTL=0.9d0)
c /*                      higher */      
      parameter (MGTU=1.1d0)
c /* Courant Number */      
      parameter (CCFL=0.3d0)
      parameter (CDYN=0.1d0)
      parameter (CGRAV=0.4d0)
c /*** For Time Interval ***/
      double precision DTGMAX
c /* Max of Time Interval Gradient */      
      parameter (DTGMAX=1.05d0)

c /*** For Arithmetic Value ***/
      double precision THIRD,INF,MININF,V34
      parameter (THIRD=1.0d0/3.0d0)
      parameter (V34=3.0d0/4.0d0)
      parameter (INF=1000000.0d0)
      parameter (MININF=1.0d-6)
      integer MAXFILENAME,MAXSTEP
      parameter (MAXFILENAME=25)
      parameter (MAXSTEP=1000000)

c /*** Normalization Unit ***/
      double precision MUSM,LUPC,LUKPC,LUCM,DU,PU,TUK,TMUYR,
     &  TMUGYR,TMUS,VUKMS,ERU,K_MU
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
c /* Density */      
      parameter (DU=6.77e-26)
c /* Pressure */
      parameter (PU=2.91e-11)
c /* tempature ( K ) */  
      parameter (TUK=1.0e4)
c /* time ( yr ) */      
      parameter (TMUYR=4.71e8)
c /* time ( Gyr ) */      
      parameter (TMUGYR=0.471d0)
c /* time (s) */      
      parameter (TMUS=1.488e16)
c /* km/s */      
      parameter (VUKMS=207.4d0)
c /* cooling (energy) rate (erg s^-1 cm^-3) */      
      parameter (ERU=1.96e-27)
c /* k/m unit /( cm^2 s^-2 K^-1) */       
      parameter (K_MU=4.3e10)
c /* erg */      
c      parameter (EU=8.56e59)

c /*** For Physical Value ***/
      double precision G,KCGS,MYU,MP,GCGS,TEQ,ZEQ1,TPRHO,H0_1,HUB0
c /* Gravitatrional Constant */      
      parameter (G=1.0d0)
c /* Boltzmann constarnt */    
      parameter (KCGS=1.381e-16)
c /* average molecular weight */      
      parameter (MYU=0.6d0)
c /* proton mass */      
      parameter (MP=1.67265e-24)
c /* G(cgs) */      
      parameter (GCGS=6.672e-8)
c /* Tempature (at t=teq) (OMG h^2 K) */      
      parameter (TEQ=10.7184e4)
c /* Redshift (at t=teq) (OMG h^2) */    
      parameter (ZEQ1=3.9e4)
c /* (KCGS/(MYU*MP))/K_MU)))*TUK */
      parameter (TPRHO=0.00320014d0)
c /* Hubble constant 9.78x10^9xh^-1 yr h = 1.0*/      
      parameter (H0_1=(1.0d0/((100.0d0*1.0d5)/(3.086d24)))/
     & ((((3600.0d0*24.0d0)*365.0d0)*4.0d0+(3600.0d0*24.0d0))/4.0d0))
      parameter (HUB0=100.0d0/(VUKMS*10.0d0))

c /*** for lookup table of yields ***/
c /*** Number of Metallicities and Time ***/
      integer NYTZ,NYTT
      parameter (NYTZ=60,NYTT=200)
c /*** for lookup table of cooling ***/
c /*** Number of Metallicities of Temperature ***/
      integer NCTZ,NCTT
      parameter (NCTZ=8,NCTT=140)
      double precision DT_CT,LT_CT,ZEROZ_CT
      parameter (DT_CT=0.05d0,LT_CT=2.0d0,ZEROZ_CT=-20.0d0)

      integer MNMGR,MNZMGR
      parameter (MNMGR=100,MNZMGR=200)
c /***** for orbit table *****/
      integer MNORBTB
      parameter (MNORBTB=10000)

c/*** For Parallel ***/
      integer NCPU
      parameter (NCPU=41)
