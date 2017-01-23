c /************************************************
c      global.f for PBTREESPH ver.p12
c  21 July, 2000    produced by D.KAWATA
c ************************************************/

      include 'define.f'

c /*****   for Particle   *****/
c /* Mass */      
      double precision m_p,ms0_p,eps_p
c /* Virtual Position Xn+1,Xn */      
      double precision x_p,y_p,z_p
c /* Correct position */      
      double precision xc_p,yc_p,zc_p
c /* Velosity Vn-1/2,Vn+1/2 */      
      double precision vx_p,vy_p,vz_p
c /* Vn+1 */
      double precision vnx_p,vny_p,vnz_p
c /* Virtual Vn+1 */      
      double precision vvnx_p,vvny_p,vvnz_p
c /* Internal Energy Un+1,Un */       
      double precision u_p,pu_p
c /* Smoothing Length hn+1,hn */      
      double precision h_p
c /* Density rhon+1,rhon */       
      double precision rho_p
c /* Presure Pn+1,Pn */      
      double precision p_p
c /* Acceleration */       
      double precision ax_p,ay_p,az_p
c /* Sound Velocity */      
      double precision cs_p
c /* div V */       
      double precision div_v_p
c /* |rot V| */      
      double precision arot_v_p
c /* Max Viscosity Coefficient */      
      double precision maxmyu_p
c /* dv/dt */      
      double precision dvx_p,dvy_p,dvz_p
c /* du/dt_n,du/dt_n+1 */      
      double precision pdu_p,ndu_p
c /* Number of neigbor within 2xh */      
c      integer numnb_p
      double precision numnb_p
c /* For Star Formation */
c /* energy of supernova */      
      double precision Gsn_p
c /* For Heavy Elements unit Solar Mass */
      double precision mzHe_p,mzC_p,mzN_p,mzO_p,mzNe_p,mzMg_p,
     & mzSi_p,mzFe_p,mzZ_p
c /* For Total Ejected mass by SNe unit Solar Mass */
      double precision tnsn_p,tmej_p,tmzHe_p,tmzC_p,tmzN_p,tmzO_p,
     & tmzNe_p,tmzMg_p,tmzSi_p,tmzFe_p,tmzZ_p
c /*** For Individual time step ***/
c /* Individual time bin */      
      double precision dt_p
c /* system time bin */      
      double precision lt_p
c /* Virtual Time bin */      
      double precision vdt_p
c /* for cooling */
      double precision ram_p
c /* star formtion time */
      double precision ts_p
c /* flag for SNIa */
      integer flagc_p      

c /***** For Dark Matter *****/
c /* mass */      
      double precision m_dm,eps_dm
c /* Virtual Position */      
      double precision x_dm,y_dm,z_dm
c /* Correct Position */      
      double precision xc_dm,yc_dm,zc_dm
c /* Velocity V(n+1/2) * /     
      double precision vx_dm,vy_dm,vz_dm
c /* Virtual Velocity V(n+1/2) */      
      double precision vnx_dm,vny_dm,vnz_dm
c /* Acceration */      
      double precision dvx_dm,dvy_dm,dvz_dm
c /* for Individual time step */
c /* Individual time bin */      
      double precision dt_dm
c /* time in system time step */      
      double precision lt_dm
c /* virtual time bin */      
      double precision vdt_dm

c /***** for Tree ****/
c /* Number of nodes */      
      integer num_tr
c /* number of contained particle */      
       integer np_tr
c /* name of Particle */      
      integer pn_tr
c /* length of side */       
      double precision l_tr
c /* first child node */      
      integer daughter_tr
c /* next node */      
      integer next_tr
c /* parent node */
      integer pare_tr
c /* level info */
c      integer level_tr,llist_tr      
c /* Coordinate of center */      
      double precision cx_tr,cy_tr,cz_tr
c /* center of mass */
      double precision cmx_tr,cmy_tr,cmz_tr
c /* total of mass */      
      double precision mass_tr
c /* distance between cm and center */      
      double precision delta_tr
c /* for Multipole Expansion */
c      double precision mx_tr,my_tr,mz_tr
      double precision mxx_tr,myy_tr,mzz_tr
      double precision mxy_tr,myz_tr,mzx_tr

c /***** For Dark Matter Tree *****/
c /* number of contained particle */       
      integer np_dmtr
c /* name of Particle */      
      integer pn_dmtr
cc /* length of side */       
      double precision l_dmtr
c /* first child node */      
      integer daughter_dmtr
c /* next node */      
      integer next_dmtr
      double precision cx_dmtr,cy_dmtr,cz_dmtr
c /* Coordinate of center */ 
      double precision cmx_dmtr,cmy_dmtr,cmz_dmtr
c /* center of mass */
c /* total of mass */      
      double precision mass_dmtr
c /* distance between cm and center */      
      double precision delta_dmtr
c /* for Multipole Expansion */
c      double precision mx_dmtr,my_dmtr,mz_dmtr
      double precision mxx_dmtr,myy_dmtr,mzz_dmtr
      double precision mxy_dmtr,myz_dmtr,mzx_dmtr

c /***** For low resolutino Dark Matter Tree *****/
c /* number of contained particle */       
      integer np_ldmtr
c /* name of Particle */      
      integer pn_ldmtr
cc /* length of side */       
      double precision l_ldmtr
c /* first child node */      
      integer daughter_ldmtr
c /* next node */      
      integer next_ldmtr
      double precision cx_ldmtr,cy_ldmtr,cz_ldmtr
c /* Coordinate of center */ 
      double precision cmx_ldmtr,cmy_ldmtr,cmz_ldmtr
c /* center of mass */
c /* total of mass */      
      double precision mass_ldmtr
c /* distance between cm and center */      
      double precision delta_ldmtr
c /* for Multipole Expansion */
c      double precision mx_ldmtr,my_ldmtr,mz_ldmtr
      double precision mxx_ldmtr,myy_ldmtr,mzz_ldmtr
      double precision mxy_ldmtr,myz_ldmtr,mzx_ldmtr

c /***** For Individual Time Step *****/
c /* list of active particles for gas and star */      
      integer list_ap
c /* list of active particles for DM */      
      integer list_adm

c /***** Tree for SPH particles _gtr *****/
c /* number of contained particle */       
      integer np_gtr
c /* name of Particle */
      integer pn_gtr
c /* length of side */      
      double precision l_gtr,hm_gtr
      double precision cx_gtr,cy_gtr,cz_gtr
c /* first child node */      
      integer daughter_gtr
c /* next node */      
      integer next_gtr

c /***** for Work ****/
      integer list,node
c /* build tree */
c /* Particle No. of not-finished Particle */      
      integer pn_nfp
c /* Node  No. included this particle */      
      integer nd_nfp
c /* Which subcell is this particle in? */      
      integer c_nfp
c /* label for not-finished particle */
      integer label_nfp
c /* number of particle in subcell */      
      integer np_sc
c /* node No. of subcell */      
      integer nd_sc
      integer pare_sc
      integer c_sc
c /* for parent cell */
c /* Is parent have a child ? no 0 yes pre-node */
      integer flag_pc
c /* for work */
      double precision tx,ty,tz
      double precision tcx,tcy,tcz
c /* for time */
c /* time criterion 1,2 */      
      double precision dtcr1,dtcr2
c /* |dv/dt| */    /* |v| */
      double precision adv,av
c /* Individual time */      
      double precision dti
c /* Individual time for DM */      
      double precision dtdmi
c /* Non-active particle list */      
      integer nalist
c /* temporary active particle list */      
      integer talist
c /* for updateu */
      double precision du0,ram,hu,lu,fx,fhx,ulim
 
c /***** Definition of Structure *****/
c /* TimeStep DTn DTn-1*/      
      double precision TM_dt
c /* Hubble time */      
      double precision TM_t0
c /* Total Time */      
      double precision TM_tot
c /* System time bin */      
      double precision TM_sdt
c /* Minimum local time bin */      
      double precision TM_ldtmin
c /* local Time in system */      
      double precision TM_lt
      integer TM_k,TM_kn
      
c /***** For OFPara *****/
c /* This structure have Number of output file. */
c /* Make file for display by step */
      integer OF_count,OF_scount
      double precision OF_odt,OF_todt

c /***** For System Information *****/
c /*** 0: live halo, otherwise: fixed halo ***/
      integer SI_nh,SI_nd
c /* number of high resolution DM ***/
      integer SI_ndm1
c /* 0 = include cooling, 1 = without */      
      integer SI_flagc
c /* 0 = include starformation 1 = without */      
      integer SI_flags
c /* 0 = include hydrodynamics otherwise = no */      
      integer SI_flaghyd
c /* flag for making output file */
      integer SI_flagout
c /* Cosmological Value */
c      double precision SI_h0,SI_omg0,SI_lam0,SI_a
c/* table for scale factor and time */
c      double precision at_ttb,ta_ttb
c /* for memory of feedback history */
      double precision SI_esnIa,SI_esnII,SI_snzg,SI_snzs
      double precision SI_esng,SI_esns

c /***** For Yields lookup table *****/
      double precision z_ytb,t_ytb,mej_ytb,nsn_ytb,mzHe_ytb,
     &  mzC_ytb,mzN_ytb,mzO_ytb,mzNe_ytb,mzMg_ytb,mzSi_ytb,
     &  mzFe_ytb,mzZ_ytb
c /***** For Cooling Lookup Table *****/
      double precision z_ctb,cr_ctb

c /***** for MPI *****/
      integer nprocs,myrank
      integer jsta,jend            
      integer istatus      
      integer idisp,jjlen
      integer ireqs,ireqr
      integer tivs,tivr
      double precision tdvs,tdvr

c /**** for Random Numbrer ***/
      integer idum,ierr
