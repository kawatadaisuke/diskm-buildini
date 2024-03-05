

## pgen_dbulhalo

 F90 version of dbulhalo
 20/12/2016 Miyamoto-Nagai disk only

dbulhalo/
 generate a disk model see Kawata et al. (2014) for more detail.
http://adsabs.harvard.edu/abs/2014MNRAS.443.2757K


## externally created files

dbulhalo/ini/eos.dat
 from radproc/cloudy/dkmod/chrate/data/eos/eos.dat

ini/mgr.dat
 from /Users/dkawata/work/fsim/metal/WW+vdHG/tzycmdiff

## Requirements:

RAN1.FOR from Numerical Recipes in Fortran 77, which should be copied to ran.f. 

## Oputput

### gas.dat (unformatted)

xp,yp,zp,vxp,vyp,vzp,mpart,rhop,ug,0.0d0,0

xp, yp, zp ... positions (100 kpc)
vxp, vyp, vzp ... velocities (207.4 km/s)
mpart ... particle mass (1e12 Msun)
rhop ... density (1e12 Msun/(100 kpc)^3)
ug ... internal energy (P=(gam-1)*rho*ug)

### gas-metal.dat

mzHeg,mzCg,mzNg,mzOg,mzNeg,mzMgg,mzSig,mzFeg,mzZg

mass of 4He, 12C, 14N, 16O, 20Ne, 24Mg, 28Si, 56Fe, metal in the particle (Msun)

We assume Solar Abundance (meteorites from Woosley and Weaver 1995)
c /* 1H,4He,12C,14N,16O,20Ne,23Na,24Mg,27Al,28Si,32S,36Ar.40Ca,56Fe,58Ni */
data xsol/0.706d0,0.275d0,3.03d-3,1.11d-3,9.59d-3,1.62d-3,
  & 3.34d-5,5.15d-4,5.81d-5,6.53d-4,3.96d-4,7.74d-5,5.99d-5,
  & 1.17d-3,4.94d-5/
c /* solar abundance */
parameter (XSZ=0.019d0)

Using these you can compute different abundance ratios for each particle as follows.

mzHg = mpart*1.0e12 - mzZg ... Hydrogen mass

[Fe/H] = log10(mzFeg/mzHg)-log10(xsol(14)/xsol(1))
[O/Fe] = log10(mzOg/MzFeg)-log10(xsol(5)/xsol(14))

You the following constansts, you can also get temperature in K.
      parameter (GAM=5.0d0/3.0d0)
c /* tempature ( K ) */  
      parameter (TUK=1.0e4)
c /* k/m unit /( cm^2 s^-2 K^-1) */       
      parameter (K_MU=4.3e10)
c /* Boltzmann constarnt */    
      parameter (KCGS=1.381e-16)
c /* average molecular weight */      
      parameter (MYU=0.6d0)
c /* proton mass */      
      parameter (MP=1.67265e-24)      
c /* (KCGS/(MYU*MP))/K_MU)))*TUK */
      parameter (TPRHO=0.00320014d0)

temperature = ug*(GAM-1.0d0)*TUK*myup/(MYU*TPRHO) (K)



