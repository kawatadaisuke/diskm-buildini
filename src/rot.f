c /***************************************
c      rot.f  ver.9
c  21 Nov. 2011    written by D.KAWATA
c ***************************************/ 
c70*******************************************************************

      subroutine rot(xp,yp,zp,vxp,vyp,vzp,ig,igal,omgal,flagrot)
      include 'define.f'
      integer MNGAL
      parameter (MNGAL=2)
      integer ig,flagrot
      double precision xp,yp,zp,vxp,vyp,vzp
      double precision igal(MNGAL),omgal(MNGAL)
      double precision tx,ty,tz

      if(flagrot.eq.0) then
c /*** rotation using y-axis ***/
c *** position ***
        tx=xp
        ty=yp
        tz=zp
        xp=dsin(igal(ig))*tz+dcos(igal(ig))*tx
        yp=ty
        zp=dcos(igal(ig))*tz-dsin(igal(ig))*tx
c *** velocity ***
        tx=vxp
        ty=vyp
        tz=vzp
        vxp=dsin(igal(ig))*tz+dcos(igal(ig))*tx
        vyp=ty
        vzp=dcos(igal(ig))*tz-dsin(igal(ig))*tx
c /*** rotation using z-axis ***/
c *** position ***
        tx=xp
        ty=yp
        tz=zp
        xp=dcos(omgal(ig))*tx-dsin(omgal(ig))*ty
        yp=dsin(omgal(ig))*tx+dcos(omgal(ig))*ty
        zp=tz
c *** velocity ***
        tx=vxp
        ty=vyp
        tz=vzp
        vxp=dcos(omgal(ig))*tx-dsin(omgal(ig))*ty
        vyp=dsin(omgal(ig))*tx+dcos(omgal(ig))*ty
        vzp=tz
      else if(flagrot.lt.0) then
c /*** rotation using y-axis ***/
c *** position ***
        xp=-xp
        zp=-zp
c *** velocity ***
        vxp=-vxp
        vzp=-vzp
      endif

      return
      end
