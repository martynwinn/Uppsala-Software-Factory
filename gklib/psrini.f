c
c
c
      subroutine psrini (iunit,psfile,prognm)
c
c ... initialise PostScript Ramachandran plot with our
c     latest definition of core regions
c
      implicit none
c
      real phi,psi
c
      integer coregn(37,37)
      integer iunit,i,j
c
      character psfile*(*),prognm*(*)
      character labx*80,laby*80
c
code ...
c
      call defcor (coregn)
c
      call xps_init ()
      call xps_open (iunit,psfile,prognm)
      call xps_scale (-180.0,180.0,-180.0,180.0)
      call xps_stroke ()
      call xps_ps_comment ('Boxes for core regions')
c
      do i=1,36
        do j=1,36
          if (coregn(i,j) .eq. 1) then
            phi = float(i-1)*10.0 - 180.0
            psi = float(j-1)*10.0 - 180.0
            call xps_light_box (phi,phi+10.0,psi,psi+10.0)
          end if
        end do
      end do
c
      call xps_move (-180.,-180.)
      call xps_draw (+180.,-180.)
      call xps_draw ( 180., 180.)
      call xps_draw (-180., 180.)
      call xps_draw (-180.,-180.)
c
c --- phi axis ticks
c
      call xps_ps_comment ('Phi axis ticks')
      do 300 i=1,11
        phi = -180.+i*30.
        call xps_move (phi,-180.)
        call xps_draw (phi,-175.)
        call xps_move (phi, 180.)
        call xps_draw (phi, 175.)
300   continue
c
c --- psi axis ticks
c
      call xps_ps_comment ('Psi axis ticks')
      do 310 i=1,11
        psi = -180+i*30.
        call xps_move (-180.,psi)
        call xps_draw (-175.,psi)
        call xps_move ( 180.,psi)
        call xps_draw ( 175.,psi)
310   continue
c
      labx = 'PHI mapped to [-180,180>'
      laby = 'PSI mapped to [-180,180>'
      call xps_label (labx,laby)
      call xps_ps_comment ('Ramachandran initialisation done')
c
      return
      end
