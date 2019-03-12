c
c
c
      subroutine psrinp (lunhp,psfile,prognm,how)
c
c ... set up PostScript Ramachandran file with our old definition
c
      implicit none
c
      integer lunhp,i
c
      real x1(2,42),phi,psi
c
      character labx*40,laby*40,psfile*(*),prognm*(*),how*1
c
      data x1 /
     $	-176.0,-67.1 ,-176.0,-20.,
     $	-147.9,-35.3,-129.5,-31.5,
     $	-122.1,-19.5,-129.4,-13.2,
     $	-176.,40.4,-176.,180.,
     $	-53.,180.,-39.8,151.6,
     $	-38.7,96.9,-68.8,77.1,
     $	-86.,66.2,-84.8,38.1,
     $	-57.0,-1.1,-35.,-18.3,
     $	-35.5,-64.8,-176.,-67.,
     +	34.7,87.4,58.2,110.4,
     $	58.2,11.2,35.0,27.5,
     $	34.7,87.4,
     +	-175.4,-180.,-175.4,-169.5,
     $	-61.6,-169.5,-57.0,-180.,
     $	-175.4,-180.,
     +	-161.,146.,-161.,98.,
     $	-59.,98.,-59.,177.,
     $	-138.,177.,-138.,160.,
     $	-161.,146.,
     +	-161.6,-52.7,-61.9,-51.,
     $	-61.9,-34.1,-111.2,-34.1,
     $	-127.5,-43.3,-161.6,-43.3,
     $	-161.6,-52.7/
c
code ...
c
      call xps_init ()
      call xps_open (lunhp,psfile,prognm)
      call xps_ps_comment ('Initialise old Ramachandran plot')
c
      if (how .eq. 'P') then
        call xps_polar (0.0,360.0)
        do i=1,42
          call fix360 (x1(1,i))
        end do
c
      else
        call xps_scale (-180.0,180.0,-180.0,180.0)
        do i=1,42
          call fixang (x1(1,i))
        end do
c
        call xps_move (-180.,-180.)
        call xps_draw (+180.,-180.)
        call xps_draw ( 180., 180.)
        call xps_draw (-180., 180.)
        call xps_draw (-180.,-180.)
c        call xps_move (   0.,-180.)
c        call xps_draw (   0., 180.)
c        call xps_move (-180.,   0.)
c        call xps_draw ( 180.,   0.)
c
c ---	phi axis ticks
c
        do 300 i=1,11
          phi = -180.+i*30.
          call xps_move (phi,-180.)
          call xps_draw (phi,-175.)
300     continue
c
c ---	psi axis ticks
c
        do 310 i=1,11
          psi = -180+i*30.
          call xps_move (-180.,psi)
          call xps_draw (-175.,psi)
310     continue
c
      end if
c
c ---	text
c
      if (how .eq. 'P') then
        labx = 'PHI mapped to [0,360>'
        laby = 'PSI mapped to [0,360>'
      else
        labx = 'PHI mapped to [-180,180>'
        laby = 'PSI mapped to [-180,180>'
      end if
      call xps_label (labx,laby)
c
c ---	boundary
c
      call xps_dash ()
      call xps_dash ()
c
      call xps_move (x1(1,1),x1(2,1))
      do 340 i=2,18
340     call xps_draw (x1(1,i),x1(2,i))
      call xps_move (x1(1,19),x1(2,19))
      do 350 i=20,23
350     call xps_draw (x1(1,i),x1(2,i))
      call xps_move (x1(1,24),x1(2,24))
      do 360 i=25,28
360     call xps_draw (x1(1,i),x1(2,i))
      call xps_move (x1(1,29),x1(2,29))
      do 370 i=30,35
370     call xps_draw (x1(1,i),x1(2,i))
      call xps_move (x1(1,36),x1(2,36))
      do 380 i=37,42
380     call xps_draw (x1(1,i),x1(2,i))
c
      call xps_solid ()
      call xps_ps_comment ('Ramachandran initialisation done')
c
      return
      end
