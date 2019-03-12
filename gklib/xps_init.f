c
c ... XPS_GRAF.F ... subroutines for producing PostScript output
c
c ... Gerard J Kleywegt @ 920724
c
c ... Update 920724
c ... Update 920803,12,24
c ... Update 920907
c ... Update 930223
c ... Update 930305
c ... Update 930805
c ... Update 931220
c ... Update 940415
c ... Update 950215
c ... Update 950529,30
c ... Update 950705
c ... Update 960209,23
c ... Update 961031
c ... Update 961101
c ... Update 961209
c ... Update 19981117,24,26,29
c ... Update 19981207,08,09,16
c ... Update 20051205
c
c ================================================================
c
      subroutine xps_init ()
c
c ... initialise the subroutines
c
      include 'xps.incl'
c
code ...
c
      psunit = 99
      psfile = 'out.ps'
      psopen = .false.
      psinit = .true.
      pspola = .false.
      pshide = .false.
      psxlab = 'No_label'
      psylab = 'No_label'
c
      psxsca = 1.0
      psxoff = 0.0
      psysca = 1.0
      psyoff = 0.0
c
      psnuml = 0
      psdash = 0
      pscolr = 0
      psntxt = 0
c
      poxmin = 0.0
      poxmax = 1.0
c
c ... initialise SGI colours
c
      call xvrml_init ()
c
c      write (*,*) 'XPS routines initialised'
c
      return
      end
