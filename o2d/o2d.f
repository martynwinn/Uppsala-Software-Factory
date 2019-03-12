      program o2d
c
c ... O2D === 2D graphics for O
c
c ... Gerard Kleywegt @ 920511
c
c ... Modified 920512,13,14,15
c ... Modified 920717,24,29
c ... Modified 920810,18,24
c ... Modified 920901,22
c ... Modified 921110
c ... Modified 930202,03,23
c ... Modified 930414,19
c ... Modified 930607
c ... Modified 931028
c ... Modified 931129
c ... Modified 931220
c ... Modified 940214,15,941013,941230
c ... Modified 951013,951020
c ... Modified 981129,981130,981216
c ... Modified 990206,990824,26,990922
c ... Modified 000628
c ... Modified 010327,010725,011011
c
c ... TO DO LIST:
c
c     DONE * implement 1D & 2D integration in cursor mode
c     DONE * add 1D integration
c     DONE * add 2D integration
c     DONE * add PIE charts
c     DONE * simple use of colour in PostScript files
c     DONE * add REMARK lines as legend to figures
c     DONE * fixed offset bug in 2D contour plotting
c     DONE * support * format and MORE in CricketGraph option
c     DONE * GRID in 1D plot files
c     DONE * "naive" option to quickly plot O datablock
c     DONE * allow for PS conversion without graphics
c     DONE * do conversion to CricketGraph file
c     DONE * auto-update command names and parameters
c     DONE * allow for multiple plots in one window/PS file
c     DONE * create PostScript file from 2D contour data
c     DONE * popup option to list points above/below a certain Y-value
c     DONE * picking in scatter and 1D-line plots
c     DONE * read 'END   ' records from data files
c     DONE * store colours per window
c     DONE * clear window
c     DONE * reset original viewport in popup (store plot limits)
c     DONE * get data from O data blocks
c     DONE * highlighting
c     DONE * draw hori line/vert line/both at cursor (in popup)
c     DONE * default label is data point number
c     DONE * display one or all labels
c     DONE * PostScript files
c     DONE * axis labels in PostScript files
c     DONE * add a few PostScript comments here and there
c     DONE * RAMP_2D to produce colour-ramped filled squares
c     DONE * LINFIT keyword for 1D line and scatter plot
c     DONE * XYSTAT keyword for 1D line and scatter plot
c     DONE * implement box plot command and BOXPAR 1D keyword
c     DONE * implement spike plot command
c     DONE * implement CD plot command
c     DONE * implement command to count in 1D file and produce 2D file
c
      character prgnam*10,version*20
      parameter (prgnam = 'O2D', version = '070118/5.6.2')
c
#ifdef ALPHA
      integer*8 iaptr,ibptr,icptr,idptr
      integer*8 fmalloc
#else
      integer iaptr,ibptr,icptr,idptr
      integer fmalloc
#endif
c
      integer maxbuf,nb
code ...
c
      call gkinit (prgnam,version)
c
      maxbuf = 500000
      call extint ('GKBUFFER',maxbuf)
      maxbuf = max ( maxbuf , 10000 )
      call jvalut (' Allocate 2D plot buffers of size :',1,maxbuf)
c
      nb = 4 * maxbuf
      iaptr = fmalloc (nb)
      ibptr = fmalloc (nb)
      icptr = fmalloc (nb)
c
      nb = 2 * maxbuf
      idptr = fmalloc (nb)
c
#ifdef LINUX
      if (iaptr .eq. 0 .or. ibptr .eq. 0 .or. icptr .eq. 0 .or.
     +    idptr .eq. 0) then
#else
      if (iaptr .le. 0 .or. ibptr .le. 0 .or. icptr .le. 0 .or.
     +    idptr .le. 0) then
#endif
        call errcon ('Could not allocate sufficient memory')
        goto 99
      end if
c
      call doo2d (maxbuf,%val(iaptr),%val(ibptr),%val(icptr),
     +            %val(idptr))
c
      call ffree (iaptr)
      call ffree (ibptr)
      call ffree (icptr)
      call ffree (idptr)
c
   99 continue
      call gkquit()
c
      end
c
c
c
      subroutine doo2d (maxbuf,xdata,ibuff1,ibuff2,bits)
c
      include 'o2d.incl'
c
      integer maxbuf
      integer ibuff1(maxbuf),ibuff2(maxbuf)
      real xdata(maxbuf)
      character bits(2*maxbuf)*1
c
c ... PARAMETERS
c
      integer f1,f2
      parameter (f1=11, f2=12)
c
      integer maxopt
      parameter (maxopt = 10)
c
c ... VARIABLES
c
      real ddx,ddy,scx,scy,xlo,xhi,ylo,yhi,total,user,sys
c
      integer length,iindex,mcurs
c
      integer ierr,nopt,ndim,inext,icol,irow,iwin,idum
      integer iredrwmode,i,j,ic,ir,iw,npnt,nobj,ipnt,leng1
c
      logical once,dops,xtx_implemented,dovrml,lsolid,doramp
c
      character optpar(maxopt)*80
      character file1*128,file2*128,file3*128,file4*128
      character line*128,key*6,option*1024
      character winame*80,obname*40,gktype*10
c      character bits(2*maxbuf)*1
c
c ... DATA
c
      data once /.false./, ddx,ddy /10.0, 10.0/
c
code ...
c
      call gkdcpu (total,user,sys)
c
      linter = (xinter() .and. xtx_implemented())
c
      call xvrml_init ()
      dovrml = .false.
      doramp = .false.
      lsolid = .false.
c
      optpar (1) = 'open_window'
      if (.not. linter)  optpar (1) = 'cricketgraph'
      nopt = 1
      file1 = 'in.pl1'
      file2 = 'in.pl2'
      file3 = 'out.ps'
      file4 = 'o2d.wrl'
      ndim = 1
      ipnt = 1
      winame = 'No_name'
c
      xlo = 0.0
      xhi = 100.0
      ylo = 0.0
      yhi = 100.0
      last2w = -9
c
      if (linter) call grfini ()
      call xps_init ()
c
   10 continue
      if (linter) then
        call grafix_dims ()
        write (*,'(/99(a/:))')
     +  ' Valid options are:',
     +  ' ?',
     +  ' !',
     +  ' open_window    Dim Name Bkgr_colr Draw_colr',
     +  ' select_window  [Nr]',
     +  ' clear_window',
     +  ' cursor_mode    Dx Dy',
     +  ' naive_1D_plot  O_datablock_file [PostScript_file]',
     +  ' 1D_line_plot   Plot_file [PostScript_file]',
     +  ' scatter_plot   Plot_file [PostScript_file]',
     +  ' spike_plot     Plot_file [PostScript_file]',
     +  ' box_plot       Plot_file [PostScript_file]',
     +  ' histogram      Plot_file [PostScript_file]',
     +  ' pie_chart      Plot_file [PostScript_file]',
     +  ' 2D_contour     Plot_file [PostScript_file]',
     +  ' ramp_2D        Plot_file PostScript_file',
     +  ' CD_plot        Plot_file PostScript_file',
     +  ' topology       Plot_file [PostScript_file]',
     +  ' focus          Data_nr',
     +  ' integrate      Xlo Xhi [Ylo Yhi]',
     +  ' close_window',
     +  ' cricketgraph   Plot_file      CricketGraph_file',
     +  ' count_1d_to_2d 1D_Plot_file   2D_Plot_file',
     +  ' vrml_2D        Plot_file   VRML_file   line_solid',
     +  ' quit'
      else
        write (*,'(/99(a/:))')
     +  ' Valid options are:',
     +  ' ?',
     +  ' !',
     +  ' naive_1D_plot   O_datablock_file PostScript_file',
     +  ' 1D_line_plot    Plot_file        PostScript_file',
     +  ' scatter_plot    Plot_file        PostScript_file',
     +  ' spike_plot      Plot_file        PostScript_file',
     +  ' box_plot        Plot_file        PostScript_file',
     +  ' histogram       Plot_file        PostScript_file',
     +  ' pie_chart       Plot_file        PostScript_file',
     +  ' 2D_contour      Plot_file        PostScript_file',
     +  ' ramp_2D         Plot_file        PostScript_file',
     +  ' CD_plot         Plot_file        PostScript_file',
     +  ' cricketgraph    Plot_file        CricketGraph_file',
     +  ' count_1d_to_2d  1D_Plot_file     2D_Plot_file',
     +  ' vrml_2D         Plot_file   VRML_file   line_solid',
     +  ' quit'
      end if
c
c ... main event loop
c
   11 continue
c
      if (linter .and. numwin .gt. 0) call xtx_update_windows ()
c
      call gkdcpu (total,user,sys)
      if (total .ge. 1.0)
     +  write (*,'(a,3f10.1)') ' CPU total/user/sys :',total,user,sys
c
      write (option,'(20(a,1x))')
     +  (optpar(i)(1:leng1(optpar(i))),i=1,nopt)
      call pretty (option)
c
      call textin (' Option ?',option)
c
      if (length(option) .lt. 1) goto 11
      if (option(1:1) .eq. '!') goto 11
      if (option(1:1) .eq. '?') goto 10
c
      call extrop (option,nopt,maxopt,optpar,ierr)
      if (nopt .lt. 1 .or. ierr .ne. no_error) then
        call textut (' ERROR - No valid option :',option)
        goto 11
      end if
      call upcase (optpar(1))
c
cc      call asciut (' PARSED :',nopt,optpar)
c
      if (optpar(1)(1:1) .eq. 'Q') then
c
c ... QUIT
c
        goto 9999
c
      else if (optpar(1)(1:1) .eq. 'O') then
c
        if (.not. linter) then
          call errcon ('Non-interactive run; no on-screen graphics')
          goto 11
        end if
c
        optpar (1) = 'open_window'
c
c ... OPEN_WINDOW
c
        inext = iindex (not_open,nought,maxwin,wistat)
        if (inext .le. nought) then
          call errcon ('No windows available')
          call ivalut (' Maximum :',1,maxwin)
          goto 11
        end if
c
        if (nopt .gt. 1) then
          read (optpar(2),*,err=1002) ndim
          goto 1003
        end if
 1002   call ivalin (' Dimension [1,2] ?',1,ndim)
 1003   if (ndim .lt. 1 .or. ndim .gt. 2) then
          call ivalut (' Invalid dimension :',1,ndim)
          goto 10
        end if
        write (optpar(2),*) ndim
c
        if (nopt .gt. 2) then
          winame = optpar(3)
        else
          call textin (' Name ?',winame)
          optpar (3) = winame
        end if
c
        if (nopt .gt. 3) then
          read (optpar(4),*,err=1005) wicol1(inext)
          goto 1006
        end if
 1005   call ivalin (' Background colour ?',1,wicol1(inext))
        write (optpar(4),*) wicol1(inext)
 1006   continue
c
        if (nopt .gt. 4) then
          read (optpar(5),*,err=1008) wicol2(inext)
          goto 1009
        end if
 1008   call ivalin (' Drawing colour ?',1,wicol2(inext))
        write (optpar(5),*) wicol2(inext)
 1009   continue
c
        nopt = 5
c
        write (key,'(i3,a1)') inext,':'
        call remspa(key)
        winame = key(1:leng1(key))//winame
c
cc        print *,'window being opened'
c
        call xtx_open_window (winame,inext,ndim,ierr)
        if (ierr .ne. no_error) goto 11
c
cc        print *,'window opened'
c
        numwin = numwin + 1
        wistat (inext) = open
        windim (inext) = ndim
        nowwin = inext
        nowcol = 1
        nowrow = 1
c
c        print *,'select window'
c
        call xtx_select_window (nowwin,nowcol,nowrow)
c
c        print *,'erase window'
c
        call xtx_erase_window (wicol1(nowwin))
c
c        print *,'erase objects'
c
        call xtx_objects_erase (nowwin)
c
c        print *,'set viewports'
c
        call xtx_set_viewports (nowwin,1,1,wicol2(nowwin))
c
c        print *,'labels on'
c
        call xtx_labels_on_off (1)
        wincol (nowwin) = 1
        winrow (nowwin) = 1
        wintyp (nowwin) = no_plot
        winpnt (nowwin) = nought
        wobnum (nowwin) = nought
c
c        print *,'redraw window'
c
        call xtx_redraw_force ()
c
        call ivalut (' Opened window     :',1,inext)
        call textut (' Name              :',winame)
        call ivalut (' Graphics dim      :',1,ndim)
        call ivalut (' Background colour :',1,wicol1(nowwin))
        call ivalut (' Drawing colour    :',1,wicol2(nowwin))
c
      else if (optpar(1)(1:1) .eq. 'I') then
c
        if (.not. linter) then
          call errcon ('Non-interactive run; no integration')
          goto 11
        end if
c
        optpar (1) = 'integrate'
c
c ... INTEGRATE
c
        if (numwin .le. nought) then
          call errcon ('No windows open')
          goto 11
        end if
c
        if (wistat(nowwin) .ne. open) then
          call errcon ('Window not open')
          goto 11
        end if
c
        if (windim(nowwin) .eq. 2) then
          if (nowwin .ne. last2w) then
            call errcon (' Only data for *last* 2D plot is stored')
            goto 11
          end if
        end if
c
        if (nopt .gt. 1) then
          call str2r (optpar(2),xlo,ierr)
          if (ierr .ne. 0) goto 11
        else
          call fvalin (' Lower X ?',1,xlo)
        end if
c
        if (nopt .gt. 2) then
          call str2r (optpar(3),xhi,ierr)
          if (ierr .ne. 0) goto 11
        else
          call fvalin (' Upper X ?',1,xhi)
        end if
c
        call rlohi (xlo,xhi)
c
c ... 1D integration
c
        if (windim(nowwin) .eq. 1) then
          call gint1d (xlo,xhi)
          goto 11
        end if
c
c ... 2D integration
c
        if (nopt .gt. 3) then
          call str2r (optpar(4),ylo,ierr)
          if (ierr .ne. 0) goto 11
        else
          call fvalin (' Lower Y ?',1,ylo)
        end if
c
        if (nopt .gt. 4) then
          call str2r (optpar(5),yhi,ierr)
          if (ierr .ne. 0) goto 11
        else
          call fvalin (' Upper Y ?',1,yhi)
        end if
c
        call rlohi (ylo,yhi)
c
        call gint2d (xlo,xhi,ylo,yhi,maxbuf,xdata)
c
        goto 11
c
      else if (optpar(1)(1:3) .eq. 'CLO') then
c
        if (.not. linter) then
          call errcon ('Non-interactive run; no on-screen graphics')
          goto 11
        end if
c
        optpar (1) = 'close_window'
c
c ... CLOSE_WINDOW
c
        if (numwin .le. nought) then
          call errcon ('No windows open')
          goto 11
        end if
c
        if (wistat(nowwin) .eq. open) then
          call xtx_close_window ()
          wistat (nowwin) = not_open
          call ivalut (' Closed window :',1,nowwin)
          numwin = numwin - 1
          if (numwin .gt. nought) then
            nowwin = iindex (open,nought,maxwin,wistat)
            nowcol = 1
            nowrow = 1
            call xtx_select_window (nowwin,nowcol,nowrow)
c            call ivalut (' Selected window :',1,nowwin)
            call xtx_redraw_force ()
          else
            nowwin = nought
            nowcol = 1
            nowrow = 1
          end if
        else
          call errcon ('Window not open')
        end if
c
      else if (optpar(1)(1:3) .eq. 'CLE') then
c
        if (.not. linter) then
          call errcon ('Non-interactive run; no on-screen graphics')
          goto 11
        end if
c
        optpar (1) = 'clear_window'
c
c ... CLEAR_WINDOW
c
        if (numwin .le. nought) then
          call errcon ('No windows open')
          goto 11
        end if
c
        call xtx_select_window (nowwin,nowcol,nowrow)
        call xtx_erase_window (wicol1(nowwin))
        call xtx_objects_erase (nowwin)
        call xtx_set_viewports (nowwin,1,1,wicol2(nowwin))
        call xtx_labels_on_off (1)
        wincol (nowwin) = 1
        winrow (nowwin) = 1
        wintyp (nowwin) = no_plot
        winpnt (nowwin) = nought
        call xtx_redraw_force ()
c
      else if (optpar(1)(1:2) .eq. 'SE') then
c
        if (.not. linter) then
          call errcon ('Non-interactive run; no on-screen graphics')
          goto 11
        end if
c
        optpar (1) = 'select_window'
c
c ... SELECT_WINDOW
c
        if (numwin .le. nought) then
          call errcon ('No windows open')
          goto 11
        end if
c
        if (nopt .lt. 2) then
          call ivalut (' Selected window :',1,nowwin)
cc          call ivalut (' Selected column :',1,nowcol)
cc          call ivalut (' Selected row    :',1,nowrow)
          goto 11
        end if
c
        if (nopt .gt. 1) then
          read (optpar(2),*,err=3002) iwin
          goto 3003
        end if
 3002   call ivalin (' Window ?',1,iwin)
        write (optpar(2),*) iwin
        nopt = 2
 3003   if (iwin .le. nought .or. iwin .gt. maxwin) goto 11
        if (wistat(iwin) .eq. not_open) goto 11
c
        nowwin = iwin
        call ivalut (' Selected window :',1,nowwin)
c
        nowcol = 1
        nowrow = 1
c
cc        if (newpar .ge. 2) then
cc          call tcigpa (2,icol,ierr)
cc          if (ierr .ne. no_error) return
cc          if (icol .ge. 0 .and. icol .le. wincol(nowwin))
cc     +      nowcol = icol
cc        end if
c
cc        if (newpar .ge. 3) then
cc          call tcigpa (3,irow,ierr)
cc          if (ierr .ne. no_error) return
cc          if (icol .ge. 0 .and. icol .le. wincol(nowwin))
cc     +      nowrow = irow
cc        end if
c
        call xtx_select_window (iwin,nowcol,nowrow)
        call xtx_redraw_force ()
c
cc        call ivalut (' Selected column :',1,nowcol)
cc        call ivalut (' Selected row    :',1,nowrow)
c
c
      else if (optpar(1)(1:2) .eq. 'CU') then
c
        if (.not. linter) then
          call errcon ('Non-interactive run; no on-screen graphics')
          goto 11
        end if
c
        optpar (1) = 'cursor_mode'
c
c ... CURSOR_MODE
c
        if (numwin .le. nought) then
          call errcon ('No windows open')
          goto 11
        end if
c
        if (o2dbug .gt. no_debug) write (o2dout,*) ' win dim ',
     +    windim(nowwin)
        if (windim(nowwin) .gt. 2) then
          call errcon ('Not implemented for 3D graphics windows')
          goto 11
        end if
c
        if (nopt .gt. 1) then
          read (optpar(2),*,err=4002) ddx
          goto 4003
        end if
 4002   call fvalin (' Arrow key delta-X ?',1,ddx)
        write (optpar(2),*) ddx
 4003   ddx = max (0.0001, abs(ddx))
c
        if (nopt .gt. 2) then
          read (optpar(3),*,err=4005) ddy
          goto 4006
        end if
 4005   call fvalin (' Arrow key delta-Y ?',1,ddy)
        write (optpar(3),*) ddy
 4006   ddy = max (0.0001, abs(ddy))
c
        call gvalut (' Horizontal arrow step :',1,ddx)
        call gvalut (' Vertical arrow step   :',1,ddy)
c
        nopt = 3
c
c ... define appropriate popup menus
c
        idum = mcurs (define_popups)
c
        call xtx_get_redraw(iredrwmode)
        call xtx_set_redraw(1)
c
        call xtx_pan_zoom (ddx,ddy)
c
        call xtx_set_redraw(iredrwmode)
c
c ... delete appropriate popup menus
c
        idum = mcurs (delete_popups)
c
c ... tell which window/vpt is now active
c
cc        call ivalut (' Selected window :',1,nowwin)
cc        call ivalut (' Selected column :',1,nowcol)
cc        call ivalut (' Selected row    :',1,nowrow)
c
      else if (optpar(1)(1:2) .eq. 'CR') then
c
        optpar (1) = 'cricketgraph'
c
c ... CRICKETGRAPH
c
        if (nopt .gt. 1) then
          file1 = optpar(2)
        else
          call textin (' Plot file ?',file1)
          optpar (2) = file1
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
        else
          call textin (' CricketGraph file ?',file3)
          optpar (3) = file3
        end if
c
        nopt = 3
c
        gktype = 'CRIC'
        call plt2cg (f1,file1,file3,gktype,maxbuf,xdata,ierr)
c
      else if (optpar(1)(1:2) .eq. 'CO') then
c
        optpar (1) = 'count_1d_to_2d'
c
c ... COUNT_1D_TO_2D
c
        if (nopt .gt. 1) then
          file1 = optpar(2)
        else
          call textin (' 1D Plot file ?',file1)
          optpar (2) = file1
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
        else
          call textin (' 2D Plot file ?',file3)
          optpar (3) = file3
        end if
c
        nopt = 3
c
        gktype = 'COUN'
        call plt2cg (f1,file1,file3,gktype,maxbuf,xdata,ierr)
c
      else if (optpar(1)(1:1) .eq. 'V') then
c
        optpar (1) = 'vrml_2d'
c
c ... VRML_2D
c
        if (linter) then
c
          if (numwin .le. nought) then
            call errcon ('No graphics window')
            goto 11
          end if
c
          if (windim(nowwin) .ne. 2) then
            call errcon ('Not a 2D graphics window')
            goto 11
          end if
c
        end if
c
        if (nopt .gt. 1) then
          file2 = optpar(2)
        else
          call textin (' Plot file ?',file2)
          optpar (2) = file2
        end if
c
        if (nopt .gt. 2) then
          file4 = optpar(3)
        else
          call textin (' VRML file ?',file4)
          optpar (3) = file4
        end if
c
        if (nopt .lt. 4) then
          optpar (4) = 'L'
          call textin (' Lines or Solid faces (L/S) ?',optpar(4))
        end if
        call upcase (optpar(4))
        lsolid = (optpar(4)(1:1) .eq. 'S')
c
        nopt = 4
c
        dovrml = .true.
        dops = .false.
        doramp = .false.
c
        call plot2d (f2,file2,npnt,doramp,
     +               dovrml,lsolid,dops,file4,obname,nobj,ierr,
     +               maxbuf,ibuff1,ibuff2,xdata,bits)
c
        if (linter) then
c
          call xtx_redraw_force ()
c
c ... do window book-keeping
c
          wintyp (nowwin) = plot_2d
          winpnt (nowwin) = npnt
          wobnum (nowwin) = nobj
c
        end if
c
      else if (optpar(1)(1:1) .eq. 'N') then
c
        optpar (1) = 'naive'
c
c ... NAIVE
c
        if (linter) then
          if (numwin .le. nought) then
            call errcon ('No graphics window')
            goto 11
          end if
c
          if (windim(nowwin) .ne. 1) then
            call errcon ('Not a 1D graphics window')
            goto 11
          end if
        end if
c
        if (nopt .gt. 1) then
          file1 = optpar(2)
        else
          call textin (' Filename ?',file1)
          optpar (2) = file1
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
          dops = .true.
        else
          dops = .false.
        end if
c
        nopt = max (2, nopt)
c
        obname = '1D_plot'
c
        if (.not. linter .and. .not. dops) then
          call errcon ('Nothing to do')
          goto 11
        end if
c
        call naive ('LINE',
     +    f1,file1,file3,dops,npnt,obname,nobj,ierr)
c
        if (linter) then
c
          call xtx_redraw_force ()
c
c ... do window book-keeping
c
          wintyp (nowwin) = plot_1d
          winpnt (nowwin) = npnt
          wobnum (nowwin) = nobj
c
        end if
c
      else if (optpar(1)(1:1) .eq. '1') then
c
        optpar (1) = '1d_plot'
c
c ... 1D_PLOT
c
        if (linter) then
          if (numwin .le. nought) then
            call errcon ('No graphics window')
            goto 11
          end if
c
          if (windim(nowwin) .ne. 1) then
            call errcon ('Not a 1D graphics window')
            goto 11
          end if
        end if
c
        if (nopt .gt. 1) then
          file1 = optpar(2)
        else
          call textin (' Filename ?',file1)
          optpar (2) = file1
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
          dops = .true.
        else
          dops = .false.
        end if
c
        nopt = max (2, nopt)
c
        obname = '1D_plot'
c
        if (.not. linter .and. .not. dops) then
          call errcon ('Nothing to do')
          goto 11
        end if
c
        gktype = 'LINE'
        call plot1d (gktype,
     +    f1,file1,file3,dops,npnt,obname,nobj,ierr)
c
        if (linter) then
c
          call xtx_redraw_force ()
c
c ... do window book-keeping
c
          wintyp (nowwin) = plot_1d
          winpnt (nowwin) = npnt
          wobnum (nowwin) = nobj
c
        end if
c
      else if (optpar(1)(1:2) .eq. 'CD') then
c
        optpar (1) = 'CD_plot'
c
c ... CD_PLOT
c
        if (linter) then
          if (numwin .le. nought) then
            call errcon ('No graphics window')
            goto 11
          end if
c
          if (windim(nowwin) .ne. 1) then
            call errcon ('Not a 1D graphics window')
            goto 11
          end if
        end if
c
        if (nopt .gt. 1) then
          file1 = optpar(2)
        else
          call textin (' Filename ?',file1)
          optpar (2) = file1
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
          dops = .true.
        else
          dops = .false.
        end if
c
        nopt = max (2, nopt)
c
        obname = 'CD_plot'
c
        if (.not. linter .and. .not. dops) then
          call errcon ('Nothing to do')
          goto 11
        end if
c
        gktype = 'CDPLOT'
        call plot1d (gktype,
     +    f1,file1,file3,dops,npnt,obname,nobj,ierr)
c
        if (linter) then
c
          call xtx_redraw_force ()
c
c ... do window book-keeping
c
          wintyp (nowwin) = plot_1d
          winpnt (nowwin) = npnt
          wobnum (nowwin) = nobj
c
        end if
c
      else if (optpar(1)(1:1) .eq. '2') then
c
        optpar (1) = '2d_plot'
c
c ... 2D_PLOT
c
        if (linter) then
c
          if (numwin .le. nought) then
            call errcon ('No graphics window')
            goto 11
          end if
c
          if (windim(nowwin) .ne. 2) then
            call errcon ('Not a 2D graphics window')
            goto 11
          end if
c
        end if
c
        if (nopt .gt. 1) then
          file2 = optpar(2)
        else
          call textin (' Filename ?',file2)
          optpar (2) = file2
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
          dops = .true.
        else
          dops = .false.
        end if
c
        nopt = max (nopt, 2)
c
        obname = '2D_plot'
c
        if (.not. linter .and. .not. dops) then
          call errcon ('Nothing to do')
          goto 11
        end if
c
        dovrml = .false.
        doramp = .false.
        call plot2d (f2,file2,npnt,doramp,
     +               dovrml,lsolid,dops,file3,obname,nobj,ierr,
     +               maxbuf,ibuff1,ibuff2,xdata,bits)
c
        if (linter) then
c
          call xtx_redraw_force ()
c
c ... do window book-keeping
c
          wintyp (nowwin) = plot_2d
          winpnt (nowwin) = npnt
          wobnum (nowwin) = nobj
c
        end if
c
      else if (optpar(1)(1:1) .eq. 'R') then
c
        optpar (1) = 'ramp_2d'
c
c ... RAMP_2D
c
        if (linter) then
c
          if (numwin .le. nought) then
            call errcon ('No graphics window')
            goto 11
          end if
c
          if (windim(nowwin) .ne. 2) then
            call errcon ('Not a 2D graphics window')
            goto 11
          end if
c
        end if
c
        if (nopt .gt. 1) then
          file2 = optpar(2)
        else
          call textin (' Filename ?',file2)
          optpar (2) = file2
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
          dops = .true.
        else
          dops = .false.
        end if
c
        nopt = max (nopt, 2)
c
        obname = 'ramp_2D'
c
        if (.not. linter .and. .not. dops) then
          call errcon ('Nothing to do')
          goto 11
        end if
c
        if (.not. dops) then
          call errcon ('Ramp_2D only works in PostScript !')
          goto 11
        end if
c
        dovrml = .false.
        doramp = .true.
        call plot2d (f2,file2,npnt,doramp,
     +               dovrml,lsolid,dops,file3,obname,nobj,ierr,
     +               maxbuf,ibuff1,ibuff2,xdata,bits)
c
        if (linter) then
c
          call xtx_redraw_force ()
c
c ... do window book-keeping
c
          wintyp (nowwin) = plot_2d
          winpnt (nowwin) = npnt
          wobnum (nowwin) = nobj
c
        end if
c
      else if (optpar(1)(1:2) .eq. 'SC') then
c
        optpar (1) = 'scatter_plot'
c
c ... SCATTER_PLOT
c
        if (linter) then
c
          if (numwin .le. nought) then
            call errcon ('No graphics window')
            goto 11
          end if
c
          if (windim(nowwin) .ne. 1) then
            call errcon ('Not a 1D graphics window')
            goto 11
          end if
c
        end if
c
        if (nopt .gt. 1) then
          file1 = optpar(2)
        else
          call textin (' Filename ?',file1)
          optpar (2) = file1
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
          dops = .true.
        else
          dops = .false.
        end if
c
        nopt = max (2, nopt)
c
        obname = 'scatter_plot'
c
        if (.not. linter .and. .not. dops) then
          call errcon ('Nothing to do')
          goto 11
        end if
c
        gktype = 'SCAT'
        call plot1d (gktype,
     +    f1,file1,file3,dops,npnt,obname,nobj,ierr)
c
        if (linter) then
c
          call xtx_redraw_force ()
c
c ... do window book-keeping
c
          wintyp (nowwin) = scatter
          winpnt (nowwin) = npnt
          wobnum (nowwin) = nobj
c
        end if
c
      else if (optpar(1)(1:2) .eq. 'SP') then
c
        optpar (1) = 'spike_plot'
c
c ... SPIKE_PLOT
c
        if (linter) then
c
          if (numwin .le. nought) then
            call errcon ('No graphics window')
            goto 11
          end if
c
          if (windim(nowwin) .ne. 1) then
            call errcon ('Not a 1D graphics window')
            goto 11
          end if
c
        end if
c
        if (nopt .gt. 1) then
          file1 = optpar(2)
        else
          call textin (' Filename ?',file1)
          optpar (2) = file1
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
          dops = .true.
        else
          dops = .false.
        end if
c
        nopt = max (2, nopt)
c
        obname = 'spike_plot'
c
        if (.not. linter .and. .not. dops) then
          call errcon ('Nothing to do')
          goto 11
        end if
c
        gktype = 'SPIKE'
        call plot1d (gktype,
     +    f1,file1,file3,dops,npnt,obname,nobj,ierr)
c
        if (linter) then
c
          call xtx_redraw_force ()
c
c ... do window book-keeping
c
          wintyp (nowwin) = scatter
          winpnt (nowwin) = npnt
          wobnum (nowwin) = nobj
c
        end if
c
      else if (optpar(1)(1:2) .eq. 'BO') then
c
        optpar (1) = 'box_plot'
c
c ... BOX_PLOT
c
        if (linter) then
c
          if (numwin .le. nought) then
            call errcon ('No graphics window')
            goto 11
          end if
c
          if (windim(nowwin) .ne. 1) then
            call errcon ('Not a 1D graphics window')
            goto 11
          end if
c
        end if
c
        if (nopt .gt. 1) then
          file1 = optpar(2)
        else
          call textin (' Filename ?',file1)
          optpar (2) = file1
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
          dops = .true.
        else
          dops = .false.
        end if
c
        nopt = max (2, nopt)
c
        obname = 'box_plot'
c
        if (.not. linter .and. .not. dops) then
          call errcon ('Nothing to do')
          goto 11
        end if
c
        gktype = 'BOX'
        call plot1d (gktype,
     +    f1,file1,file3,dops,npnt,obname,nobj,ierr)
c
        if (linter) then
c
          call xtx_redraw_force ()
c
c ... do window book-keeping
c
          wintyp (nowwin) = scatter
          winpnt (nowwin) = npnt
          wobnum (nowwin) = nobj
c
        end if
c
      else if (optpar(1)(1:2) .eq. 'HI') then
c
        optpar (1) = 'histogram'
c
c ... HISTOGRAM
c
        if (linter) then
c
          if (numwin .le. nought) then
            call errcon ('No graphics window')
            goto 11
          end if
c
          if (windim(nowwin) .ne. 1) then
            call errcon ('Not a 1D graphics window')
            goto 11
          end if
c
        end if
c
        if (nopt .gt. 1) then
          file1 = optpar(2)
        else
          call textin (' Filename ?',file1)
          optpar (2) = file1
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
          dops = .true.
        else
          dops = .false.
        end if
c
        nopt = max (2, nopt)
c
        obname = 'histogram'
c
        if (.not. linter .and. .not. dops) then
          call errcon ('Nothing to do')
          goto 11
        end if
c
        gktype = 'HIST'
        call plot1d (gktype,
     +    f1,file1,file3,dops,npnt,obname,nobj,ierr)
c
        if (linter) then
c
          call xtx_redraw_force ()
c
c ... do window book-keeping
c
          wintyp (nowwin) = histo
          winpnt (nowwin) = npnt
          wobnum (nowwin) = nobj
c
        end if
c
c
      else if (optpar(1)(1:1) .eq. 'P') then
c
        optpar (1) = 'pie_chart'
c
c ... PIE CHART
c
        if (linter) then
c
          if (numwin .le. nought) then
            call errcon ('No graphics window')
            goto 11
          end if
c
          if (windim(nowwin) .ne. 1) then
            call errcon ('Not a 1D graphics window')
            goto 11
          end if
c
        end if
c
        if (nopt .gt. 1) then
          file1 = optpar(2)
        else
          call textin (' Filename ?',file1)
          optpar (2) = file1
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
          dops = .true.
        else
          dops = .false.
        end if
c
        nopt = max (2, nopt)
c
        obname = 'pie'
c
        if (.not. linter .and. .not. dops) then
          call errcon ('Nothing to do')
          goto 11
        end if
c
        gktype = 'PIE'
        call plot1d (gktype,
     +    f1,file1,file3,dops,npnt,obname,nobj,ierr)
c
        if (linter) then
c
          call xtx_redraw_force ()
c
c ... do window book-keeping
c
          wintyp (nowwin) = pie
          winpnt (nowwin) = npnt
          wobnum (nowwin) = nobj
c
        end if
c
      else if (optpar(1)(1:1) .eq. 'T') then
c
        optpar (1) = 'topology'
c
c ... TOPOLOGY
c
        if (linter) then
c
          if (numwin .le. nought) then
            call errcon ('No graphics window')
            goto 11
          end if
c
          if (windim(nowwin) .ne. 2) then
            call errcon ('Not a 2D graphics window')
            goto 11
          end if
c
        end if
c
        if (nopt .gt. 1) then
          file1 = optpar(2)
        else
          call textin (' Filename ?',file1)
          optpar (2) = file1
        end if
c
        if (nopt .gt. 2) then
          file3 = optpar(3)
          dops = .true.
        else
          dops = .false.
        end if
c
        nopt = max (2, nopt)
c
        obname = 'topology'
c
        if (.not. linter .and. .not. dops) then
          call errcon ('Nothing to do')
          goto 11
        end if
c
        call plotop (
     +    f1,file1,file3,dops,obname,nobj,ierr)
c
        if (linter) then
c
          call xtx_redraw_force ()
c
c ... do window book-keeping
c
          wintyp (nowwin) = topo
          winpnt (nowwin) = 0
          wobnum (nowwin) = nobj
c
        end if
c
      else if (optpar(1)(1:1) .eq. 'F') then
c
        if (.not. linter) then
          call errcon ('Non-interactive run; no on-screen graphics')
          goto 11
        end if
c
        optpar (1) = 'focus'
c
c ... FOCUS
c
        if (numwin .le. nought) then
          call errcon ('No windows open')
          goto 11
        end if
c
        if (windim(nowwin) .ne. 1) then
          call errcon ('Not a 1D window')
          goto 11
        end if
c
        if (winpnt(nowwin) .lt. 1) then
          call errcon ('No points plotted')
          goto 11
        end if
c
        if (nopt .gt. 1) then
          read (optpar(2),*,err=5002) ipnt
          goto 5003
        end if
 5002   call ivalin (' Data point number ?',1,ipnt)
        write (optpar(2),*) ipnt
c
 5003   if (ipnt .le. nought .or. ipnt .gt. winpnt(nowwin)) then
          call errcon ('Invalid data point number')
          goto 11
        end if
c
        nopt = 2
c
        scx = 0.01*(winlim(2,nowwin)-winlim(1,nowwin))
        scy = 0.01*(winlim(4,nowwin)-winlim(3,nowwin))
c
        call xtx_set_2d_world (1,1,
     +    x1d(ipnt,nowwin)-scx,x1d(ipnt,nowwin)+scx,
     +    y1d(ipnt,nowwin)-scy,y1d(ipnt,nowwin)+scy)
c
        call xtx_redraw_force ()
c
      else
c
        optpar (1) = '?'
        nopt = 1
c
        call textut (' ERROR - Invalid or ambiguous option :',option)
c
      end if
c
      goto 11
c
c ... all done
c
 9999 continue
      call gkquit
c
      stop
c
      end
c
c ========================================================================
c
      subroutine grfini
c
c --- GRFINI () - initialise graphics-related data structures
c
c --- Gerard J Kleywegt @ 910621
c     Modified @ 910712
c     Modified @ 920411
c     Modified @ 920512
c
      include 'o2d.incl'
c
      integer   defcol,defint
      parameter (defcol=2, defint=1.0e6)
c
      integer i,iw,ir,ic,iwin,icol,irow
c
code ...
c
      do i=1,maxlev
        levint (i) = defint
        levcol (i) = defcol
      end do
c
      nlevel = 1
c
      call xtx_initialise ()
      call xtx_dimensions (iw,ic,ir)
c
      maxwin = min (remwin, iw)
      maxcol = min (remcol, ic)
      maxrow = min (remrow, ir)
c
      call jvalut (' Max nr of graphics windows         :',1,maxwin)
      call jvalut (' Max nr of viewport columns         :',1,maxcol)
      call jvalut (' Max nr of viewport rows            :',1,maxrow)
c
      call xtx_list_dims ()
c
      numwin = nought
      nowwin = nought
      nowcol = nought
      nowrow = nought
      do iwin=1,remwin
        wistat (iwin) = not_open
        wincol (iwin) = nought
        winrow (iwin) = nought
        windim (iwin) = nought
        wintyp (iwin) = no_plot
        winpnt (iwin) = nought
        wicol1 (iwin) = 0
        wicol2 (iwin) = 1
        wobnum (iwin) = 0
        winlim (1,iwin) = 0.0
        winlim (2,iwin) = 1.0
        winlim (3,iwin) = 0.0
        winlim (4,iwin) = 1.0
      end do
c
      return
      end
c
c ========================================================================
c
      subroutine grafix_dims
c
c --- GRAFIX_DIMS () - list dimensions
c
c --- Gerard J Kleywegt @ 911030
c     Modified 920511
c
      include 'o2d.incl'
c
code ...
c
      write (*,*)
      write (*,*) '... GRAPHICS ...'
      call jvalut (' Max number of graphics windows     :',1,maxwin)
      call jvalut (' Max number of viewport columns     :',1,maxcol)
      call jvalut (' Max number of viewport rows        :',1,maxrow)
c
      call xtx_list_dims ()
c
      return
      end
c
c ========================================================================
c
      integer function mcurs (iopt)
c
c --- INT FUNC MCURS (IOPT) - handle popup menus
c
c --- Gerard J Kleywegt @ 910816
c     Modified @ 911022
c     Modified @ 920411
c     Modified @ 920512
c     Modified @ 920810
c
      include 'o2d.incl'
c
      external ic1spc
c
      integer iopt,npspc1,npspc2
c
      save npspc1,npspc2
c
code ...
c
c ... handle option
c
      if (iopt .eq. define_popups) then
c
c ... Manipulate POPUP
c
        call xtx_popup_new ('Manipulate',npspc1)
        call xtx_popup_add (npspc1,(
     +      'Manipulate              %t   %F |'
     +    //'ID data point ^         %x4     |'
     +    //'location ^              %x1  %l |'
     +    //'reset full viewport     %x3     |'
     +    //'next window             %x2  %l |'
     +    //'display one label ^     %x8     |'
     +    //'display all labels      %x9  %l |'
     +    //'draw horizontal line ^  %x5     |'
     +    //'draw vertical line ^    %x6     |'
     +    //'draw cross lines ^      %x7  %l |'
     +    //'quit cursor mode       %x10  %l |'
     +    //'QUIT O2D               %x11     |'
     +  ),ic1spc)
c
c ... Analyse POPUP
c
        call xtx_popup_new ('Analyse',npspc2)
        call xtx_popup_add (npspc2,(
     +      'Analyse                 %t   %F  |'
     +    //'List data  Y > cursor ^ %x12     |'
     +    //'List data  Y < cursor ^ %x13  %l |'
     +    //'quit cursor mode        %x10  %l |'
     +    //'QUIT O2D                %x11     |'
     +  ),ic1spc)
c
      else if (iopt .eq. delete_popups) then
c
        call xtx_popup_delete (npspc1)
        call xtx_popup_delete (npspc2)
c
      end if
c
      mcurs = nought
c
      return
      end
c
c ========================================================================
c
      integer function ic1spc (iopt)
c
c --- INT FUNC IC1SPC (IOPT) - handle 1D spectrum cursor
c
c --- Gerard J Kleywegt @ 911022
c     Modified 920512
c
      include 'o2d.incl'
c
      integer iopt,ierr,iw,ic,ir,i,iwin,nmin,length,leng1
c
      real wx,wy,dmin,dist,scx,scy
c
code ...
c
      ic1spc = nought
c
      if (iopt .lt. 1) then
        call errcon ('IC1SPC - Invalid option')
        return
      end if
c
c ... get cursor position
c
      call xtx_last_2d_cursor (iw,ic,ir,wx,wy)
c
c ... see which data point this is (for some options)
c
      if (iopt .eq. 4 .or. iopt .eq. 8) then
        scx = 1.0/(winlim(2,nowwin)-winlim(1,nowwin))
        scy = 1.0/(winlim(4,nowwin)-winlim(3,nowwin))
        dmin = (scx*(x1d(1,nowwin)-wx))**2 +
     +         (scy*(y1d(1,nowwin)-wy))**2
        nmin = 1
        do i=2,winpnt(nowwin)
          dist = (scx*(x1d(i,nowwin)-wx))**2 +
     +           (scy*(y1d(i,nowwin)-wy))**2
          if (dist .lt. dmin) then
            dmin = dist
            nmin = i
          end if
        end do
        dmin = sqrt (dmin)
      end if
c
c ... IOPT = 1 : location
c
      if (iopt .eq. 1) then
c
        write (o2dout,*) wx,wy
        return
c
c ... IOPT = 2 : next window
c
      else if (iopt .eq. 2) then
c
        iwin = nowwin
        do i=iwin+1,maxwin
          if (wistat(i) .eq. open) then
            iwin = i
            goto 100
          end if
        end do
        do i=1,iwin
          if (wistat(i) .eq. open) then
            iwin = i
            goto 100
          end if
        end do
c
  100   nowwin = iwin
        call ivalut (' Selected window :',1,nowwin)
        nowcol = 1
        nowrow = 1
        call xtx_select_window (iwin,nowcol,nowrow)
c
c ... IOPT = 3 : reset full viewport
c
      else if (iopt .eq. 3) then
c
        call xtx_set_2D_world (1,1,
     +    winlim(1,nowwin),winlim(2,nowwin),
     +    winlim(3,nowwin),winlim(4,nowwin))
        call xtx_redraw_force ()
c
c ... IOPT = 4 : ID data point
c
      else if (iopt .eq. 4) then
c
        if (windim(nowwin) .ne. 1) then
          call errcon ('Not a 1D window')
          return
        end if
c
        if (winpnt(nowwin) .le. 0) then
          call errcon ('No points plotted')
          return
        end if
c
        write (o2dout,6000) nmin,dmin,
     +    x1d(nmin,nowwin),y1d(nmin,nowwin),
     +    labels(nmin,nowwin)(1:leng1(labels(nmin,nowwin)))
c
 6000 format (' Point # ',i6,' (D = ',1pe12.4,')'/
     +  ' X,Y,Label = ',e12.4,',',e12.4,' |',a,'|')
c
c ... IOPT = 5,6,7 : draw line(s) at cursor position
c
      else if (iopt .ge. 5 .and. iopt. le. 7) then
c
        if (wobnum(nowwin) .le. 0) then
          call errcon ('No graphics object in window yet')
          return
        end if
c
        call xtx_object_edit (wobnum(nowwin),ierr)
        if (ierr .ne. no_error) return
c
        call xtx_set_colour (wicol2(nowwin))
c
        if (iopt .eq. 5 .or. iopt .eq. 7) then
          call xtx_2D_move (winlim(1,nowwin),wy)
          call xtx_2D_draw (winlim(2,nowwin),wy)
          call rvalut (' Horizontal line at Y =',1,wy)
        end if
c
        if (iopt .eq. 6 .or. iopt .eq. 7) then
          call xtx_2D_move (wx,winlim(3,nowwin))
          call xtx_2D_draw (wx,winlim(4,nowwin))
          call rvalut (' Vertical   line at X =',1,wx)
        end if
c
        call xtx_object_close (wobnum(nowwin),ierr)
        call xtx_redraw_force ()
c
c ... IOPT = 8 : plot one label
c
      else if (iopt .eq. 8) then
c
        if (windim(nowwin) .ne. 1) then
          call errcon ('Not a 1D window')
          return
        end if
c
        if (wobnum(nowwin) .le. 0) then
          call errcon ('No graphics object in window yet')
          return
        end if
c
        call xtx_object_edit (wobnum(nowwin),ierr)
        if (ierr .ne. no_error) return
c
        call xtx_set_colour (wicol2(nowwin))
c
        call xtx_2D_string (x1d(nmin,nowwin),y1d(nmin,nowwin),
     +    labels(nmin,nowwin))
c
        call xtx_object_close (wobnum(nowwin),ierr)
        call xtx_redraw_force ()
c
c ... IOPT = 9 : plot all labels
c
      else if (iopt .eq. 9) then
c
        if (windim(nowwin) .ne. 1) then
          call errcon ('Not a 1D window')
          return
        end if
c
        if (wobnum(nowwin) .le. 0) then
          call errcon ('No graphics object in window yet')
          return
        end if
c
        call xtx_object_edit (wobnum(nowwin),ierr)
        if (ierr .ne. no_error) return
c
        call xtx_set_colour (wicol2(nowwin))
c
        do i=1,winpnt(nowwin)
          call xtx_2D_string (x1d(i,nowwin),y1d(i,nowwin),
     +      labels(i,nowwin))
        end do
c
        call xtx_object_close (wobnum(nowwin),ierr)
        call xtx_redraw_force ()
c
c ... IOPT = 10 : quit cursor mode
c
      else if (iopt .eq. 10) then
c
        ic1spc = -999
c
c ... IOPT = 11 : quit program
c
      else if (iopt .eq. 11) then
c
        call gkquit
        stop ' DONE-'
c
c ... IOPT = 12 : list data with Y > Y-cursor
c
      else if (iopt .eq. 12) then
c
        write (o2dout,*) 'Data points with Y > ',wy
c
        do i=1,winpnt(nowwin)
          if (y1d(i,nowwin) .gt. wy) then
            write (o2dout,6010) i,
     +       x1d(i,nowwin),y1d(i,nowwin),
     +       labels(i,nowwin)(1:leng1(labels(i,nowwin)))
          end if
        end do
c
 6010 format (' Point # ',i6,' X,Y,Label = ',1pe12.4,',',
     +  e12.4,' |',a,'|')
c
c ... IOPT = 13 : list data with Y < Y-cursor
c
      else if (iopt .eq. 13) then
c
        write (o2dout,*) 'Data points with Y < ',wy
c
        do i=1,winpnt(nowwin)
          if (y1d(i,nowwin) .lt. wy) then
            write (o2dout,6010) i,
     +       x1d(i,nowwin),y1d(i,nowwin),
     +       labels(i,nowwin)(1:leng1(labels(i,nowwin)))
          end if
        end do
c
      else
c
        call errcon ('IC1SPC - Invalid option')
        return
c
      end if
c
      return
      end
c
c ========================================================================
c
      subroutine aplot(x,y,ipt)
c     =========================
c
c ... RB's magic bit
c
c ... join IPT coordinate pairs X,Y of contour crossings, in grid units
c     for m*n grid , coordinates are in range 0,0 to m,n
c
      implicit none
c
      integer ipt,i
c
      real x(ipt),y(ipt),x1,y1,xl,yl,xlo,xhi,ylo,yhi,xd,yd
      real scx,scy
c
      logical lps,linter
c
      common /LPOST/ lps,linter,xlo,xhi,ylo,yhi,scx,scy
c
code ...
c
c ... save 1st and last points before transformation
c
      x1=x(1)
      y1=y(1)
      xl=x(ipt)
      yl=y(ipt)
c
      if (linter) call xtx_2d_polyline (x,y,ipt)
c
      if (lps) then
cc        xd = (x(1)+xlo)*scx
cc        yd = (y(1)+ylo)*scy
        xd = x(1)*scx+xlo
        yd = y(1)*scy+ylo
        call xps_move (xd,yd)
        do i=2,ipt
cc          xd = (x(i)+xlo)*scx
cc          yd = (y(i)+ylo)*scy
          xd = x(i)*scx+xlo
          yd = y(i)*scy+ylo
          call xps_draw (xd,yd)
        end do
        call xps_stroke ()
      end if
c
c ... restore 1st & last points
c
      x(1)=x1
      y(1)=y1
      x(ipt)=xl
      y(ipt)=yl
c
      return
      end
c
c ========================================================================
c
c      subroutine axaoff(xoff,yoff)
c     ============================
c
c      implicit none
c
c      real xoff,yoff
c
c      logical lps
c
c      common /LPOST/ lps
c
code ...
c
c      call xtx_set_2d_offset (xoff,yoff)
c
c      return
c      end
c
c ========================================================================
c
      subroutine acolor (icol)
c     ========================
c
c ... set colour index
c
      implicit none
c
      real xlo,xhi,ylo,yhi,scx,scy
c
      integer icol
c
      logical lps,linter
c
      common /LPOST/ lps,linter,xlo,xhi,ylo,yhi,scx,scy
c
code ...
c
      if (linter) call xtx_set_colour (icol)
c
      if (lps) then
        call xps_colour (icol)
      end if
c
      return
      end
c
c ===============================================================
c
      subroutine gkcntr (func,bits,m,n,cont,lcol,ncont)
c
      implicit none
c
c ... max nr of points in any one contour curve
c
      integer npt
      parameter (npt = 1000)
c
      integer m,n,ncont
c
      real func(0:m-1,0:n-1),cont(0:ncont-1)
      real xp(0:npt-1),yp(0:npt-1)
      real cl,xs,ys
c
      integer lcol(0:ncont-1)
      integer ilevel,i,j,k,l,n1,n2,ipt,i1,j1,idir,ierr,ip,jp
      integer in,jn,ino,jno,ioldir
c      integer ips,jps
c
      logical nbrs(4,0:100,0:100)
      logical lmore,linner
c
      character*1 bits(0:1,0:m-1,0:n-1)
      character line *80
c
code ...
c
      if (m .ge. 100 .or. n .ge. 100) then
        call errcon ('Logical array NBRS dimensioned too small')
        return
      end if
c
      do ilevel=0,ncont-1
c
c ... next level
c
        write (*,*)
        call ivalut (' Level nr :',1,ilevel)
        call ivalut (' Colour   :',1,lcol(ilevel))
        call rvalut (' Level    :',1,cont(ilevel))
c
        call acolor (lcol(ilevel))
        cl = cont(ilevel)
c
c ... label points greater/equal level
c
        n1 = 0
        n2 = 0
c
        call xtx_set_colour (2)
        do i=0,m-1
          do j=0,n-1
            bits (0,i,j) = '0'
            bits (1,i,j) = '0'
            nbrs (1,i,j) = .false.
            nbrs (2,i,j) = .false.
            nbrs (3,i,j) = .false.
            nbrs (4,i,j) = .false.
            if (func(i,j) .ge. cl) then
              bits (0,i,j) = '1'
              n1 = n1 + 1
      if (ilevel .eq. 0) then
              xs = float(i)
              ys = float(j)
              call xtx_2d_symbol (1,xs-0.2,xs+0.2,ys-0.2,ys+0.2)
      end if
            end if
          end do
        end do
c
        call jvalut (' N1 :',1,n1)
        if (n1 .eq. 0) goto 3000
c
c ... label points that are less than level, but neighbours
c     to at least one point greater/equal level
c
        call xtx_set_colour (3)
        do i=0,m-1
          do j=0,n-1
            if (bits(0,i,j) .ne. '0') goto 2000
c
            if (j .gt. 0) then
              if (bits(0,i,j-1) .eq. '1') then
                nbrs (2,i,j) = .true.
                if (bits(0,i,j) .eq. '0') then
                  bits(0,i,j) = '2'
                  n2 = n2 + 1
      if (ilevel .eq. 0) then
              xs = float(i)
              ys = float(j)
              call xtx_2d_symbol (2,xs-0.2,xs+0.2,ys-0.2,ys+0.2)
      end if
                end if
ccc                goto 2000
              end if
            end if
c
            if (j .lt. n-1) then
              if (bits(0,i,j+1) .eq. '1') then
                nbrs (1,i,j) = .true.
                if (bits(0,i,j) .eq. '0') then
                  bits(0,i,j) = '2'
                  n2 = n2 + 1
      if (ilevel .eq. 0) then
              xs = float(i)
              ys = float(j)
              call xtx_2d_symbol (2,xs-0.2,xs+0.2,ys-0.2,ys+0.2)
      end if
                end if
ccc                goto 2000
              end if
            end if
c
            if (i .gt. 0) then
              if (bits(0,i-1,j) .eq. '1') then
                nbrs (3,i,j) = .true.
                if (bits(0,i,j) .eq. '0') then
                  bits(0,i,j) = '2'
                  n2 = n2 + 1
      if (ilevel .eq. 0) then
              xs = float(i)
              ys = float(j)
              call xtx_2d_symbol (2,xs-0.2,xs+0.2,ys-0.2,ys+0.2)
      end if
                end if
ccc                goto 2000
              end if
            end if
c
            if (i .lt. m-1) then
              if (bits(0,i+1,j) .eq. '1') then
                nbrs (4,i,j) = .true.
                if (bits(0,i,j) .eq. '0') then
                  bits(0,i,j) = '2'
                  n2 = n2 + 1
      if (ilevel .eq. 0) then
              xs = float(i)
              ys = float(j)
              call xtx_2d_symbol (2,xs-0.2,xs+0.2,ys-0.2,ys+0.2)
      end if
                end if
ccc                goto 2000
              end if
            end if
c
 2000       continue
          end do
        end do
c
        call jvalut (' N2 :',1,n2)
        if (n2 .eq. 0) goto 3000
c
        call acolor (lcol(ilevel))
c
c ... first handle the borders
c
c==>        do i=0,m-1
c==>          if (nbrs(3,i,0)) then
c==>            idir = 3
c==>            ipt = 0
c==>            call cntget (i,0,idir,func,m,n,cl,xpt(ipt),ypt(ipt),ierr)
c
c
c ... find the next point labelled '2' (use a while loop with LMORE)
c
        lmore = .true.
        i = 0
        j = 0
 2010   continue
c
        if (.not. (nbrs(1,i,j) .or. nbrs(2,i,j) .or.
     +             nbrs(3,i,j) .or. nbrs(4,i,j))) then
          j = j + 1
          if (j .gt. n-1) then
            j = 0
            i = i + 1
            if (i .gt. m-1) then
              lmore = .false.
              goto 2020
            end if
          end if
        end if
c
        if (.not. (nbrs(1,i,j) .or. nbrs(2,i,j) .or.
     +             nbrs(3,i,j) .or. nbrs(4,i,j))) goto 2020
c
c        if (bits(0,i,j) .ne. '2') goto 2020
c
        idir = 0
        if (nbrs(1,i,j)) then
          idir = 1
        else if (nbrs(2,i,j)) then
          idir = 2
        else if (nbrs(3,i,j)) then
          idir = 3
        else if (nbrs(4,i,j)) then
          idir = 4
        end if
c
c        if (j .gt. 0) then
c          if (bits(0,i,j-1) .eq. '1') idir = 1
c        end if
c        if (idir .eq. 0 .and. j .lt. n-1) then
c          if (bits(0,i,j+1) .eq. '1') idir = 2
c        end if
c        if (idir .eq. 0 .and. i .gt. 0) then
c          if (bits(0,i-1,j) .eq. '1') idir = 3
c        end if
c        if (idir .eq. 0 .and. i .lt. n-1) then
c          if (bits(0,i+1,j) .eq. '1') idir = 4
c        end if
c
c ... this shouldn't happen
c
        if (idir .eq. 0) then
          write (*,*) ' ? No nbr at ',i,j,' ?'
          goto 2020
        end if
c
c ... do interpolation and store in polyline arrays XP,YP for
c     plotting later on (or dump in display buffer)
c
        ipt = 0
        i1 = i
        j1 = j
        ioldir = idir
        call cntget (i,j,idir,func,m,n,cl,xs,ys,ierr)
        if (ierr .ne. 0) then
          call errcon ('In CNTGET 1 ??? Should not happen !!!')
          goto 2020
        end if
        xp (ipt) = xs
        yp (ipt) = ys
        ip = 1
        jp = 1
        ino = 1
        jno = 1
        bits (1,i1,j1) = '1'
        nbrs (ioldir,i1,j1) = .false.
c
ccc        write (*,6000) ipt,xp(ipt),yp(ipt),idir,ip,jp
        write (*,6010) ipt,i1,j1,idir,ip,jp
c
c ... another while loop, now with LINNER
c
        linner = .true.
        in = i
        jn = j
c
 2030   continue
c
c        ips = ip
c        jps = jp
c
c ... check for boundary
c
c        if (in .eq. 0 .or. in .eq. (m-1) .or.
c     +      jn .eq. 0 .or. jn .eq. (n-1)) then
c
c ... encountered border - first plot current points, then
c     restart with a new array of contour points
c
c          call aplot (xp(0),yp(0),ipt)
c          
c          call cntbor (in,jn,idir,bits,nbrs,100,m,n,ip,jp,ino,jno,ierr)
c
c          ipt = 0
c          in = ip
c          jn = jp
c          call cntget (in,jn,idir,func,m,n,cl,xp(ipt),yp(ipt),ierr)
c          goto 2040
c        else
          call cntnxt (in,jn,idir,bits,nbrs,100,m,n,ip,jp,ino,jno,ierr)
c        end if
c
        if (ierr .ne. 0) then
          call errcon (' In CNTNXT 1 - Should not happen !!!')
c
c ... FUDGE !!!
c
          linner = .false.
          goto 2040
        end if
c
c ... once you've done 2 points, the very first one is allowed
c     again so as to be able to come full circle
c
         if (ipt .eq. 1) nbrs (ioldir,i1,j1) = .true.
c
c ... check if we're not alternating
c
c        if (ips .eq. 0 .and .ip .eq. 0 .and. jps .eq. 0 .and.
c     +      jp .eq. 0) then
c          call prompt (' Alternate ?')
c          ip = -1
c          jp = -1
c          call cntnxt (in,jn,idir,bits,nbrs,100,m,n,ip,jp,ino,jno,ierr)
c          if (ierr .ne. 0) then
c            call errcon (' In CNTNXT 2 - Should not happen !!!')
c          end if
c        end if
c
        in = in + ip
        jn = jn + jp
        ipt = ipt + 1
c
c ... if max nr of plot points, plot the points and restart array
c
        if (ipt .ge. npt-1) then
          call errcon ('Max nr of points in contour - dump points')
          call aplot (xp(0),yp(0),ipt)
          xp(0) = xp(ipt-1)
          yp(0) = yp(ipt-1)
          ipt = 1
ccc          linner = .false.
        end if
c
        nbrs (idir,in,jn) = .false.
        bits (1,in,jn) = '1'
c
        if (bits(0,in,jn) .ne. '2') then
          call errcon (' Not a 2 point ??? Should not happen !!!')
        end if
c
        call cntget (in,jn,idir,func,m,n,cl,xp(ipt),yp(ipt),ierr)
        if (ierr .ne. 0) then
          call errcon ('In CNTGET 2 ??? Should not happen !!!')
          goto 2040
        end if
c
c ... update directions in X and Y in which we are currently moving
c
        if (xp(ipt) .ge. xp(ipt-1)) then
          ino = 1
        else
          ino = -1
        end if
c
        if (yp(ipt) .ge. yp(ipt-1)) then
          jno = 1
        else
          jno = -1
        end if
c
c        if (ip .ne. 0) ino = ip
c        if (jp .ne. 0) jno = jp
c
ccc        write (*,6000) ipt,xp(ipt),yp(ipt),idir,ip,jp
        write (*,6010) ipt,in,jn,idir,ip,jp
c
c        if (abs(in-37).le.1 .and. abs(jn-29).le.1) then
c          write (*,6000) ipt,xp(ipt),yp(ipt),idir,ip,jp
c        end if
c
c ... have we hit a boundary ?
c
        if (i .eq. 0 .or. i .eq. (m-1) .or. j .eq. 0 .or.
     +      j .eq. (n-1)) then
          call prompt ('Found a boundary !')
          linner = .false.
        end if
c
c ... done ?
c
        if (in .eq. i1 .and. jn .eq. j1 .and. idir .eq. ioldir) then
          call prompt (' Back at the start !')
          linner = .false.
        end if
c
c        if (ipt .gt. 2) then
c          if ( (iabs(in-i1)+iabs(jn-j1)) .eq. 1) then
c            call prompt (' Adding connection back to first point')
c            linner = .false.
c            ipt = ipt + 1
c            xp (ipt) = xs
c            yp (ipt) = ys
c        write (*,6000) ipt,xp(ipt),yp(ipt),idir,ip,jp
c        write (*,6010) ipt,in,jn,idir,ip,jp
c          end if
c        end if
c
c ... we're staying in the same point ???
c
ccc        if ( abs(xp(ipt)-xp(ipt-1)) .le. 0.001 .and.
ccc     +       abs(yp(ipt)-yp(ipt-1)) .le. 0.001) linner = .false.
c
c ... end of while loop with LINNER
c
 2040   continue
        if (linner) goto 2030
c
        nbrs (ioldir,i1,j1) = .false.
c
c ... plot & label the points
c
        call aplot (xp(0),yp(0),ipt+1)
c
c ... label points on screen (debugging)
c
      if (ilevel .eq. 0) then
        do k=0,ipt
          write (line,'(i10)') k+1
          call remspa (line)
          call xtx_2d_string (xp(k),yp(k),line)
        end do
      end if
c
c ... reset the '2' points that were used up in this contour
c
        do k=0,m-1
          do l=0,n-1
            if (bits(1,k,l) .eq. '1') then
              if (.not. (nbrs(1,k,l) .or. nbrs(2,k,l) .or.
     +                   nbrs(3,k,l) .or. nbrs(4,k,l))) then
                bits(0,k,l) = '0'
              end if
              bits(1,k,l) = '0'
            end if
          end do
        end do
c
c ... end of while loop with LMORE (loop over 2D grid)
c
 2020   continue
        if (lmore) goto 2010
c




c
c ... next contour level
c
 3000   continue
      end do
c
 6000 format (' Pt nr ',i8,' X,Y = ',2f10.3,' IDIR, IP, JP = ',3i4)
 6010 format (' Pt nr ',i8,' X,Y = ',2i6,' IDIR, IP, JP = ',3i4)
c
      return
      end
c
c
c
      subroutine cntnxt (i,j,idir,bits,nbrs,nbdim,m,n,ip,jp,
     +                   ino,jno,ierr)
c
      implicit none
c
      integer m,n,i,j,idir,ip,jp,ierr,ino,jno,nbdim
      integer loopy
c
      logical nbrs(4,0:nbdim,0:nbdim)
c
      character*1 bits(0:1,0:m-1,0:n-1)
c
code ...
c
      ierr = -1
      loopy = 1
c
 1000 continue
c
      if (idir .eq. 1) then
        if (ino .ge. 0 .and. (i+1).le.(m-1) ) then
          if (bits(0,i+1,j) .eq. '1' .and.
     +        nbrs(4,i,j)) then
            ip = 0
            jp = 0
            idir = 4
            ierr = 0
          else if ( (j+1) .le. (n-1) ) then
            if (bits(0,i+1,j+1) .eq. '1' .and.
     +          nbrs(1,i+1,j)) then
              ip = 1
              jp = 0
              idir = 1
              ierr = 0
            else if (nbrs(3,i+1,j+1)) then
              ip = 1
              jp = 1
              idir = 3
              ierr = 0
            end if
          end if
        else if ( (i-1) .ge. 0 ) then
          if (bits(0,i-1,j) .eq. '1' .and.
     +        nbrs(3,i,j)) then
            ip = 0
            jp = 0
            idir = 3
            ierr = 0
          else if ( (j+1) .le. (n-1) ) then
            if (bits(0,i-1,j+1) .eq. '1' .and.
     +          nbrs(1,i-1,j)) then
              ip = -1
              jp = 0
              idir = 1
              ierr = 0
            else if (nbrs(4,i-1,j+1)) then
              ip = -1
              jp = 1
              idir = 4
              ierr = 0
            end if
          end if
        end if
      else if (idir .eq. 2) then
        if (ino .ge. 0  .and. (i+1).le.(m-1) ) then
          if (bits(0,i+1,j) .eq. '1' .and.
     +        nbrs(4,i,j)) then
            ip = 0
            jp = 0
            idir = 4
            ierr = 0
          else if ( (j-1) .ge. 0 ) then
            if (bits(0,i+1,j-1) .eq. '1' .and.
     +          nbrs(2,i+1,j)) then
              ip = 1
              jp = 0
              idir = 2
              ierr = 0
            else if (nbrs(3,i+1,j-1)) then
              ip = 1
              jp = -1
              idir = 3
              ierr = 0
            end if
          end if
        else if ( (i-1) .ge. 0 ) then
          if (bits(0,i-1,j) .eq. '1' .and.
     +        nbrs(3,i,j)) then
            ip = 0
            jp = 0
            idir = 3
            ierr = 0
          else if ( (j-1) .ge. 0 ) then
            if (bits(0,i-1,j-1) .eq. '1' .and.
     +          nbrs(2,i-1,j)) then
              ip = -1
              jp = 0
              idir = 2
              ierr = 0
            else if (nbrs(4,i-1,j-1)) then
              ip = -1
              jp = -1
              idir = 4
              ierr = 0
            end if
          end if
        end if
      else if (idir .eq. 3) then
        if (jno .ge. 0 .and. (j+1) .le. (n-1) ) then
          if (bits(0,i,j+1) .eq. '1' .and.
     +        nbrs(1,i,j)) then
            ip = 0
            jp = 0
            idir = 1
            ierr = 0
          else if ( (i-1) .ge. 0 ) then
            if (bits(0,i-1,j+1) .eq. '1' .and.
     +          nbrs(3,i,j+1)) then
              ip = 0
              jp = 1
              idir = 3
              ierr = 0
            else if (nbrs(2,i-1,j+1)) then
              ip = -1
              jp = 1
              idir = 2
              ierr = 0
            end if
          end if
        else if ( (j-1) .ge. 0 ) then
          if (bits(0,i,j-1) .eq. '1' .and.
     +        nbrs(2,i,j)) then
            ip = 0
            jp = 0
            idir = 2
            ierr = 0
          else if ( (i-1) .ge. 0 ) then
            if (bits(0,i-1,j-1) .eq. '1' .and.
     +          nbrs(3,i,j-1)) then
              ip = 0
              jp = -1
              idir = 3
              ierr = 0
            else if (nbrs(1,i-1,j-1)) then
              ip = -1
              jp = -1
              idir = 1
              ierr = 0
            end if
          end if
        end if
      else if (idir .eq. 4) then
        if (jno .ge. 0 .and. (j+1) .le. (n-1) ) then
          if (bits(0,i,j+1) .eq. '1' .and.
     +        nbrs(1,i,j)) then
            ip = 0
            jp = 0
            idir = 1
            ierr = 0
          else if ( (i+1) .le. (m-1) ) then
            if (bits(0,i+1,j+1) .eq. '1' .and.
     +          nbrs(4,i,j+1)) then
              ip = 0
              jp = 1
              idir = 4
              ierr = 0
            else if (nbrs(2,i+1,j+1)) then
              ip = 1
              jp = 1
              idir = 2
              ierr = 0
            end if
          end if
        else if ( (j-1) .ge. 0 ) then
          if (bits(0,i,j-1) .eq. '1' .and.
     +        nbrs(2,i,j)) then
            ip = 0
            jp = 0
            idir = 2
            ierr = 0
          else if ( (i+1) .le. (m-1) ) then
            if (bits(0,i+1,j-1) .eq. '1' .and.
     +         nbrs(4,i,j-1)) then
              ip = 0
              jp = -1
              idir = 4
              ierr = 0
            else if (nbrs(1,i+1,j-1)) then
              ip = 1
              jp = -1
              idir = 1
              ierr = 0
            end if
          end if
        end if
      end if
c
c ... don't visit the same point twice (unless you stand still
c     or get back to the starting point)
c
ccc      if (.not. (ip.eq.0 .and. jp.eq.0)) then
ccc        if (bits(0,i+ip,j+jp) .ne. '2' .and. loopy) then
        if (ierr .ne. 0 .and. loopy .eq. 1) then
          call prompt (' Yeuch - loopy 1')
          if (idir .eq. 3) ino = -ino
          if (idir .eq. 4) ino = -ino
          if (idir .eq. 1) jno = -jno
          if (idir .eq. 2) jno = -jno
          loopy = 2
          goto 1000
        end if
ccc      end if
c
        if (ierr .ne. 0 .and. loopy .eq. 2) then
          call prompt (' Yeuch - loopy 2')
          if (idir .eq. 1) ino = -ino
          if (idir .eq. 2) ino = -ino
          if (idir .eq. 3) jno = -jno
          if (idir .eq. 4) jno = -jno
          loopy = 3
          goto 1000
        end if
c
        if (ierr .ne. 0) then
          call errcon (' OUCH - no continuation found ???')
        end if
c
c      if (ierr .ne. 0 .and. loopy) then
c        call prompt (' Yeuch - loopy')
c        ip = -ip
c        jp = -jp
c        loopy = .false.
c        goto 1000
c      end if
c
      return
      end
c
c
c
      subroutine cntget (i,j,idir,func,m,n,cl,xs,ys,ierr)
c
      implicit none
c
      integer m,n,i,j,idir,ierr
c
      real func(0:m-1,0:n-1),cl,xs,ys
c
code ...
c
      ierr = 0
c
      xs = i
      ys = j
c
      if (idir .eq. 1) then
        ys = j + (cl-func(i,j))/(func(i,j+1)-func(i,j))
      else if (idir .eq. 2) then
        ys = j - (cl-func(i,j))/(func(i,j-1)-func(i,j))
      else if (idir .eq. 3) then
        xs = i - (cl-func(i,j))/(func(i-1,j)-func(i,j))
      else if (idir .eq. 4) then
        xs = i + (cl-func(i,j))/(func(i+1,j)-func(i,j))
      else
        ierr = -1
      end if
c
c      if (abs(xs-i) .le. 0.001 .and. abs(ys-j) .le. 0.001) then
c        print *,'FISHY ! ',i,j,xs,ys,idir,cl
c        print *,func(i,j),func(i,j-1),func(i,j+1)
c        print *,func(i-1,j),func(i+1,j)
c      end if        
c
c ... this shouldn't happen
c
c      if (abs(xs-i) .ge. 0.999 .or. abs(ys-j) .ge. 0.999) then
c        print *,'Huh ? ',i,j,xs,ys,idir,cl
c        print *,func(i,j),func(i,j-1),func(i,j+1)
c        print *,func(i-1,j),func(i+1,j)
c      end if
c
      return
      end
C
C
C
      SUBROUTINE hpcntr (A,BITS,M,N,CONT,LCOL,NCONT)
C     ==============================================
C
      PARAMETER(IPTMAX=1024)
C
C***************************************
C
C                  contour routine
C
C***************************************
C
C  a is function to be contoured, dimensioned a(m,n) or a(m*n)
C  m = number of rows of a
C  n = number of columns of a
C  bits = logical*1 array dimensioned at least 2*m*n in calling routine
C  cont = array of ncont contour levels
C  lcol = array of ncont contour level colours
c  ncont = number of contour levels
c
C  x,y  arrays of dimension iptmax used in subroutine to store contour l
C
C  idiv = number of line segments between contour level crossings in uni
C     cell.  if >= 1, then routine finds idiv-1 coordinates between cros
C     by linear interpolation on the unit cell.
C  contour output via calls to subroutine aplot (x,y,ipt)
C     x,y are coordinates with 0<=x<m-1 and 0<=y<n-1.
C  subroutine constl(cl) is called to set up plotting style for contour
C
C
C
      DIMENSION IDIR(4), ISIDE(5), IVCT(4), JVCT(4), KVCT(4), LVCT(4),
     *SIDE(4), NVCT(4)
      CHARACTER*1 BITS(*)
      DIMENSION A(*),C(4),CONT(*)
      INTEGER LCOL(*)
      DIMENSION X(IPTMAX),Y(IPTMAX)
      DATA  IDIR/3,4,1,2/, ISIDE/1,2,3,4,1/, IVCT/0,-1,0,1/,
     *JVCT/-1,0,1,0/, SIDE/0.,0.,1.0,1.0/, KVCT/4,1,2,3/, LVCT/2,3,4,1/
C
Cc
Cc
      NX=M-1
      NY=N-1
      NCOUNT=0
      CLOLD=1.E28
      IDIV=1
      MN=M*N
      MN2=2*MN
      DIV=IDIV
C
C loop contour levels cl
      DO 190 ICONT=1,NCONT
      CL=CONT(ICONT)
C
C set colour index for this level
      CALL ACOLOR(LCOL(ICONT))
C
C set up plotting style for this contour level
C     call constl(cl)
C
      DO 10 I=1,MN2
   10 BITS(I)='f'
      NVCT(1)=0
      NVCT(2)=MN
      NVCT(3)=M
      NVCT(4)=MN+1
      IPT=1
      MM=M-1
      NN=N-1
C     search for contour crossing between adjacent column of array a(i,j
      I=0
      J=1
      ISUB=0
      JSUB=0
      IRTN=1
  100 IF (J .GT. N) GO TO 140
  110 IF (I .GE. MM) GO TO 130
      I=I+1
      ISUB=ISUB+1
      JSUB=JSUB+1
      IF (A(ISUB)-CL) 115,600,120
  115 IF (A(ISUB+1)-CL) 110,110,125
  120 IF (A(ISUB+1)-CL) 125,110,110
  125 IF (BITS(JSUB+NVCT(1)).EQ.'t')  GO TO 110
      XSTART=(CL-A(ISUB))/(A(ISUB+1)-A(ISUB))
      YSTART=0
      GO TO 200
  130 I=0
      ISUB=ISUB+1
      JSUB=JSUB+1
      J=J+1
      GO TO 100
C     search for contour crossing between adjacent rows of array a(i,j)
  140 I=0
      J=1
      JSUB=0
      ISUB=0
      IRTN=2
  150 IF (J .GT. NN) GO TO 190
  160 IF (I .GE. M) GO TO 180
      I=I+1
      ISUB=ISUB+1
      JSUB=JSUB+1
      IF (A(ISUB)-CL) 165,160,170
  165 IF (A(ISUB+M)-CL) 160,160,175
  170 IF (A(ISUB+M)-CL) 175,160,160
  175 IF (BITS(JSUB+NVCT(2)).EQ.'t')  GO TO 160
      YSTART=(CL-A(ISUB))/(A(ISUB+M)-A(ISUB))
      XSTART=0
      GO TO 200
  180 I=0
      J=J+1
      GO TO 150
C
C     begin following contour line... save indices for return to search
  200 ISAVE=I
      JSAVE=J
      ISUBSV=ISUB
      JSUBSV=JSUB
      XSAVE=XSTART
      YSAVE=YSTART
      X(1)=XSTART+FLOAT(I-1)
      Y(1)=YSTART+FLOAT(J-1)
      IENT=IRTN
      IRS=0
      GO TO 250
C     dump line and follow contour line on opposite side of starting pio
C     when used a second time this entry returns to search
  205 IRS=1
  210 IF (IPT .GT. 1) CALL  APLOT (X,Y,IPT)
      IPT=1
      I=ISAVE
      J=JSAVE
      ISUB=ISUBSV
      JSUB=JSUBSV
      XSTART=XSAVE
      YSTART=YSAVE
      X(1)=XSTART+FLOAT(I-1)
      Y(1)=YSTART+FLOAT(J-1)
      IF (IRS.NE.0) GO TO (110,160), IRTN
      IEXIT=IRTN
      IRS=1
      GO TO 240
C     return from following contour line through a cell
  230 IF (BITS(JSUB+NVCT(IEXIT)).EQ.'t')  GO TO 205
  240 I=I+IVCT(IEXIT)
      J=J+JVCT(IEXIT)
      JSUB=I+(J-1)*M
      ISUB=JSUB
      IENT=IDIR(IEXIT)
  250 BITS(JSUB+NVCT(IENT))='t'
      IF (I.LT.1 .OR. I.GT.MM .OR. J.LT.1 .OR. J.GT.NN)  GO TO 210
C     find contour crossing in new cell
  260 IF (ISUB+1.GT.MN .OR. ISUB+M.GT.MN
     1     .OR. ISUB+1+M.GT.MN)  GO TO 210
      C(1)=A(ISUB+1)
      C(2)=A(ISUB)
      C(3)=A(ISUB+M)
      C(4)=A(ISUB+1+M)
      JRTN=1
      ICNT=1
      JCNT=1
      DO 290 IROUND=1,4
      IF (IROUND .EQ. IENT) GO TO 290
      I1=ISIDE(IROUND)
      I2=ISIDE(IROUND+1)
      IF (C(I1)-CL) 270,285,275
  270 IF (C(I2)-CL) 290,290,280
  275 IF (C(I2)-CL) 280,290,290
  280 IEXIT=IROUND
      ICNT=ICNT+1
      GO TO 290
  285 JEXIT=IROUND
      JCNT=JCNT+1
  290 CONTINUE
      GO TO (300,310,700,210), JCNT
  300 GO TO (210,320,210,800), ICNT
  310 GO TO (710,320,210,210), ICNT
  320 GO TO (330,340,350,360), IENT
  330 GO TO (210,410,500,410), IEXIT
  340 GO TO (510,210,510,400), IEXIT
  350 GO TO (500,410,210,410), IEXIT
  360 GO TO (510,400,510,210), IEXIT
C     follow contour line across a cell to a side
  400 XSTART=SIDE(IENT)
  410 XFIN=SIDE(IEXIT)
      XINC=(XFIN-XSTART)/DIV
      XBASE=FLOAT(I-1)
      YBASE=FLOAT(J-1)
      A1=CL-C(2)
      A2=C(1)-C(2)
      A3=C(3)-C(2)
      A4=C(2)-C(1)+C(4)-C(3)
      DO 440 INTERP=1,IDIV
      XSTART=XSTART+XINC
      YSTART=(A1-A2*XSTART)/(A3+A4*XSTART)
      IF (IPT.LT.IPTMAX)  GO TO 430
      CALL APLOT(X, Y, IPT)
      X(1)=X(IPT)
      Y(1)=Y(IPT)
      IPT=1
  430 IPT=IPT+1
      X(IPT)=XBASE+XSTART
      Y(IPT)=YBASE+YSTART
  440 CONTINUE
      GO TO (230,210,615,635), JRTN
  500 YSTART=SIDE(IENT)
C     follow contour line across a cell to a top or bottom
  510 YFIN=SIDE(IEXIT)
      XBASE=FLOAT(I-1)
      YINC=(YFIN-YSTART)/DIV
      YBASE=FLOAT(J-1)
      A1=CL-C(2)
      A2=C(3)-C(2)
      A3=C(1)-C(2)
      A4=C(2)-C(1)+C(4)-C(3)
      DO 540 INTERP=1,IDIV
      YSTART=YSTART+YINC
      XSTART=(A1-A2*YSTART)/(A3+A4*YSTART)
      IF (IPT.LT.IPTMAX) GO TO 530
      CALL APLOT(X, Y, IPT)
      X(1)=X(IPT)
      Y(1)=Y(IPT)
      IPT=1
  530 IPT=IPT+1
      X(IPT)=XBASE+XSTART
      Y(IPT)=YBASE+YSTART
  540 CONTINUE
      GO TO (230,210,615,635), JRTN
C     follow contour line from corner to corner
  600 K1=ISUB-M
      K2=ISUB+1-M
      K3=ISUB+1
      K4=ISUB+1+M
      K5=ISUB+M
      K6=ISUB-1+M
      K7=ISUB-1
      C1=A(K1)
      C2=A(K2)
      C3=A(K3)
      C4=A(K4)
      C5=A(K5)
      C6=A(K6)
      C7=A(K7)
      C8=A(ISUB)
      IF (ISUB.LT.1 .OR. ISUB.GT.MN) GO TO 640
      X(1)=FLOAT(I-1)
      Y(1)=FLOAT(J-1)
      IF (J .EQ. 1 .OR. J .EQ. M)  GO TO 610
      IF (K1.LT.1 .OR. K1.GT.MN)  GO TO 610
      IF (K2.LT.1 .OR. K2.GT.MN)  GO TO 610
      IF (K3.LT.1 .OR. K3.GT.MN)  GO TO 610
      IF (K4.LT.1 .OR. K4.GT.MN)  GO TO 610
      IF (K5.LT.1 .OR. K5.GT.MN)  GO TO 610
      IF (C3.NE.CL)  GO TO 610
      IF (C1 .EQ. CL .AND. C2 .EQ. CL .AND.
     *    C4 .EQ. CL .AND. C5 .EQ. CL)  GO TO 610
      X(2)=X(1)+1.
      Y(2)=Y(1)
      CALL APLOT (X, Y, 2)
      GO TO 620
  610 IF (J .EQ. 1)  GO TO 620
      IF (K1.LT.1 .OR. K1.GT.MN)  GO TO 620
      IF (K2.LT.1 .OR. K2.GT.MN)  GO TO 620
      IF (K3.LT.1 .OR. K3.GT.MN)  GO TO 620
      IF (C2 .NE. CL)  GO TO 620
      IF (C1 .EQ. CL .OR. C3 .EQ. CL)  GO TO 620
      IF (C1 .GT. CL .AND. C3 .GT. CL .OR.
     *    C1 .LT. CL .AND. C3 .LT. CL)  GO TO 620
      C(1)=C2
      C(2)=C1
      C(3)=C8
      C(4)=C3
      J=J-1
      JRTN=3
      IENT=3
      IEXIT=1
      GO TO 500
  615 IF (IPT .GT. 1)  CALL APLOT (X, Y, IPT)
      IPT=1
      J=J+1
      X(1)=FLOAT(I-1)
      Y(1)=FLOAT(J-1)
  620 IF (J .EQ. M .OR. I .EQ. 1)  GO TO 630
      IF (K3.LT.1 .OR. K3.GT.MN)  GO TO 630
      IF (K4.LT.1 .OR. K4.GT.MN)  GO TO 630
      IF (K5.LT.1 .OR. K5.GT.MN)  GO TO 630
      IF (K6.LT.1 .OR. K6.GT.MN)  GO TO 630
      IF (K7.LT.1 .OR. K7.GT.MN)  GO TO 630
      IF (C5 .NE. CL)  GO TO 630
      IF (C3 .EQ. CL .AND. C4 .EQ. CL .AND.
     *    C6 .EQ. CL .AND. C7 .EQ. CL)   GO TO 630
      X(2)=X(1)
      Y(2)=Y(1)+1.
      CALL APLOT (X, Y, 2)
      GO TO 640
  630 IF (J .EQ. M)  GO TO 640
      IF (K3.LT.1 .OR. K3.GT.MN)  GO TO 640
      IF (K4.LT.1 .OR. K4.GT.MN)  GO TO 640
      IF (K5.LT.1 .OR. K5.GT.MN)  GO TO 640
      IF (C4 .NE. CL)  GO TO 640
      IF (C3 .EQ. CL .OR. C5 .EQ. CL)  GO TO 640
      IF (C3 .GT. CL .AND. C5 .GT. CL .OR.
     *    C3 .LT. CL .AND. C5 .LT. CL)  GO TO 640
      C(1)=C3
      C(2)=C8
      C(3)=C5
      C(4)=C4
      JRTN=4
      IENT=1
      IEXIT=3
      GO TO 500
  635 IF (IPT .GT. 1) CALL APLOT (X, Y, IPT)
      IPT=1
      X(1)=FLOAT(I-1)
      Y(1)=FLOAT(J-1)
  640 GO TO (110,160), IRTN
C    follow contour line from side to corner or corners
  700 JRTN=2
      IOPP=IDIR(IENT)
      I1=ISIDE(IOPP)
      I2=ISIDE(IOPP+1)
      IEXIT=IOPP
      C(I1)=C(KVCT(I1))
      C(I2)=C(LVCT(I2))
      GO TO 320
  710 JRTN=2
      IEXIT=JEXIT
      GO TO 320
C     follow contour line through saddle point
  800 IOPP=IDIR(IENT)
      I1=ISIDE(IENT)
      C1=C(I1)
      I2=ISIDE(IENT+1)
      C2=C(I2)
      I3=ISIDE(IOPP)
      C3=C(I3)
      I4=ISIDE(IOPP+1)
      C4=C(I4)
      IF ((C1-CL)/(C1-C2) .EQ. (C4-CL)/(C4-C3))  GO TO 820
      IF ((C1-CL)/(C1-C4) .GT. (C2-CL)/(C2-C3))  GO TO 810
      IEXIT=I4
      GO TO 320
  810 IEXIT=I2
      GO TO 320
  820 C(I3)=C(I2)
      C(I4)=C(I1)
      IEXIT=I3
      GO TO 320
C
C
  190 CONTINUE
      RETURN
      END
c
c ========================================================================
c
      subroutine plot1d (type,iunit,filnam,psfile,dops,npnt,
     +                   obname,nobj,ierr)
c
c --- PLOT1D - plot 1D graph
c
c --- Gerard J Kleywegt @ 920512
c     Modified @ 920717
c
      include 'o2d.incl'
c
      integer iunit,npnt,nobj,ierr,i,icol,mtyp,junit,ndata
      integer idum(maxpnt),nlabel,iprev,inew,length,next,leng1
c
      real xbuf(maxpnt),ybuf(maxpnt)
      real xlo,xint,ylo,yint,labdat(3,maxtxt)
      real dx,dy,xmin,xmax,ymin,ymax,gx1,gx2,gdx,gy1,gy2,gdy
      real x,y,sum,slope,icept,boxlow,boxmax,boxbin
      real cdvert,cdymax,cdgrey,rxmin,rxmax,rymin,rymax,x1,x2
c
      logical lnew,dops,lview,lnext,first,lmore,lnaive,lgrid,lpie
      logical lhide,ldash,linfit,lstat
c
      character obname*(*),filnam*(*),key*6,line*250,format*80
      character xlabel*80,ylabel*80,type*(*)
      character odata(4)*40,psfile*(*)
      character labtxt(maxtxt)*80
c
code ...
c
      lnaive = .false.
      lpie = (type(1:1) .eq. 'P')
      goto 6501
c
      entry naive (type,iunit,filnam,psfile,dops,npnt,
     +             obname,nobj,ierr)
      lnaive = .true.
      lpie = .false.
c
 6501 continue
      linfit = .false.
      lstat = .false.
      lhide = .false.
      ldash = .true.
      call textut (' Plot type :',type)
      call prompt (' Initialising ...')
      close (iunit)
      call xopxoa (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
      ierr = -1
c
      lnext = .false.
      lmore = .false.
      first = .true.
      lgrid = .false.
c
c ... initialise data
c
      npnt = 0
      do i=1,maxpnt
        x1d (i,nowwin) = float (i)
        y1d (i,nowwin) = x1d (i,nowwin)
        if (linter) then
          write (labels (i,nowwin),'(a,i10)') 'Data_pt_',i
          call remspa (labels(i,nowwin))
        end if
      end do
c
cc      if (labels(1,nowwin)(1:9) .ne. 'Data_pt_1') then
cc        do i=1,maxpnt
cc          write (labels (i,nowwin),'(a,i10)') 'Data_pt_',i
cc          call remspa (labels(i,nowwin))
cc        end do
cc      end if
c
      xlabel = 'no_label'
      ylabel = 'no_label'
      dx = 0.0
      dy = 0.0
      mtyp = 2
      lview = .false.
      nlabel = 0
      icol = 1
      boxlow = 0.0
      boxmax = 0.0
      boxbin = -1.0
      cdvert = 0.1
c
      if (linter) then
        call xtx_set_colour (wicol2(nowwin))
        icol = wicol2(nowwin)
      end if
c
      call prompt (' Reading ...')
c
      if (lnaive) then
        junit = iunit
        read (junit,'(a)') line
        call extrop (line,ndata,4,odata,ierr)
        ylabel = odata(1)
        call remspa (odata(2))
        call upcase (odata(2))
        read (odata(3),*) npnt
        format = odata(4)
        call remspa (format)
c
        call textut (' Label  :',ylabel)
        call textut (' Type   :',odata(2))
        call textut (' Format :',format)
        call jvalut (' Number :',1,npnt)
c
        if (odata(2)(1:1) .eq. 'I') then
          read (junit,format) (idum(i),i=1,npnt)
          do i=1,npnt
            y1d (i,nowwin) = float (idum(i))
          end do
        else if (odata(2)(1:1) .eq. 'F') then
          read (junit,format) (y1d(i,nowwin),i=1,npnt)
        else
          read (junit,*) (y1d(i,nowwin),i=1,npnt)
        end if
c
        close (junit)
        goto 999
      end if
c
c ... number of extras (BOXES, LINES)
c
      next = 0
c
   10 read (iunit,'(a6,1x,a)',end=999,err=10) key,line
      call upcase (key)
      if (key(1:1) .eq. '!') goto 10
c
      if (key .eq. 'END   ') then
        goto 999
      else if (key .eq. 'REMARK') then
        call xps_legend (line)
        call textut (' >',line)
      else if (key .eq. 'NPOINT') then
        read (line,*) npnt
      else if (key .eq. 'COLOUR') then
        read (line,*) icol
      else if (key .eq. 'GRID  ') then
        read (line,*) gx1,gx2,gdx,gy1,gy2,gdy
        lgrid = .true.
      else if (key .eq. 'XYVIEW') then
        read (line,*) xmin,xmax,ymin,ymax
        lview = .true.
      else if (key .eq. 'MRKTYP') then
        read (line,*) mtyp
      else if (key .eq. 'HIDE  ') then
        lhide = .true.
      else if (key .eq. 'NODASH') then
        ldash = .false.
      else if (key .eq. 'DASH  ') then
        ldash = .true.
      else if (key .eq. 'LINFIT') then
        linfit = .true.
      else if (key .eq. 'XYSTAT') then
        lstat = .true.
      else if (key .eq. 'NOHIDE') then
        lhide = .false.
      else if (key .eq. 'CDVERT') then
        read (line,*) cdvert
      else if (key .eq. 'MRKSIZ') then
        read (line,*) dx,dy
      else if (key .eq. 'BOXPAR') then
        read (line,*) boxlow,boxbin,boxmax
      else if (key .eq. 'LINE  ') then
        if (next .lt. maxext) then
          next = next + 1
          typext (next) = 'L'
          read (line,*) (extras(i,next),i=1,4)
          call fvalut (' Line :',4,extras(1,next))
        end if
      else if (key .eq. 'BOX   ') then
        if (next .lt. maxext) then
          next = next + 1
          typext (next) = 'B'
          read (line,*) (extras(i,next),i=1,4)
          call fvalut (' Box  :',4,extras(1,next))
        end if
      else if (key .eq. 'CIRCLE') then
        if (next .lt. maxext) then
          next = next + 1
          typext (next) = 'C'
          read (line,*) (extras(i,next),i=1,3)
          call fvalut (' Circle :',3,extras(1,next))
        end if
      else if (key .eq. 'ELLIPS') then
        if (next .lt. maxext) then
          next = next + 1
          typext (next) = 'E'
          read (line,*) (extras(i,next),i=1,4)
          call fvalut (' Ellipse :',4,extras(1,next))
        end if
      else if (key .eq. 'TEXT  ') then
        if (nlabel .lt. maxtxt) then
          nlabel = nlabel + 1
          read (line,*) labdat(1,nlabel),labdat(2,nlabel),
     +      labdat(3,nlabel),labtxt(nlabel)
        else
          call errcon ('Too many text labels')
        end if
      else if (key .eq. 'XLIMIT') then
        read (line,*) xlo,xint
        x1d (1,nowwin) = xlo
        do i=2,npnt
          x1d (i,nowwin) = x1d (i-1,nowwin) + xint
        end do
      else if (key .eq. 'YLIMIT') then
        read (line,*) ylo,yint
        y1d (1,nowwin) = ylo
        do i=2,npnt
          y1d (i,nowwin) = y1d (i-1,nowwin) + yint
        end do
      else if (key .eq. 'XFROMO') then
c
        junit = iunit + 1
        close (junit)
        call xopxoa (junit,line,xinter(),ierr)
        if (ierr .ne. 0) return
c
c ... read and "parse" header of O data block
c
        read (junit,'(a)') line
        call extrop (line,ndata,4,odata,ierr)
        xlabel = odata(1)
        call remspa (odata(2))
        call upcase (odata(2))
        read (odata(3),*) npnt
        format = odata(4)
        call remspa (format)
c
        call textut (' Label  :',xlabel)
        call textut (' Type   :',odata(2))
        call textut (' Format :',format)
        call jvalut (' Number :',1,npnt)
c
        if (odata(2)(1:1) .eq. 'I') then
          read (junit,format) (idum(i),i=1,npnt)
          do i=1,npnt
            x1d (i,nowwin) = float (idum(i))
          end do
        else if (odata(2)(1:1) .eq. 'F') then
          read (junit,format) (x1d(i,nowwin),i=1,npnt)
        else
          read (junit,*) (x1d(i,nowwin),i=1,npnt)
        end if
c
        close (junit)
c
      else if (key .eq. 'YFROMO') then
c
        junit = iunit + 1
        close (junit)
        call xopxoa (junit,line,xinter(),ierr)
        if (ierr .ne. 0) return
c
c ... read and "parse" header of O data block
c
        read (junit,'(a)') line
        call extrop (line,ndata,4,odata,ierr)
        ylabel = odata(1)
        call remspa (odata(2))
        call upcase (odata(2))
        read (odata(3),*) npnt
        format = odata(4)
        call remspa (format)
c
        call textut (' Label  :',ylabel)
        call textut (' Type   :',odata(2))
        call textut (' Format :',format)
        call jvalut (' Number :',1,npnt)
c
        if (odata(2)(1:1) .eq. 'I') then
          read (junit,format) (idum(i),i=1,npnt)
          do i=1,npnt
            y1d (i,nowwin) = float (idum(i))
          end do
        else if (odata(2)(1:1) .eq. 'F') then
          read (junit,format) (y1d(i,nowwin),i=1,npnt)
        else
          read (junit,*) (y1d(i,nowwin),i=1,npnt)
        end if
c
        close (junit)
c
      else if (key .eq. 'LFROMO') then
c
        junit = iunit + 1
        close (junit)
        call xopxoa (junit,line,xinter(),ierr)
        if (ierr .ne. 0) return
c
c ... read and "parse" header of O data block
c
        read (junit,'(a)') line
        call extrop (line,ndata,4,odata,ierr)
        ylabel = odata(1)
        call remspa (odata(2))
        call upcase (odata(2))
        read (odata(3),*) npnt
        format = odata(4)
        call remspa (format)
c
        call textut (' Label  :',ylabel)
        call textut (' Type   :',odata(2))
        call textut (' Format :',format)
        call jvalut (' Number :',1,npnt)
c
        read (junit,format) (labels(i,nowwin),i=1,npnt)
        close (junit)
c
      else if (key .eq. 'XLABEL') then
        xlabel = line
      else if (key .eq. 'YLABEL') then
        ylabel = line
      else if (key .eq. 'LABELS') then
        format = line
        call remspa (format)
        read (iunit,format) (labels(i,nowwin),i=1,npnt)
      else if (key .eq. 'XVALUE') then
        format = line
        call remspa (format)
        if (format(1:1) .eq. '*') then
          call prompt (' Free format read')
          read (iunit,*) (x1d(i,nowwin),i=1,npnt)
        else
          call textut (' Format :',format)
          read (iunit,format) (x1d(i,nowwin),i=1,npnt)
        end if
      else if (key .eq. 'YVALUE') then
        format = line
        call remspa (format)
        if (format(1:1) .eq. '*') then
          call prompt (' Free format read')
          read (iunit,*) (y1d(i,nowwin),i=1,npnt)
        else
          call textut (' Format :',format)
          read (iunit,format) (y1d(i,nowwin),i=1,npnt)
        end if
      else if (key .eq. 'VALUES') then
        format = line
        call remspa (format)
        call remspa (format)
        if (format(1:1) .eq. '*') then
          call prompt (' Free format read')
          read (iunit,*) 
     +      (x1d(i,nowwin),y1d(i,nowwin),i=1,npnt)
        else
          call textut (' Format :',format)
          read (iunit,format)
     +      (x1d(i,nowwin),y1d(i,nowwin),i=1,npnt)
        end if
      else if (key .eq. 'VALABS') then
        format = line
        call remspa (format)
        read (iunit,format) (x1d(i,nowwin),y1d(i,nowwin),
     +    labels(i,nowwin),i=1,npnt)
      else if (key .eq. 'NEXTFN') then
        filnam = line
        lnext = .true.
      else if (key .eq. 'MORE  ') then
        lmore = .true.
        goto 999
c
c ... XCOUNT not used in plotting (only for COUNT_1D_TO_2D)
c
      else if (key .eq. 'XCOUNT') then
        call textut (' Keyword ignored :',key)
c
c ... YCOUNT not used in plotting (only for COUNT_1D_TO_2D)
c
      else if (key .eq. 'YCOUNT') then
        call textut (' Keyword ignored :',key)
c
c ... any other unrecognised keyword generates an error message
c
      else
        call errcon ('Invalid keyword !')
        call textut (' Key   :',key)
        call textut (' Value :',line)
      end if
c
      goto 10
c
c ... done reading
c
  999 continue
      call prompt (' Plotting ...')
      if (.not. lmore) close (iunit)
c
      if (npnt .lt. 2) return
c
c ... fit data to least-squares line if required
c
      if (linfit) then
        call fitlin (npnt,x1d(1,nowwin),y1d(1,nowwin),slope,icept)
      end if
c
c ... print statistics if required
c
      if (lstat) then
        call dostat (npnt,x1d(1,nowwin),y1d(1,nowwin))
      end if
c
c ... get extremes
c
      if (.not. lview) then
        xmin = x1d(1,nowwin)
        xmax = x1d(1,nowwin)
        ymin = y1d(1,nowwin)
        ymax = y1d(1,nowwin)
        do i=1,npnt
          xmin = min (x1d(i,nowwin),xmin)
          xmax = max (x1d(i,nowwin),xmax)
          ymin = min (y1d(i,nowwin),ymin)
          ymax = max (y1d(i,nowwin),ymax)
        end do
      end if
c
      if (lpie) then
        xmin = -70.0
        xmax =  70.0
        ymin = -70.0
        ymax =  70.0
      end if
c
c .. make sure min and max are not identical !
c
ccc      call rlohi (xmin,xmax)
      if (xmin .eq. xmax) then
        xmin = 0.99*xmin
        xmax = 1.01*xmax
        if (xmin .eq. xmax) then
          xmin = -0.01
          xmax =  0.01
        end if
      end if
c
ccc      call rlohi (ymin,ymax)
      if (ymin .eq. ymax) then
        ymin = 0.99*ymin
        ymax = 1.01*ymax
        if (ymin .eq. ymax) then
          ymin = -0.01
          ymax =  0.01
        end if
      end if
c
ccc      if (type(1:3) .eq. 'BOX') then
ccc        if (boxbin .le. 0.0) then
ccc          boxlow = xmin
ccc          boxmax = xmax
ccc          boxbin = 0.1 * (xmax-xmin)
ccc        end if
ccc      end if
c
      if (linter) then
c
        call xtx_set_2d_world (1,1,xmin,xmax,ymin,ymax)
c
        call xtx_labels_define (xlabel,ylabel,'dummy')
c
        if (first) then
          call xtx_object_alloc (nobj,obname,lnew,ierr)
          if (ierr .ne. 0) return
        end if
c
        call xtx_object_edit (nobj,ierr)
c
        call xtx_2d_box (xmin,xmax,ymin,ymax)
        call xtx_2d_move (xmin,0.0)
        call xtx_2d_draw (xmax,0.0)
        call xtx_2d_move (0.0,ymin)
        call xtx_2d_draw (0.0,ymax)
        call xtx_set_colour (icol)
c
        if (lgrid) then
          do x=gx1,gx2,gdx
            call xtx_2d_move (x,gy1)
            call xtx_2d_draw (x,gy2)
          end do
          do y=gy1,gy2,gdy
            call xtx_2d_move (gx1,y)
            call xtx_2d_draw (gx2,y)
          end do
        end if
c
        if (linfit) then
          if (type(1:4) .eq. 'SCAT' .or.
     +        type(1:4) .eq. 'LINE' .or.
     +        type(1:4) .eq. 'SPIK') then
            call xtx_2d_move (xmin,slope*xmin+icept)
            call xtx_2d_draw (xmax,slope*xmax+icept)
          end if
        end if
c
      end if
c
      if (dops) then
c
        if (first) then
          call xps_open  (-1,psfile,'O2D')
          call xps_inquire (rxmin,rxmax,rymin,rymax)
        end if
c 
c ... next line type
c
        if (ldash) call xps_dash ()
c
        if (first) then
          call xps_scale (xmin,xmax,ymin,ymax)
        end if
        call xps_label (xlabel,ylabel)
c
        call xps_hide (lhide)
c
        if (lgrid) then
          call xps_ps_comment ('Here comes your grid')
          do x=gx1,gx2,gdx
            call xps_move (x,gy1)
            call xps_draw (x,gy2)
          end do
          do y=gy1,gy2,gdy
            call xps_move (gx1,y)
            call xps_draw (gx2,y)
          end do
        end if
c
        if (linfit) then
          if (type(1:4) .eq. 'SCAT' .or.
     +        type(1:4) .eq. 'LINE' .or.
     +        type(1:4) .eq. 'SPIK') then
            call xps_move (xmin,slope*xmin+icept)
            call xps_draw (xmax,slope*xmax+icept)
          end if
        end if
c
c ... set colour
c
        call xps_colour (icol)
c
      end if
c
c ... SCAT
c
      if (type(1:4) .eq. 'SCAT') then
        if (dx .le. 0.0) dx = 0.01*(xmax-xmin)
        if (dy .le. 0.0) dy = 0.01*(ymax-ymin)
        if (dops) call xps_ps_comment ('Here comes the scatter plot')
        do i=1,npnt
          if (linter) call xtx_2d_symbol (mtyp,
     +      x1d(i,nowwin)-dx,x1d(i,nowwin)+dx,
     +      y1d(i,nowwin)-dy,y1d(i,nowwin)+dy)
c
          if (dops) call xps_symbol (mtyp,
     +      x1d(i,nowwin)-dx,x1d(i,nowwin)+dx,
     +      y1d(i,nowwin)-dy,y1d(i,nowwin)+dy)
        end do
c
c ... SPIKE
c
      else if (type(1:4) .eq. 'SPIK') then
        if (dx .le. 0.0) dx = 0.01*(xmax-xmin)
        if (dy .le. 0.0) dy = 0.01*(ymax-ymin)
        if (dops) call xps_ps_comment ('Here comes the spike plot')
        do i=1,npnt
          if (linter) then
            call xtx_2d_move (x1d(i,nowwin),0.0)
            call xtx_2d_draw (x1d(i,nowwin),y1d(i,nowwin))
          end if
          if (dops) then
            call xps_move (x1d(i,nowwin),0.0)
            call xps_draw (x1d(i,nowwin),y1d(i,nowwin))
          end if
        end do
c
c ... LINE
c
      else if (type(1:4) .eq. 'LINE') then
ccc        print *,'linter, dops, npnt ',linter, dops, npnt
        if (linter) call xtx_2d_move (x1d(1,nowwin),y1d(1,nowwin))
        if (dops) call xps_ps_comment ('Here comes the line plot')
        if (dops) call xps_move (x1d(1,nowwin),y1d(1,nowwin))
        do i=1,npnt
          if (linter) call xtx_2d_draw (x1d(i,nowwin),y1d(i,nowwin))
          if (dops) call xps_draw (x1d(i,nowwin),y1d(i,nowwin))
        end do
c
c ... CDPLOT
c
      else if (type(1:6) .eq. 'CDPLOT') then
ccc        print *,'linter, dops, npnt ',linter, dops, npnt
        if (linter) call xtx_2d_move (x1d(1,nowwin),y1d(1,nowwin))
        if (dops) then
          call xps_ps_comment ('Here comes the CD plot')
          cdvert = min (1.0, max (0.01, cdvert))
          cdymax = ymin + cdvert * (ymax - ymin)
          x1d(npnt+1,nowwin) = x1d(npnt,nowwin) +
     +      ( (x1d(npnt,nowwin) - x1d(1,nowwin)) / (npnt-1) )
          call fvalut (' CDVERT :',1,cdvert)
          call fvalut (' CDYMAX :',1,cdymax)
        end if
        do i=1,npnt
          if (linter) call xtx_2d_draw (x1d(i,nowwin),y1d(i,nowwin))
          if (dops) then
            cdgrey = 1.0 - ((y1d(i,nowwin) - ymin) / (ymax - ymin))
            call xps_filled_box (x1d(i,nowwin),x1d(i+1,nowwin),
     +        ymin,cdymax,cdgrey,cdgrey,cdgrey)
          end if
        end do
c
c ... add shaded box as legend
c
        if (dops) then
          call xps_inquire (rxmin,rxmax,rymin,rymax)
          call xps_scale (rxmin,rxmax,rymin,rymax)
          rxmin = 350.0
          rxmax = 500.0
          rymin = 50.0
          rymax = 75.0
          do i=1,30
            x1 = rxmin + float(i-1)*5.0
            x2 = x1 + 5.0
            cdgrey = 1.0 - (float(i-1)/30.0)
            call xps_filled_box (x1,x2,rymin,rymax,cdgrey,cdgrey,cdgrey)
          end do
          call xps_colour (0)
          call xps_move (rxmin,rymin)
          call xps_draw (rxmax,rymin)
          call xps_draw (rxmax,rymax)
          call xps_draw (rxmin,rymax)
          call xps_draw (rxmin,rymin)
          call xps_colour (0)
          write (line,'(f10.3)') ymin
          call pretty (line)
          call xps_text (rxmin,rymin-20.0,12.0,line)
          write (line,'(f10.3)') ymax
          call pretty (line)
          call xps_text (rxmin+125.0,rymin-20.0,12.0,line)
          call xps_scale (xmin,xmax,ymin,ymax)
        end if
c
c ... BOX
c
      else if (type(1:3) .eq. 'BOX') then
        if (dops) call xps_ps_comment ('Here comes the box plot')
        if (boxbin .gt. 0.0) then
          call doboxy (boxlow,boxmax,boxbin,x1d(1,nowwin),
     +      y1d(1,nowwin),xbuf,ybuf,npnt,linter,dops,mtyp,dx,dy)
        else
          call dobox2 (boxlow,boxmax,boxbin,x1d(1,nowwin),
     +      y1d(1,nowwin),xbuf,ybuf,npnt,linter,dops,
     +      mtyp,dx,dy)
        end if
c
c ... HIST
c
      else if (type(1:4) .eq. 'HIST') then
        if (dops) call xps_ps_comment ('Here comes the histogram')
        do i=1,npnt-1
c
          if (linter) then
            call xtx_2d_move (x1d(i,nowwin),   0.0)
            call xtx_2d_draw (x1d(i,nowwin),   y1d(i,nowwin))
            call xtx_2d_draw (x1d(i+1,nowwin), y1d(i,nowwin))
            call xtx_2d_draw (x1d(i+1,nowwin), 0.0)
          end if
c
          if (dops) then
            call xps_move (x1d(i,nowwin),   0.0)
            call xps_draw (x1d(i,nowwin),   y1d(i,nowwin))
            call xps_draw (x1d(i+1,nowwin), y1d(i,nowwin))
            call xps_draw (x1d(i+1,nowwin), 0.0)
          end if
        end do
        npnt = npnt - 1
c
c ... PIE
c
      else if (type(1:3) .eq. 'PIE') then
c
ccc        print *,'PIE'
c
        if (linter) then
          call xtx_2d_move (50.0,0.0)
          do i=1,361
            x = 50.0 * cos (float(i)*degtor)
            y = 50.0 * sin (float(i)*degtor)
            call xtx_2d_draw (x,y)
          end do
          call xtx_2d_move (0.0,0.0)
          call xtx_2d_draw (50.0,0.0)
        end if
c
        if (dops) then
          call xps_ps_comment ('Here comes the pie chart')
          call xps_move (50.0,0.0)
          do i=1,361
            x = 50.0 * cos (float(i)*degtor)
            y = 50.0 * sin (float(i)*degtor)
            call xps_draw (x,y)
          end do
          call xps_move (0.0,0.0)
          call xps_draw (50.0,0.0)
        end if
c
        sum = 0.0
        do i=1,npnt
          sum = sum + y1d(i,nowwin)
        end do
c
ccc        print *,'SUM ',sum
c
        iprev = 0
        do i=1,npnt-1
          if (y1d(i,nowwin) .gt. 0.0) then
c
            if (dops) then
              write (line,*) char(ichar('A')-1+i),' ',
     +          y1d(i,nowwin),' between ',
     +          x1d(i,nowwin),' and ',x1d(i+1,nowwin)
              call pretty (line)
              call xps_legend (line)
c
              print *,line(1:leng1(line))
c
            end if
c
            inew = iprev + nint (360.0*y1d(i,nowwin)/sum)
            inew = min (360, inew)
c
            x = 50.0 * cos (float(inew)*degtor)
            y = 50.0 * sin (float(inew)*degtor)
c
ccc       print *,'i, iprev, inew ',i,iprev,inew,x,y
c
            if (linter) then
              call xtx_2d_move (0.0,0.0)
              call xtx_2d_draw (x,y)
            end if
c
            if (dops) then
              call xps_move (0.0,0.0)
              call xps_draw (x,y)
            end if
c
            x = 60.0 * cos (0.5*float(inew+iprev)*degtor)
            y = 60.0 * sin (0.5*float(inew+iprev)*degtor)
c
            line = char(ichar('A')-1+i)
            if (linter) call xtx_2d_string (x,y,line)
            if (dops) call xps_text (x,y,10.0,line)
c
            iprev = inew
c
          end if
        end do
c
        npnt = npnt - 1
c
      else
c
        call errcon ('Not a valid plot type ? Tell Gerard !')
        call textut (' Type :',type)
c
      end if
c
c ... text labels
c
      if (nlabel .gt. 0) then
        if (linter) call xtx_set_colour (wicol2(nowwin))
        if (dops)   call xps_colour (0)
        if (dops) call xps_ps_comment ('Here come the labels')
        do i=1,nlabel
          if (linter) call xtx_2d_string (labdat(1,i),labdat(2,i),
     +      labtxt(i))
          if (dops) then
            call xps_text (labdat(1,i),labdat(2,i),
     +        labdat(3,i),labtxt(i))
          end if
        end do
        if (linter) call xtx_set_colour (icol)
        if (dops)   call xps_colour (icol)
      end if
c
      first = .false.
c
      if (lnext) then
c
        if (linter) call xtx_object_close (nobj,ierr)
c
        close (iunit)
        call xopxoa (iunit,filnam,xinter(),ierr)
        if (ierr .ne. 0) goto 8372
c
        lnext = .false.
        nlabel = 0
c
        goto 10
c
      end if
c
      if (lmore) then
c
        if (linter) call xtx_object_close (nobj,ierr)
c
        lmore = .false.
        nlabel = 0
c
ccc        if (dops) call xps_dash()
c
        goto 10
c
      end if
c
 8372 continue
c
c ... extras (LINEs, BOXes)
c
      call ivalut (' Nr of extras :',1,next)
      if (next .gt. 0) then
        if (linter) call xtx_set_colour (wicol2(nowwin))
        if (dops) then
          call xps_solid ()
          call xps_colour (0)
          call xps_ps_comment ('Here come the extras')
        end if
        do i=1,next
          if (typext(i) .eq. 'L') then
            if (linter) call xtx_2d_line (extras(1,i),extras(2,i),
     +        extras(3,i),extras(4,i))
            if (dops) then
              call xps_move (extras(1,i),extras(2,i))
              call xps_draw (extras(3,i),extras(4,i))
            end if
          else if (typext(i) .eq. 'B') then
            if (linter) call xtx_2d_box (extras(1,i),extras(2,i),
     +        extras(3,i),extras(4,i))
            if (dops) then
              call xps_move (extras(1,i),extras(3,i))
              call xps_draw (extras(2,i),extras(3,i))
              call xps_draw (extras(2,i),extras(4,i))
              call xps_draw (extras(1,i),extras(4,i))
              call xps_draw (extras(1,i),extras(3,i))
            end if
          else if (typext(i) .eq. 'C') then
            if (linter) call xtx_2d_circle (extras(1,i),extras(2,i),
     +        extras(3,i),extras(4,i))
            if (dops) then
              call xps_circle (extras(1,i),extras(2,i),extras(3,i))
            end if
          else if (typext(i) .eq. 'E') then
            if (linter) call xtx_2d_ellipse (extras(1,i),extras(2,i),
     +        extras(3,i),extras(4,i))
            if (dops) then
              call xps_ellipse (extras(1,i),extras(2,i),
     +                          extras(3,i),extras(4,i))
            end if
          end if
        end do
      end if
c
      if (linter) call xtx_object_close (nobj,ierr)
      if (dops) call xps_close ()
c
      call jvalut (' Number of data points :',1,npnt)
      if (.not. lpie) then
        call rvalut (' Lowest X-value        :',1,xmin)
        call rvalut (' Highest X-value       :',1,xmax)
        call rvalut (' Lowest Y-value        :',1,ymin)
        call rvalut (' Highest Y-value       :',1,ymax)
      end if
c
      if (linter) then
        winlim (1,nowwin) = xmin
        winlim (2,nowwin) = xmax
        winlim (3,nowwin) = ymin
        winlim (4,nowwin) = ymax
      end if
c
      ierr = 0
c
      return
      end
c
c ========================================================================
c
      subroutine plot2d (iunit,filnam,npnt,doramp,
     +                   dovrml,lsolid,dops,filenm,obname,nobj,ierr,
     +                   maxbuf,ibuff1,ibuff2,xdata,bits)
c
c --- PLOT2D - plot 2D contour
c
c --- Gerard J Kleywegt @ 920512
c     Modified @ 920824
c     Modified @ 951020
c
      include 'o2d.incl'
c
      integer maxbuf
      integer ibuff1(maxbuf),ibuff2(maxbuf)
      real xdata(maxbuf)
      character bits(2*maxbuf)*1
c
c      integer ibuff1(maxbuf),ibuff2(maxbuf)
      integer iunit,npnt,nobj,ierr,i,icol,nx,ny,ioff
      integer j,length,next,leng1
c
      real xlo,ylo,xhi,yhi,scx,scy,zmin,zmax
      real rgbi(3),rgbf(3)
c
      logical lnew,dops,lps,ginter,dovrml,lsolid,doramp
c
      character obname*(*),filnam*(*),key*6,line*1024,format*80
c      character xlabel*80,ylabel*80,bits(2*maxbuf)*1,filenm*(*)
      character xlabel*80,ylabel*80,filenm*(*)
      character ramode*80,ratype*80
c
      data rgbi /0.0,0.0,1.0/, rgbf /1.0,0.0,0.0/
c
      common /LPOST/ lps,ginter,xlo,xhi,ylo,yhi,scx,scy
c
code ...
c
      ginter = linter
      close (iunit)
      call xopxoa (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
      ierr = -1
c
c ... initialise data
c
      npnt = 0
      nx = 0
      ny = 0
c
      xlabel = 'no_label'
      ylabel = 'no_label'
c
      zmin = 0.0
      zmax = 0.0
      ramode = 'linear'
      ratype = 'hue'
c
      if (linter) call xtx_set_colour (wicol2(nowwin))
c
c ... number of extras (BOXES, LINES)
c
      next = 0
c
   10 read (iunit,'(a6,1x,a)',end=999,err=10) key,line
      call upcase (key)
      if (key(1:1) .eq. '!') goto 10
c
cc      write (*,'(a6,a)') key,line(1:leng1(line))
c
      if (key .eq. 'END   ') then
        goto 999
      else if (key .eq. 'REMARK') then
        call xps_legend (line)
        call textut (' >',line)
      else if (key .eq. 'XPOINT') then
        read (line,*) nx
      else if (key .eq. 'YPOINT') then
        read (line,*) ny
      else if (key .eq. 'XLIMIT') then
        read (line,*) xlo,xhi
      else if (key .eq. 'YLIMIT') then
        read (line,*) ylo,yhi
      else if (key .eq. 'XLABEL') then
        xlabel = line
      else if (key .eq. 'YLABEL') then
        ylabel = line
      else if (key .eq. 'LINE  ') then
        if (next .lt. maxext) then
          next = next + 1
          typext (next) = 'L'
          read (line,*) (extras(i,next),i=1,4)
          call fvalut (' Line :',4,extras(1,next))
        end if
      else if (key .eq. 'BOX   ') then
        if (next .lt. maxext) then
          next = next + 1
          typext (next) = 'B'
          read (line,*) (extras(i,next),i=1,4)
          call fvalut (' Box  :',4,extras(1,next))
        end if
      else if (key .eq. 'NLEVEL') then
        read (line,*) nlevel
      else if (key .eq. 'LEVELS') then
        read (iunit,*) (levint(i),i=1,nlevel)
        write (line,'(a,1p,50(1x,e10.1))')
     +     'Levels: ',(levint(i),i=1,nlevel)
        call pretty (line)
        call xps_legend (line)
      else if (key .eq. 'COLOUR') then
        read (iunit,*) (levcol(i),i=1,nlevel)
      else if (key .eq. 'RAMLIM') then
        read (line,*) zmin,zmax
        call rlohi (zmin,zmax)
      else if (key .eq. 'RAMODE') then
        ramode = line
      else if (key .eq. 'RAMTYP') then
        ratype = line
      else if (key .eq. 'RAMCOL') then
        read (line,*) (rgbi(i),i=1,3),(rgbf(i),i=1,3)
        do i=1,3
          if (rgbi(i) .gt. 1.0) rgbi(i)=rgbi(i)/256.0
          if (rgbf(i) .gt. 1.0) rgbf(i)=rgbf(i)/256.0
          rgbi(i) = max (0.0, min (1.0, rgbi(i)))
          rgbf(i) = max (0.0, min (1.0, rgbf(i)))
        end do
      else if (key .eq. 'ZVALUE') then
        if (nx*ny .gt. maxbuf) then
          write (*,*)
          call errcon ('Too much data for 2D plot buffer')
          call jvalut (' Nr of data points :',1,(nx*ny))
          call jvalut (' Buffer size       :',1,maxbuf)
          return
        end if
        format = line
        call remspa (format)
        if (format(1:1).eq.'*') then
          print *,'FREE format read'
          read (iunit,*)
     +      ((xdata( (i-1)*nx + j ),j=1,nx),i=1,ny)
        else
          print *,'Formatted read'
          read (iunit,format)
     +      ((xdata( (i-1)*nx + j ),j=1,nx),i=1,ny)
        end if
      else
        call errcon ('Invalid keyword !')
        call textut (' Key   :',key)
        call textut (' Value :',line)
      end if
c
      goto 10
c
c ... done reading
c
  999 continue
      close (iunit)
c
      npnt = nx*ny
      if (nx .lt. 2 .or. ny .lt. 2) return
c
cc      call rvalut (' Z-values   1- 10:',10,xdata)
cc      call rvalut (' Z-values 101-110:',10,xdata(101))
cc      call rvalut (' Z-values 201-210:',10,xdata(201))
cc      call rvalut (' Z-values 301-310:',10,xdata(301))
c
      if (linter) then
c
        call xtx_set_2d_world (1,1,xlo,xhi,ylo,yhi)
        call xtx_labels_define (xlabel,ylabel,'dummy')
        call xtx_object_alloc (nobj,obname,lnew,ierr)
        if (ierr .ne. 0) return
        call xtx_object_edit (nobj,ierr)
c
        call xtx_set_colour (icol)
        call xtx_2d_box (xlo,xhi,ylo,yhi)
c
        call xtx_set_2d_offset (xlo,ylo)
c
      end if
c
      scx = (xhi-xlo)/float(nx-1)
      scy = (yhi-ylo)/float(ny-1)
      if (linter) call xtx_set_2d_scales (scx,scy)
c
      if (dops) then
        call xps_open  (-1,filenm,'O2D')
        call xps_scale (xlo,xhi,ylo,yhi)
        call xps_label (xlabel,ylabel)
c ???
cc        call xps_offset (xlo,ylo)
c 
c ... next line type
c
        call xps_dash ()
      end if
c
      if (dovrml) then
        i = iunit + 1
        call xopxua (i,filenm,linter,ierr)
        if (ierr .ne. 0) return
        call xvrml_open (i,1.0,1.0,1.0)
        call xvrml_colour (0.0,0.0,1.0)
        if (lsolid) then
          call xvrml_face_surf (nx,ny,xdata,xlo,xhi,ylo,yhi)
        else
          call xvrml_wire_surf (nx,ny,xdata,xlo,xhi,ylo,yhi)
        end if
        call xvrml_close ()
        return
      end if
c
      lps = dops
      if (dops) then
        if (doramp) then
c
          lps = .false.
          call hpcntr (xdata,bits,nx,ny,levint,levcol,nlevel)
          lps = .true.
c
          call remspa (ramode)
          call upcase (ramode)
          call remspa (ratype)
          call upcase (ratype)
          call xps_ps_comment ('Here come the ramped boxes')
          if (zmin .eq. 0.0 .and. zmax .eq. 0.0) then
            zmin = xdata (1)
            zmax = xdata (1)
            do i=2,nx*ny
              if (xdata(i).lt.zmin) zmin = xdata(i)
              if (xdata(i).gt.zmax) zmax = xdata(i)
            end do
          end if
          call ramp2d (xdata,nx,ny,rgbi,rgbf,ramode,ratype,
     +                 xlo,xhi,ylo,yhi,zmin,zmax,ibuff1,ibuff2)
        else
          call xps_ps_comment ('Here come the contours')
          call hpcntr (xdata,bits,nx,ny,levint,levcol,nlevel)
        end if
      else
        call hpcntr (xdata,bits,nx,ny,levint,levcol,nlevel)
      end if
c
c ... extras (LINEs, BOXes)
c
      call ivalut (' Nr of extras :',1,next)
      if (next .gt. 0) then
        if (linter) then
          call xtx_set_2d_offset (0.0,0.0)
          call xtx_set_2d_scales (1.0,1.0)
          call xtx_set_2d_world (1,1,xlo,xhi,ylo,yhi)
          call xtx_set_colour (wicol2(nowwin))
        end if
        if (dops) call xps_ps_comment ('Here come the extras')
        if (dops) call xps_colour (0)
        do i=1,next
          if (typext(i) .eq. 'L') then
            if (linter) call xtx_2d_line (extras(1,i),extras(2,i),
     +        extras(3,i),extras(4,i))
            if (dops) then
              call xps_move (extras(1,i),extras(2,i))
              call xps_draw (extras(3,i),extras(4,i))
            end if
          else if (typext(i) .eq. 'B') then
            if (linter) call xtx_2d_box (extras(1,i),extras(2,i),
     +        extras(3,i),extras(4,i))
            if (dops) then
              call xps_move (extras(1,i),extras(3,i))
              call xps_draw (extras(2,i),extras(3,i))
              call xps_draw (extras(2,i),extras(4,i))
              call xps_draw (extras(1,i),extras(4,i))
              call xps_draw (extras(1,i),extras(3,i))
            end if
          end if
        end do
      end if
c
      if (linter) then
c
        call xtx_set_2d_offset (0.0,0.0)
        call xtx_set_2d_scales (1.0,1.0)
c
        call xtx_object_close (nobj,ierr)
c
      end if
c
      call jvalut (' Number of data points :',1,npnt)
      call rvalut (' Lowest X-value        :',1,xlo)
      call rvalut (' Highest X-value       :',1,xhi)
      call rvalut (' Lowest Y-value        :',1,ylo)
      call rvalut (' Highest Y-value       :',1,yhi)
      call jvalut (' Number of levels      :',1,nlevel)
      call rvalut (' Level intensities     :',nlevel,levint)
      call jvalut (' Level colour indices  :',nlevel,levcol)
c
      if (dops) then
        call xps_close ()
      end if
c
      if (linter) then
        winlim (1,nowwin) = xlo
        winlim (2,nowwin) = xhi
        winlim (3,nowwin) = ylo
        winlim (4,nowwin) = yhi
        last2x = nx
        last2y = ny
        last2w = nowwin
      end if
c
      ierr = 0
c
      return
      end
c
c ========================================================================
c
      subroutine plt2cg (iunit,filnam,cgfile,type,maxbuf,xdata,ierr)
c
c --- PLT2CG - convert 1D PLT file to CricketGraph format
c
c --- Gerard J Kleywegt @ 921110
c     Modified @ 921110
c
      include 'o2d.incl'
c
      integer iunit,npnt,kunit,length,idum(maxpnt)
      integer i,ierr,junit,ndata,leng1,nxbin,nybin,maxbuf
      integer j,ix,iy,imax,nok
c
      real xdata(maxbuf)
      real xlo,xint,ylo,yint,xmin,xmax,ymin,ymax,xbin,ybin,zmax
c
      logical lcrick
c
      character filnam*(*),type*(*),key*6,line*128,format*80
      character xlabel*40,ylabel*40
      character odata(4)*40,cgfile*(*)
      character tab*1
c
code ...
c
      close (iunit)
      call xopxoa (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
      kunit = iunit + 29
      close (kunit)
      call xopxua (kunit,cgfile,xinter(),ierr)
      if (ierr .ne. 0) return
c
      ierr = -1
c
      tab = char(9)
      nxbin = -1
      nybin = -1
      lcrick = (type(1:4) .eq. 'CRIC')
      call prompt (' Reading 1D plot file ...')
c
      if (.not. lcrick) then
        call stamp (line)
        write (kunit,'(a,1x,a)') 'REMARK',line(1:leng1(line))
      end if
c
c ... initialise data
c
      npnt = 0
      do i=1,maxpnt
        x1d (i,0) = float (i)
        y1d (i,0) = x1d (i,0)
      end do
c
      xlabel = 'no_X_label'
      ylabel = 'no_Y_label'
c
   10 read (iunit,'(a6,1x,a)',end=999,err=10) key,line
      call upcase (key)
      if (key(1:1) .eq. '!') goto 10
c
      if (key .eq. 'END   ') then
        goto 999
      else if (key .eq. 'REMARK') then
        call textut (' >',line)
        if (.not. lcrick) then
          write (kunit,'(a,1x,a)') 'REMARK',line(1:leng1(line))
        end if
      else if (key .eq. 'NPOINT') then
        read (line,*) npnt
      else if (key .eq. 'XLIMIT') then
        read (line,*) xlo,xint
        x1d (1,0) = xlo
        do i=2,npnt
          x1d (i,0) = x1d (i-1,0) + xint
        end do
      else if (key .eq. 'YLIMIT') then
        read (line,*) ylo,yint
        y1d (1,0) = ylo
        do i=2,npnt
          y1d (i,0) = y1d (i-1,0) + yint
        end do
      else if (key .eq. 'XFROMO') then
c
        junit = iunit + 1
        close (junit)
        call xopxoa (junit,line,xinter(),ierr)
        if (ierr .ne. 0) return
c
c ... read and "parse" header of O data block
c
        read (junit,'(a)') line
        call extrop (line,ndata,4,odata,ierr)
        xlabel = odata(1)
        call remspa (odata(2))
        call upcase (odata(2))
        read (odata(3),*) npnt
        format = odata(4)
        call remspa (format)
c
        if (odata(2)(1:1) .eq. 'I') then
          read (junit,format) (idum(i),i=1,npnt)
          do i=1,npnt
            x1d (i,0) = float (idum(i))
          end do
        else
          read (junit,format) (x1d(i,0),i=1,npnt)
        end if
c
        close (junit)
c
      else if (key .eq. 'YFROMO') then
c
        junit = iunit + 1
        close (junit)
        call xopxoa (junit,line,xinter(),ierr)
        if (ierr .ne. 0) return
c
c ... read and "parse" header of O data block
c
        read (junit,'(a)') line
        call extrop (line,ndata,4,odata,ierr)
        ylabel = odata(1)
        call remspa (odata(2))
        call upcase (odata(2))
        read (odata(3),*) npnt
        format = odata(4)
        call remspa (format)
c
        if (odata(2)(1:1) .eq. 'I') then
          read (junit,format) (idum(i),i=1,npnt)
          do i=1,npnt
            y1d (i,0) = float (idum(i))
          end do
        else
          read (junit,format) (y1d(i,0),i=1,npnt)
        end if
c
        close (junit)
c
      else if (key .eq. 'XLABEL') then
        xlabel = line
      else if (key .eq. 'YLABEL') then
        ylabel = line
      else if (key .eq. 'LABELS') then
        format = line
        call remspa (format)
        read (iunit,format) (labels(i,0),i=1,npnt)
      else if (key .eq. 'XVALUE') then
        format = line
        call remspa (format)
        if (format(1:1) .eq. '*') then
          read (iunit,*) (x1d(i,0),i=1,npnt)
        else
          read (iunit,format) (x1d(i,0),i=1,npnt)
        end if
      else if (key .eq. 'YVALUE') then
        format = line
        call remspa (format)
        if (format(1:1) .eq. '*') then
          read (iunit,*) (y1d(i,0),i=1,npnt)
        else
          read (iunit,format) (y1d(i,0),i=1,npnt)
        end if
      else if (key .eq. 'VALUES') then
        format = line
        call remspa (format)
        if (format(1:1) .eq. '*') then
          read (iunit,*) (x1d(i,0),y1d(i,0),i=1,npnt)
        else
          read (iunit,format) (x1d(i,0),y1d(i,0),i=1,npnt)
        end if
      else if (key .eq. 'VALABS') then
        format = line
        call remspa (format)
        read (iunit,format) (x1d(i,0),y1d(i,0),
     +    labels(i,0),i=1,npnt)
      else if (key .eq. 'NEXTFN') then
        call errcon ('PLT2CG - Key NEXTFN not implemented')
      else if (key .eq. 'MORE  ' .and. lcrick) then
        call jvalut (' Nr of data points :',1,npnt)
        if (lcrick .and. npnt .gt. 2) then
          write (kunit,'(a1)') '*'
          write (kunit,'(a,a1,a)') xlabel(1:leng1(xlabel)),
     +      tab,ylabel(1:leng1(ylabel))
          do i=1,npnt
            write (line,'(f15.5,a1,f15.5)') x1d(i,0),tab,y1d(i,0)
            call remspa (line)
            write (kunit,'(a)') line(1:leng1(line))
          end do
        end if
      else if (key .eq. 'XCOUNT' .and. (.not. lcrick)) then
        read (line,*) xmin,xbin,xmax
        call rlohi (xmin,xmax)
        nxbin = 1 + int((xmax-xmin)/xbin)
      else if (key .eq. 'YCOUNT' .and. (.not. lcrick)) then
        read (line,*) ymin,ybin,ymax
        call rlohi (ymin,ymax)
        nybin = 1 + int((ymax-ymin)/ybin)
c
      else
        call textut (' Keyword ignored :',key)
ccc        call textut (' Value :',line)
      end if
c
      goto 10
c
c ... done reading
c
  999 continue
      close (iunit)
c
      call jvalut (' Nr of data points :',1,npnt)
c
      if (npnt .lt. 2) return
c
ccc      write (kunit,'(a1)') '*'
c
      if (lcrick) then
c
        write (kunit,'(a,a1,a)') xlabel(1:leng1(xlabel)),
     +    tab,ylabel(1:leng1(ylabel))
        do i=1,npnt
          write (line,'(f15.5,a1,f15.5)') x1d(i,0),tab,y1d(i,0)
          call remspa (line)
          write (kunit,'(a)') line(1:leng1(line))
        end do
c
        call prompt (' CricketGraph file written')
        close (kunit)
        return
c
      endif
c
c ... do 1D to 2D conversion by counting in bins
c
      call prompt (' Binning and counting ...')
c
      if (nxbin .le. 0) then
        xmin = x1d(1,0)
        xmax = x1d(1,0)
        do i=1,npnt
          if (x1d(i,0) .lt. xmin) xmin = x1d(i,0)
          if (x1d(i,0) .gt. xmax) xmax = x1d(i,0)
        end do
      end if
      if (nxbin .lt. 10) then
        nxbin = 25
        xbin = 0.04 * (xmax-xmin)
      end if
      write (*,'(a,1p,3e12.4,i6)') ' Xmin, Xmax, Xbin, Nxbin = ',
     +  xmin,xmax,xbin,nxbin
c
      if (nybin .le. 0) then
        ymin = y1d(1,0)
        ymax = y1d(1,0)
        do i=1,npnt
          if (y1d(i,0) .lt. ymin) ymin = y1d(i,0)
          if (y1d(i,0) .gt. ymax) ymax = y1d(i,0)
        end do
      end if
      if (nybin .lt. 10) then
        nybin = 25
        ybin = 0.04 * (ymax-ymin)
      end if
      write (*,'(a,1p,3e12.4,i6)') ' Ymin, Ymax, Ybin, Nybin = ',
     +  ymin,ymax,ybin,nybin
c
      if ((nxbin * nybin) .gt. maxbuf) then
        call errcon ('Buffer too small for 2D array')
        goto 9999
      end if
c
      do i=1,nxbin*nybin
        xdata (i) = 0.0
      end do
c
      nok = 0
      do i=1,npnt
        ix = 1 + int ((x1d(i,0)-xmin)/xbin)
        if (ix .ge. 1 .and. ix .le. nxbin) then
          iy = 1 + int ((y1d(i,0)-ymin)/ybin)
          if (iy .ge. 1 .and. iy .le. nybin) then
            j = (iy-1)*nxbin + ix
            xdata (j) = xdata (j) + 1.0
            nok = nok + 1
          end if
        end if
      end do
      call jvalut (' Nr of points binned  :',1,nok)
c
      zmax = xdata (1)
      do i=1,nxbin*nybin
        if (xdata (i) .gt. zmax) zmax = xdata (i)
      end do
      call jvalut (' Max count in any bin :',1,nint(zmax))
c
c ... create the 2D output file
c
      call prompt (' Writing 2D plot file ...')
c
      write (kunit,'(a)') 'NLEVEL 10'
c
      write (kunit,'(a)') 'LEVELS'
      write (line,'(10(1x,i10))') 
     +  (nint(float(i)*0.09*zmax),i=1,10)
      call pretty (line)
      write (kunit,'(a)') line(1:leng1(line))
c
      write (kunit,'(a)') 'COLOUR'
ccc      write (kunit,'(a)') '1 2 3 3 4 4 5 5 6 6'
      write (kunit,'(a)') '73 123 82 146 401 230 304 221 410 1'
c
      write (kunit,'(a)') 'RAMCOL 1 0 0   0 0 1'
      write (kunit,'(a)') 'RAMTYP HUE'
c
      write (kunit,'(a,1x,a)') 'XLABEL',
     +  xlabel(1:leng1(xlabel))
      write (kunit,'(a,1x,a)') 'YLABEL',
     +  ylabel(1:leng1(ylabel))
      write (kunit,'(a,1x,i10)') 'XPOINT',nxbin
      write (kunit,'(a,1x,i10)') 'YPOINT',nybin
      write (kunit,'(a,2(1x,1p,e15.4))') 'XLIMIT',xmin,xmax
      write (kunit,'(a,2(1x,1p,e15.4))') 'YLIMIT',ymin,ymax
c
      write (kunit,'(a)') 'ZVALUE *'
      do i=1,nxbin*nybin,10
        imax = min (i+9,nxbin*nybin)
        write (line,'(10i10)') (nint(xdata(j)),j=i,imax)
        call pretty (line)
        write (kunit,'(a)') line(1:leng1(line))
      end do
c
      call prompt (' 2D plot file written')
c
c ... done
c
 9999 continue
      close (kunit)
c
      return
      end
c
c ========================================================================
c
      subroutine plotop (iunit,filnam,psfile,dops,
     +                   obname,nobj,ierr)
c
c --- PLOTOP - plot 2D topology graph
c
c --- Gerard J Kleywegt @ 930302
c     Modified @ 
c
      include 'o2d.incl'
c
      integer iunit,nobj,ierr
c
      real xmin,xmax,ymin,ymax,xwid
c
      logical lnew,dops
c
      character obname*(*),filnam*(*),key*6,line*80
      character psfile*(*),mytext*6
c
code ...
c
      close (iunit)
      call xopxoa (iunit,filnam,xinter(),ierr)
      if (ierr .ne. 0) return
c
      ierr = -1
c
      call xtx_set_2d_world (1,1,xmin,xmax,ymin,ymax)
c
      call xtx_labels_define ('topo','topo','dummy')
c
      call xtx_object_alloc (nobj,obname,lnew,ierr)
      if (ierr .ne. 0) return
c
      call xtx_object_edit (nobj,ierr)
c
      call xtx_set_colour (wicol2(nowwin))
c
   10 read (iunit,'(a6,1x,a)',end=999,err=10) key,line
      call upcase (key)
      if (key(1:1) .eq. '!') goto 10
c
      if (key .eq. 'REMARK') then
ccc        call xps_legend (line)
        call textut (' >',line)
        goto 10
      end if
c
      if (key(1:6) .eq. 'XYVIEW') then
c
        read (line,*) xmin,xmax,ymin,ymax
cc        print *,key,xmin,xmax,ymin,ymax
        call xtx_set_2d_world (1,1,xmin,xmax,ymin,ymax)
        call xtx_2d_box (xmin,xmax,ymin,ymax)
        call xtx_2d_move (xmin,0.0)
        call xtx_2d_draw (xmax,0.0)
        call xtx_2d_move (0.0,ymin)
        call xtx_2d_draw (0.0,ymax)
c
        winlim (1,nowwin) = xmin
        winlim (2,nowwin) = xmax
        winlim (3,nowwin) = ymin
        winlim (4,nowwin) = ymax
c
        if (dops) then
          call xps_open  (-1,psfile,'O2D')
          call xps_scale (xmin,xmax,ymin,ymax)
          call xps_label ('topo','topo')
c 
c ... next line type
c
          call xps_dash ()
        end if
c
      else if (key(1:6) .eq. 'ALPHA ') then
c
        read (line,*) xmin,ymin,xmax,ymax,xwid,mytext
        call xtx_2d_move (xmin-0.5*xwid,ymin)
        call xtx_2d_draw (xmin+0.5*xwid,ymin)
        call xtx_2d_draw (xmin+0.5*xwid,ymax)
        call xtx_2d_draw (xmin-0.5*xwid,ymax)
        call xtx_2d_draw (xmin-0.5*xwid,ymin)
        call xtx_2d_move (xmin+0.5*xwid,ymin+0.7*(ymax-ymin))
        call xtx_2d_draw (xmin,ymax)
        call xtx_2d_draw (xmin-0.5*xwid,ymin+0.7*(ymax-ymin))
        call xtx_2d_string (xmin-0.5*xwid,
     +                      ymin+1.1*(ymax-ymin),mytext)
c
        if (dops) then
          call xps_ps_comment ('Here comes a helix')
          call xps_move (xmin-0.5*xwid,ymin)
          call xps_draw (xmin+0.5*xwid,ymin)
          call xps_draw (xmin+0.5*xwid,ymax)
          call xps_draw (xmin-0.5*xwid,ymax)
          call xps_draw (xmin-0.5*xwid,ymin)
          call xps_move (xmin+0.5*xwid,ymin+0.7*(ymax-ymin))
          call xps_draw (xmin,ymax)
          call xps_draw (xmin-0.5*xwid,ymin+0.7*(ymax-ymin))
          call xps_text (xmin-0.5*xwid,
     +                   ymin+1.1*(ymax-ymin),10.0,mytext)
        end if
c
      else if (key(1:6) .eq. 'BETA  ') then
c
        read (line,*) xmin,ymin,xmax,ymax,xwid,mytext
        call xtx_2d_move (xmin-0.3*xwid,ymin)
        call xtx_2d_draw (xmin+0.3*xwid,ymin)
        call xtx_2d_draw (xmin+0.3*xwid,ymin+0.7*(ymax-ymin))
        call xtx_2d_draw (xmin+0.5*xwid,ymin+0.7*(ymax-ymin))
        call xtx_2d_draw (xmin,ymax)
        call xtx_2d_draw (xmin-0.5*xwid,ymin+0.7*(ymax-ymin))
        call xtx_2d_draw (xmin-0.3*xwid,ymin+0.7*(ymax-ymin))
        call xtx_2d_draw (xmin-0.3*xwid,ymin)
        call xtx_2d_string (xmin-0.5*xwid,
     +                      ymin+1.1*(ymax-ymin),mytext)
c
        if (dops) then
          call xps_ps_comment ('Here comes a strand')
          call xps_move (xmin-0.3*xwid,ymin)
          call xps_draw (xmin+0.3*xwid,ymin)
          call xps_draw (xmin+0.3*xwid,ymin+0.7*(ymax-ymin))
          call xps_draw (xmin+0.5*xwid,ymin+0.7*(ymax-ymin))
          call xps_draw (xmin,ymax)
          call xps_draw (xmin-0.5*xwid,ymin+0.7*(ymax-ymin))
          call xps_draw (xmin-0.3*xwid,ymin+0.7*(ymax-ymin))
          call xps_draw (xmin-0.3*xwid,ymin)
          call xps_text (xmin-0.5*xwid,
     +                   ymin+1.1*(ymax-ymin),10.0,mytext)
        end if
c
      else if (key(1:6) .eq. 'END   ') then
c
        goto 999
c
      else
        call errcon ('Invalid keyword !')
        call textut (' Key   :',key)
        call textut (' Value :',line)
      end if
c
      goto 10
c
c ... done reading
c
  999 continue
      close (iunit)
c
      call xtx_object_close (nobj,ierr)
c
      if (dops) then
        call xps_close ()
      end if
c
      return
      end
c
c ========================================================================
c
      subroutine gint1d (xlo,xhi)
c
c --- GINT1D - integrate 1D plot in window NOWWIN between XLO and XHI
c
c --- Gerard J Kleywegt @ 951020
c
      include 'o2d.incl'
c
      integer numint,i
c
      real xlo,xhi,xmin,xmax,integr,totint,sumxy,frac,cog,aver,taver
c
code ...
c
      numint = 0
      integr = 0.0
      totint = 0.0
      xmin = xhi
      xmax = xlo
      sumxy = 0.0
c
      do i=1,winpnt(nowwin)
        totint = totint + y1d(i,nowwin)
        if (x1d(i,nowwin) .ge. xlo .and.
     +      x1d(i,nowwin) .le. xhi) then
          numint = numint + 1
          integr = integr + y1d(i,nowwin)
          xmin = min (xmin,x1d(i,nowwin))
          xmax = max (xmax,x1d(i,nowwin))
          sumxy = sumxy + x1d(i,nowwin)*y1d(i,nowwin)
        end if
      end do
c
      write (*,6000) numint,xlo,xhi,xmin,xmax
c
 6000 format (
     + ' Integrated ',i6,' points between ',1pe12.4,' and ',e12.4/
     + ' Actual integration limits found  ',e12.4,' and ',e12.4)
 6010 format (
     + ' Integral      = ',1pe12.4,' Total curve integral = ',e12.4/
     + ' Average value = ',e12.4  ,' Total curve average  = ',e12.4)
 6020 format (
     + ' Fraction of total = ',f8.4,' Centre of gravity    = ',1pe12.4)
c
      if (numint .gt. 0) then
        aver = integr/float(numint)
        taver = totint/float(winpnt(nowwin))
        write (*,6010) integr,totint,aver,taver
        if (totint .ne. 0.0) then
          frac = integr/totint
          cog  = sumxy / integr
          write (*,6020) frac,cog
        end if
      end if
c
      return
      end
c
c ========================================================================
c
      subroutine gint2d (xlo,xhi,ylo,yhi,maxbuf,xdata)
c
c --- GINT2D - integrate 2D plot in window NOWWIN between XLO-XHI and
c              YLO-YHI
c
c --- Gerard J Kleywegt @ 951020
c
      include 'o2d.incl'
c
      integer maxbuf
      real xdata(maxbuf)
c
      integer numint,i,j,i1,i2,j1,j2,nx,ny,k
c
      real xlo,xhi,ylo,yhi,integr,totint,sumxz,sumyz,frac
      real cogx,cogy,aver,taver,xx,yy,sx,sy
c
code ...
c
      nx = last2x
      ny = last2y
      sx=(winlim(2,nowwin)-winlim(1,nowwin))/float(nx-1)
      sy=(winlim(4,nowwin)-winlim(3,nowwin))/float(ny-1)
c
      i1 = int  ( 1 + (xlo-winlim(1,nowwin))/sx)
      i1 = max(1,min(nx,i1))
      i2 = nint ( 1 + (xhi-winlim(1,nowwin))/sx)
      i2 = max(i1,min(nx,i2))
c
      j1 = int  ( 1 + (ylo-winlim(3,nowwin))/sy)
      j1 = max(1,min(ny,j1))
      j2 = nint ( 1 + (yhi-winlim(3,nowwin))/sy)
      j2 = max(j1,min(ny,j2))
c
      totint = 0.0
      do i=1,(nx*ny)
        totint = totint + xdata(i)
      end do
c
      numint = 0
      integr = 0.0
      sumxz = 0.0
      sumyz = 0.0
      do i=i1,i2
        do j=j1,j2
          numint = numint + 1
          k = (j-1)*nx + i
          integr = integr + xdata(k)
          xx = winlim(1,nowwin) + sx*float(i)
          yy = winlim(3,nowwin) + sy*float(j)
          sumxz = sumxz + xx*xdata(k)
          sumyz = sumyz + yy*xdata(k)
        end do
      end do
c
      write (*,6000) numint,xlo,xhi,winlim(1,nowwin)+sx*float(i1-1),
     +  winlim(1,nowwin)+sx*float(i2-1),ylo,yhi,
     +  winlim(3,nowwin)+sy*float(j1-1),
     +  winlim(3,nowwin)+sy*float(j2-1)
c
 6000 format (
     + ' Integrated ',i6,' points'/
     + ' X range ',1pe12.4,' to ',e12.4,' (',2e12.4,')'/
     + ' Y range ',e12.4,' to ',e12.4,' (',2e12.4,')')
 6010 format (
     + ' Integral      = ',1pe12.4,' Total curve integral = ',e12.4/
     + ' Average value = ',e12.4  ,' Total curve average  = ',e12.4)
 6020 format (
     + ' Fraction of total = ',f8.4/
     + ' Centre of gravity = ',1p,2e12.4)
c
      if (numint .gt. 0) then
        aver = integr/float(numint)
        taver = totint/float(nx*ny)
        write (*,6010) integr,totint,aver,taver
        if (totint .ne. 0.0) then
          frac = integr/totint
          cogx = sumxz /integr
          cogy = sumyz /integr
          write (*,6020) frac,cogx,cogy
        end if
      end if
c
      return
      end
c
c ========================================================================
c
      subroutine gintnd (xlo,xhi,ylo,yhi,maxbuf,xdata)
c
c --- GINTnD - generic integration routine
c
c --- Gerard J Kleywegt @ 951020
c
      include 'o2d.incl'
c
      integer maxbuf
      real xdata(maxbuf)
c
      real xlo,xhi,ylo,yhi
c
code ...
c
      write (*,*)
c
c ... 1D integration
c
      if (windim(nowwin) .eq. 1) then
        call gint1d (xlo,xhi)
        return
      end if
c
c ... check if last 2D window
c
      if (nowwin .ne. last2w) then
        call errcon (' Only data for *last* 2D plot is stored')
        return
      end if
c
c ... 2D integration
c
      call gint2d (xlo,xhi,ylo,yhi,maxbuf,xdata)
c
      return
      end
c
c
c
      subroutine ramp2d (xdata,nx,ny,rgbi,rgbf,ramode,ratype,
     +                   xlo,xhi,ylo,yhi,zmin,zmax,ibuff1,ibuff2)
c
c ... RAMP2D - produce colour-ramped squares in PostScript
c
      implicit none
c
      integer nx,ny
      real xdata (nx,ny)
      integer ibuff1(nx,ny),ibuff2(nx,ny)
c
      real rgbi(3),rgbf(3),hue(3),sat(3),val(3),rgb(3)
      real xlo,xhi,ylo,yhi,zmin,zmax
      real dx,dy,x1,x2,y1,y2,h,delta,zhalf
c
      integer i,j,k,ityp,itop
c
      character ramode*(*),ratype*(*)
c
code ...
c
      if (ramode(1:3) .eq. 'LOG') then
        call prompt (' Ramp mode : LOGARITHMIC')
        call c2log (xdata,nx*ny,zmin,zmax)
      else if (ramode(1:3) .eq. 'ABS') then
        call prompt (' Ramp mode : ABSOLUTE')
        call c2abs (xdata,nx*ny,zmin,zmax)
      else
        call prompt (' Ramp mode : LINEAR')
      end if
c
      if (ratype(1:3) .eq. 'SAT') then
        call prompt (' Ramp type : SATURATION')
        ityp = 2
      else if (ratype(1:3) .eq. 'VAL') then
        call prompt (' Ramp type : VALUE')
        ityp = 3
      else if (ratype(1:3) .eq. 'WHI') then
        call prompt (' Ramp type : WHITE')
        ityp = 4
      else
        call prompt (' Ramp type : HUE')
        ityp = 1
      end if
c
      dx = (xhi-xlo) / float(nx-1)
      dy = (yhi-ylo) / float(ny-1)
c
c ... set up for ramping
c     the algorithm is to change rgb space to hue, sat, value space;
c     then ramp on hue; then turn back into rgb
c 
      call rgbhsv (rgbi(1),rgbi(2),rgbi(3),hue(1),sat(1),val(1))
      call rvalut (' Limit     :',1,zmin)
      call fvalut (' RGB start :',3,rgbi)
ccc      print *,' HSV start',hue(1), sat(1), val(1)
      call rgbhsv (rgbf(1),rgbf(2),rgbf(3),hue(2),sat(2),val(2))
      call rvalut (' Limit     :',1,zmax)
      call fvalut (' RGB end   :',3,rgbf)
ccc      print *,' HSV end  ',hue(2), sat(2), val(2)
c
      if (ityp .eq. 1) then
        if (hue(1) .eq. hue(2)) hue(2) = hue(1) + 180.0
        hue(3) = (hue(2)-hue(1))/100.0
        if (zmin .eq. zmax) then
          delta = 0.0
        else
          delta = (hue(2)-hue(1)) / (zmax-zmin)
        end if
        call fvalut (' HUE start :',1,hue(1))
        call fvalut (' HUE end   :',1,hue(2))
      else if (ityp .eq. 2) then
        if (hue(2) .gt. hue(1)) then
          hue (3) = hue (2)
          sat (1) = 0.0
          sat (2) = 1.0
        else
          hue (3) = hue (1)
          sat (2) = 0.0
          sat (1) = 1.0
        end if
        call fix360 (hue(3))
        val (3) = 1.0
        if (zmin .eq. zmax) then
          delta = 0.0
        else
          delta = (sat(2)-sat(1)) / (zmax-zmin)
        end if
        call fvalut (' SAT start :',1,sat(1))
        call fvalut (' SAT end   :',1,sat(2))
ccc        print *,'HUE, VAL :',hue(3),val(3)
      else if (ityp .eq. 3) then
        if (hue(2) .gt. hue(1)) then
          hue (3) = hue (2)
          val (1) = 0.0
          val (2) = 1.0
        else
          hue (3) = hue (1)
          val (2) = 0.0
          val (1) = 1.0
        end if
        call fix360 (hue(3))
        sat (3) = 0.0
        if (zmin .eq. zmax) then
          delta = 0.0
        else
          delta = (val(2)-val(1)) / (zmax-zmin)
        end if
        call fvalut (' VAL start :',1,val(1))
        call fvalut (' VAL end   :',1,val(2))
ccc        print *,'HUE, SAT :',hue(3),sat(3)
      else if (ityp .eq. 4) then
        zhalf = 0.5*(zmax-zmin)
        if (zmin .eq. zmax) then
          delta = 0.0
        else
          delta = 2.0 / (zmax-zmin)
        end if
      end if
c
      call prompt (' Computing colours ...')
      do i=1,nx
c
        do j=1,ny
c
          if (ityp .eq. 1) then
            if (hue(3) .ne. 0.0) then
              k = int(float(int((xdata(i,j)-zmin)*delta/hue(3)))*hue(3) 
     +            + hue(1))
            else
              k = hue(1)
            end if
            h = k
            if (hue(2) .gt. hue(1)) then
              h = max (hue(1), min (hue(2), h))
            else
              h = max (hue(2), min (hue(1), h))
            end if
            call fix360 (h)
c ---      Convert back from hue, sat, assume full intensity
            call hsvrgb (h,1.0,1.0,rgb(1),rgb(2),rgb(3))
          else if (ityp .eq. 2) then
            sat (3) = sat(1) + (xdata(i,j)-zmin)*delta
            sat (3) = max (0.0, min (1.0, sat(3)))
            sat (3) = float(nint(100.0*sat(3))) / 100.0
ccc            h = hue(3)
            call hsvrgb (hue(3),sat(3),val(3),rgb(1),rgb(2),rgb(3))
ccc         print *,hue(3),sat(3),val(3),rgb(1),rgb(2),rgb(3)
          else if (ityp .eq. 3) then
            val (3) = val(1) + (xdata(i,j)-zmin)*delta
            val (3) = max (0.0, min (1.0, val(3)))
            val (3) = float(nint(100.0*val(3))) / 100.0
ccc            h = hue(3)
            call hsvrgb (hue(3),sat(3),val(3),rgb(1),rgb(2),rgb(3))
ccc         print *,hue(3),sat(3),val(3),rgb(1),rgb(2),rgb(3)
          else if (ityp .eq. 4) then
            if (xdata(i,j) .lt. zhalf) then
              hue(3) = hue(1)
              val(3) = val(1)
              sat(3) = (zhalf-xdata(i,j))*delta
            else
              hue(3) = hue(2)
              val(3) = val(2)
              sat(3) = (xdata(i,j)-zhalf)*delta
            end if
            sat (3) = max (0.0, min (1.0, sat(3)))
            sat (3) = float(nint(50.0*sat(3))) / 50.0
            call hsvrgb (hue(3),sat(3),val(3),rgb(1),rgb(2),rgb(3))
ccc         print *,hue(3),sat(3),val(3),rgb(1),rgb(2),rgb(3)
          end if
c
c ... store colour
c
          call xvrml_encode_rgb (rgb(1),rgb(2),rgb(3),ibuff1(i,j))
          ibuff2 (i,j) = ibuff1 (i,j)
c
        end do
      end do
c
c ... sort colours to find most-often occurring one
c
      call sorcol (nx*ny,ibuff2,itop)
c
      call xvrml_decode_rgb (rgb(1),rgb(2),rgb(3),itop)
      call fvalut (' Most frequent RGB :',3,rgb)
c
c ... paint whole area in most-frequent colour
c
      call xps_filled_box (xlo,xhi,ylo,yhi,rgb(1),rgb(2),rgb(3))
c
c ... now plot the remaining boxes
c
      call prompt (' Writing PostScript boxes ...')
      do i=1,nx
        x1 = xlo + float(i-1)*dx - 0.5*dx
        x2 = x1 + dx
        x1 = max(xlo, min (xhi, x1))
        x2 = max(xlo, min (xhi, x2))
c
        do j=1,ny
c
          if (ibuff1(i,j) .ne. itop) then
c
            y1 = ylo + float(j-1)*dy - 0.5*dy
            y2 = y1 + dy
            y1 = max(ylo, min (yhi, y1))
            y2 = max(ylo, min (yhi, y2))
c
            call xvrml_decode_rgb (rgb(1),rgb(2),rgb(3),ibuff1(i,j))
            call xps_filled_box (x1,x2,y1,y2,rgb(1),rgb(2),rgb(3))
c
          end if
        end do
      end do
      call prompt (' Done')
c
      return
      end
c
c
c
      subroutine c2log (xdata,nd,zmin,zmax)
c
      implicit none
c
      integer nd
      real xdata(nd)
c
      real zmin,zmax
      real q1,q2
c
      integer i
c
code ...
c
      q1 = xdata (1)
      q2 = xdata (1)
      do i=2,nd
        if (xdata(i).lt.q1) q1 = xdata(i)
        if (xdata(i).gt.q2) q2 = xdata(i)
      end do
      q1 = min (q1,zmin)
      q2 = max (q2,zmax)
ccc      print *,'Q1, Q2',q1,q2
ccc      print *,'Z1, Z2',zmin,zmax
      do i=1,nd
        xdata(i) = xdata(i) - q1 + 1.0
        xdata(i) = alog10 (xdata(i))
      end do
      zmin = zmin - q1 + 1.0
      zmin = alog10 (zmin)
      zmax = zmax - q1 + 1.0
      zmax = alog10 (zmax)
ccc      print *,'Z1, Z2 after',zmin,zmax
c
      return
      end
c
c
c
      subroutine c2abs (xdata,nd,zmin,zmax)
c
      implicit none
c
      integer nd
      real xdata(nd)
c
      real zmin,zmax
c
      integer i
c
      logical lzero
c
code ...
c
      lzero = (zmin*zmax .le. 0.0)
      do i=1,nd
        xdata(i) = abs (xdata(i))
      end do
      zmin = abs (zmin)
      zmax = abs (zmax)
      call rlohi (zmin,zmax)
      if (lzero) zmin = 0.0
c
      return
      end
c
c
c
      subroutine sorcol (nd,icol,itop)
c
      implicit none
c
      integer nd
      integer icol(nd)
c
      real xx
c
      integer i,j,n,ncol,itop,ntop,jt
c
code ...
c
      call hsorti (nd,icol)
c
      ncol = 0
      itop = icol(1)
      ntop = 1
      i = 1
c
   10 continue
      do j=i+1,nd
        if (icol(j) .ne. icol(i)) then
          jt = j - 1
          goto 20
        end if
      end do
      jt = nd
   20 continue
      n = jt - i + 1
      ncol = ncol + 1
      if (n .gt. ntop) then
        itop = icol(i)
        ntop = n
      end if
      i = jt + 1
      if (i .lt. nd) goto 10
c
      call jvalut (' Nr of colours :',1,ncol)
      call jvalut (' Most frequent :',1,itop)
      xx = 100.0 * float(ntop) / float(nd)
      call fvalut (' Percentage    :',1,xx)
c
      return
      end
c
c
c
      subroutine fitlin (ndata,x,y,b,a)
c
      implicit none
c
      integer ndata
c
      real x(*),y(*)
      real a,b,ss,sx,sy,sxoss,t,st2
c
      integer i
c
      character line*128
c
code ...
c
      sx = 0.0
      sy = 0.0
      st2 = 0.0
      b = 0.0
c
      do i=1,ndata
        sx = sx + x(i)
        sy = sy + y(i)
      end do
      ss = float(ndata)
c
      sxoss = sx / ss
      do i=1,ndata
        t = x(i) - sxoss
        st2 = st2 + t*t
        b = b + t*y(i)
      end do
c
      b = b / st2
      a = (sy - sx*b) / ss
c
c      call prompt (' ==> Linear fit gives :')
c      call rvalut (' ==> Slope     :',1,b)
c      call rvalut (' ==> Intercept :',1,a)
c
      write (line,'(a,1p,e12.4,a,e12.4)')
     +  'Linear fit: slope ',b,' ... intercept ',a
      call pretty (line)
      call xps_legend (line)
      line = ' '//line
      call prompt (line)
c
      return
      end
c
c
c
      subroutine dostat (ndata,x,y)
c
      implicit none
c
      integer ndata
c
      real x(*),y(*)
      real xave,xsdv,xmin,xmax,xtot,yave,ysdv,ymin,ymax,ytot
      real rmsd,shap,corr,rf1,rf2,x6,x7,x8,x9
c
      character line*128
c
code ...
c
      call xstats (x,ndata,xave,xsdv,xmin,xmax,xtot)
      call xstats (y,ndata,yave,ysdv,ymin,ymax,ytot)
      call xystat (x,y,ndata,
     +             rmsd,shap,corr,rf1,rf2,x6,x7,x8,x9)
c
      write (line,'(a,1p,5e12.4)')
     +  'Stats X-values: ave, sdv, min, max, sum:',
     +  xave,xsdv,xmin,xmax,xtot
      call pretty (line)
      call xps_legend (line)
      line = ' '//line
      call prompt (line)
c
      write (line,'(a,1p,5e12.4)')
     +  'Stats Y-values: ave, sdv, min, max, sum:',
     +  yave,ysdv,ymin,ymax,ytot
      call pretty (line)
      call xps_legend (line)
      line = ' '//line
      call prompt (line)
c
      write (line,'(a,1p,e12.4,0p,a,f8.3)')
     +  'Stats X versus Y: RMSD ',rmsd,' ... CC ',corr
      call pretty (line)
      call xps_legend (line)
      line = ' '//line
      call prompt (line)
c
      return
      end
c
c
c
      subroutine doboxy (boxlow,boxmax,boxbin,xd,yd,
     +    xbuf,ybuf,npnt,linter,dops,mtyp,dx,dy)
c
      implicit none
c
      integer maxbox
      parameter (maxbox=1000)
c
      integer npnt,mtyp,i,n,j,nbox
c
      real xd(npnt),yd(npnt),xbuf(npnt),ybuf(npnt)
      real xbox(maxbox),ybox(maxbox)
      real boxlow,boxbin,boxmax,xnow,dx,dy,y10,y90,y25,y75,y50
      real x1,x2,x3,yave,ysdv,ymin,ymax,ytot,rmsd,shap,corr
      real rf1,rf2,x6,x7,x8,x9,slope,icept
c
      logical linter,dops
c
      character line*256
c
code ...
c
      xnow = boxlow
      nbox = 0
c
      write (line,'(a,i10,a,1p,e12.4)')
     +    'Box plot of ',npnt,' points with bin size ',boxbin
      call pretty (line)
      call xps_legend (line)
      line = ' '//line
      call prompt (line)
c
  100 continue
      n = 0
      do i=1,npnt
        if (xd(i) .ge. xnow .and.
     +      xd(i) .lt. (xnow+boxbin)) then
          n = n + 1
          xbuf (n) = xd(i)
          ybuf (n) = yd(i)
        end if
      end do
c
      if (n .le. 0) goto 200
c
      nbox = nbox + 1
      xbox (nbox) = 0.0
      ybox (nbox) = 0.0
c
      if (n .lt. 10) then
c
        do i=1,n
          xbox (nbox) = xbox (nbox) + xbuf(i)
c
          if (linter) call xtx_2d_symbol (mtyp,
     +      xbuf(i)-dx,xbuf(i)+dx,
     +      ybuf(i)-dy,ybuf(i)+dy)
c
          if (dops) call xps_symbol (mtyp,
     +      xbuf(i)-dx,xbuf(i)+dx,
     +      ybuf(i)-dy,ybuf(i)+dy)
        end do
c
        call xstats (ybuf,n,yave,ysdv,ymin,ymax,ytot)
c
        xbox (nbox) = xbox (nbox) / float(n)
        ybox (nbox) = yave
c
        write (*,6000) 'X-min','X-max','Nr'
        write (*,6010) xnow,xnow+boxbin,n
c
        write (*,6020) 'X-Average','Average','St.dev.',
     +    'Minimum','Maximum','Ave-2*Sdv','Ave+2*Sdv'
        write (*,6030) xbox(nbox),yave,ysdv,ymin,ymax,yave-2.0*ysdv,
     +    yave+2.0*ysdv
c
c ... DO NOT INCLUDE BINS WITH LESS THAN 10 DATA POINTS !
c
        nbox = nbox - 1
c
      else
c
ccc        call rvalut (' Before :',5,ybuf)
c
        call sort2 (n,ybuf,xbuf)
c
ccc        call rvalut (' After  :',5,ybuf)
c
        j = nint(0.1*float(n))
        i = max (1, min (n, j))
        y10 = ybuf (i)
        j = nint(0.9*float(n))
        i = max (1, min (n, j))
        y90 = ybuf (i)
        j = nint(0.25*float(n))
        i = max (1, min (n, j))
        y25 = ybuf (i)
        j = nint(0.75*float(n))
        i = max (1, min (n, j))
        y75 = ybuf (i)
        j = nint(0.5*float(n))
        i = max (1, min (n, j))
        y50 = ybuf (i)
c
        do i=1,n
          xbox (nbox) = xbox (nbox) + xbuf(i)
        end do
c
        call xstats (ybuf,n,yave,ysdv,ymin,ymax,ytot)
c
        xbox (nbox) = xbox (nbox) / float(n)
        ybox (nbox) = yave
c
        write (*,6000) 'X-min','X-max','Nr','Median','25th %',
     +    '75th %','10th %','90th %'
        write (*,6010) xnow,xnow+boxbin,n,y50,y25,y75,y10,y90
c
        write (*,6020) 'X-Average','Average','St.dev.',
     +    'Minimum','Maximum','Ave-2*Sdv','Ave+2*Sdv'
        write (*,6030) xbox(nbox),yave,ysdv,ymin,ymax,yave-2.0*ysdv,
     +    yave+2.0*ysdv
c
 6000   format (/1x,2a12,a6,5a11)
 6010   format (1x,1p,2e12.4,i6,5e11.3)
 6020   format (1x,8x,7a11)
 6030   format (1x,1p,8x,7e11.3)
c
        x1 = xnow+0.05*boxbin
        x2 = xnow+0.5*boxbin
        x3 = xnow+0.95*boxbin
c
        if (linter) then
          call xtx_2d_move (x1,y10)
          call xtx_2d_draw (x3,y10)
          call xtx_2d_move (x2,y10)
          call xtx_2d_draw (x2,y25)
          call xtx_2d_move (x1,y25)
          call xtx_2d_draw (x3,y25)
          call xtx_2d_draw (x3,y75)
          call xtx_2d_draw (x1,y75)
          call xtx_2d_draw (x1,y25)
          call xtx_2d_move (x1,y50)
          call xtx_2d_draw (x3,y50)
          call xtx_2d_move (x2,y75)
          call xtx_2d_draw (x2,y90)
          call xtx_2d_move (x1,y90)
          call xtx_2d_draw (x3,y90)
c
c ... + sign for average
c
          x6 = max (x1 + 0.05*boxbin, xbox(nbox) - 0.2*boxbin)
          x7 = min (x3 - 0.05*boxbin, xbox(nbox) + 0.2*boxbin)
          call xtx_2d_move (x6,yave)
          call xtx_2d_draw (x7,yave)
          call xtx_2d_move (xbox(nbox),yave-0.2*(y75-y25))
          call xtx_2d_draw (xbox(nbox),yave+0.2*(y75-y25))
c
c          call xtx_2d_symbol (1,x2-0.2*boxbin,x2+0.2*boxbin,
c     +      yave-0.2*(y75-y25),yave+0.2*(y75-y25))
c
          do i=1,n
            if (ybuf(i) .lt. y10 .or. ybuf(i) .gt. y90) then
              call xtx_2d_symbol (mtyp,
     +          xbuf(i)-dx,xbuf(i)+dx,
     +          ybuf(i)-dy,ybuf(i)+dy)
            end if
          end do
        end if
c
        if (dops) then
          call xps_move (x1,y10)
          call xps_draw (x3,y10)
          call xps_move (x2,y10)
          call xps_draw (x2,y25)
          call xps_move (x1,y25)
          call xps_draw (x3,y25)
          call xps_draw (x3,y75)
          call xps_draw (x1,y75)
          call xps_draw (x1,y25)
          call xps_move (x1,y50)
          call xps_draw (x3,y50)
          call xps_move (x2,y75)
          call xps_draw (x2,y90)
          call xps_move (x1,y90)
          call xps_draw (x3,y90)
c
c ... + sign for average
c
          x6 = max (x1 + 0.05*boxbin, xbox(nbox) - 0.2*boxbin)
          x7 = min (x3 - 0.05*boxbin, xbox(nbox) + 0.2*boxbin)
          call xps_move (x6,yave)
          call xps_draw (x7,yave)
          call xps_move (xbox(nbox),yave-0.2*(y75-y25))
          call xps_draw (xbox(nbox),yave+0.2*(y75-y25))
c
c          call xps_symbol (1,x2-0.2*boxbin,x2+0.2*boxbin,
c     +      yave-0.2*(y75-y25),yave+0.2*(y75-y25))
c
          do i=1,n
            if (ybuf(i) .lt. y10 .or. ybuf(i) .gt. y90) then
              call xps_symbol (mtyp,
     +          xbuf(i)-dx,xbuf(i)+dx,
     +          ybuf(i)-dy,ybuf(i)+dy)
            end if
          end do
        end if
c
      end if
c
  200 continue
      xnow = xnow + boxbin
      if (xnow .le. boxmax) goto 100
c
      if (nbox .gt. 1) then
c
        write (*,*)
        call xystat (xbox,ybox,nbox,
     +    rmsd,shap,corr,rf1,rf2,x6,x7,x8,x9)
        write (line,'(a,1p,e12.4,0p,a,f8.3,a,i6,a)')
     +    'Stats box-average X versus Y: RMSD ',rmsd,' ... CC ',
     +    corr,' ... ',nbox,' bins'
        call pretty (line)
        call xps_legend (line)
        line = ' '//line
        call prompt (line)
c
        call fitlin (nbox,xbox,ybox,slope,icept)
c
      end if
c
      write (*,*)
c
      return
      end
c
c
c
      subroutine dobox2 (boxlow,boxmax,boxbin,xd,yd,
     +    xbuf,ybuf,npnt,linter,dops,mtyp,dx,dy)
c
      implicit none
c
      integer maxbox
      parameter (maxbox=1000)
c
      integer indlo(maxbox),indhi(maxbox)
      integer npnt,mtyp,i,n,j,k,nbox,numbox
c
      real xd(npnt),yd(npnt),xbuf(npnt),ybuf(npnt)
      real xbox(maxbox),ybox(maxbox)
      real boxlow,boxbin,boxmax,xnow,dx,dy,y10,y90,y25,y75,y50
      real x1,x2,x3,yave,ysdv,ymin,ymax,ytot,rmsd,shap,corr
      real rf1,rf2,x6,x7,x8,x9,slope,icept
c
      logical linter,dops
c
      character line*256
c
code ...
c
      numbox = max (3, nint(-boxbin))
      j = nint(float(npnt)/float(numbox))
      if (j .lt. 3) then
        call errcon ('Not enough points for box plot')
        return
      end if
c
      indlo(1) = 1
      do i=1,numbox
        if (i .gt. 1) indlo(i) = indhi(i-1) + 1
        indhi (i) = nint(float(i)*float(npnt)/float(numbox))
        indhi (i) = max(indlo(i),indhi(i))
        indhi (i) = min(indhi(i),npnt)
      end do
c
      do i=1,npnt
        xbuf (i) = xd(i)
        ybuf (i) = yd(i)
      end do
c
      call sort2 (npnt,xbuf,ybuf)
c
      nbox = 0
c
      write (line,'(a,i10,a,i3,a)')
     +    'Box plot of ',npnt,' points in ',numbox,' bins'
      call pretty (line)
      call xps_legend (line)
      line = ' '//line
      call prompt (line)
c
      do j=1,numbox
c
        if (indlo(j) .gt. npnt) goto 200
c
        n = indhi(j) - indlo(j) + 1
        xnow = xbuf(indlo(j))
        boxbin = xbuf(indhi(j)) - xbuf(indlo(j))
c
ccc        write (*,6100) j,n,xnow,xnow+boxbin
c
ccc 6100 format (' Bin nr ',i3,' Nr of points ',i8/
ccc     +  ' X-lo = ',1pe12.4,' X-hi = ',e12.4)
c
        if (n .le. 0) goto 200
c
        nbox = nbox + 1
        xbox (nbox) = 0.0
c
        if (n .lt. 10) then
c
          do i=1,n
            xbox (nbox) = xbox (nbox) + xbuf(i)
c
            if (linter) call xtx_2d_symbol (mtyp,
     +        xbuf(i)-dx,xbuf(i)+dx,
     +        ybuf(i)-dy,ybuf(i)+dy)
c
            if (dops) call xps_symbol (mtyp,
     +        xbuf(i)-dx,xbuf(i)+dx,
     +        ybuf(i)-dy,ybuf(i)+dy)
          end do
c
          call xstats (ybuf(indlo(j)),n,yave,ysdv,ymin,ymax,ytot)
c
          xbox (nbox) = xbox (nbox) / float(n)
          ybox (nbox) = yave
c
          write (*,6000) 'X-min','X-max','Nr'
          write (*,6010) xnow,xnow+boxbin,n
c
          write (*,6020) 'X-Average','Average','St.dev.',
     +      'Minimum','Maximum','Ave-2*Sdv','Ave+2*Sdv'
          write (*,6030) xbox(nbox),yave,ysdv,ymin,ymax,
     +      yave-2.0*ysdv,yave+2.0*ysdv
c
c ... DO NOT INCLUDE BINS WITH LESS THAN 10 DATA POINTS !
c
          nbox = nbox - 1
c
        else
c
ccc        call rvalut (' Before :',5,ybuf(indlo(j)))
c
        call sort2 (n,ybuf(indlo(j)),xbuf(indlo(j)))
c
ccc        call rvalut (' After  :',5,ybuf(indlo(j)))
c
          k = nint(0.1*float(n))
          i = indlo(j) - 1 + max (1, min (n, k))
          y10 = ybuf (i)
          k = nint(0.9*float(n))
          i = indlo(j) - 1 + max (1, min (n, k))
          y90 = ybuf (i)
          k = nint(0.25*float(n))
          i = indlo(j) - 1 + max (1, min (n, k))
          y25 = ybuf (i)
          k = nint(0.75*float(n))
          i = indlo(j) - 1 + max (1, min (n, k))
          y75 = ybuf (i)
          k = nint(0.5*float(n))
          i = indlo(j) - 1 + max (1, min (n, k))
          y50 = ybuf (i)
c
          do i=indlo(j),indhi(j)
            xbox (nbox) = xbox (nbox) + xbuf(i)
          end do
c
          call xstats (ybuf(indlo(j)),n,yave,ysdv,ymin,ymax,ytot)
c
          xbox (nbox) = xbox (nbox) / float(n)
          ybox (nbox) = yave
c
          write (*,6000) 'X-min','X-max','Nr','Median','25th %',
     +      '75th %','10th %','90th %'
          write (*,6010) xnow,xnow+boxbin,n,y50,y25,y75,y10,y90
c
          write (*,6020) 'X-Average','Average','St.dev.',
     +      'Minimum','Maximum','Ave-2*Sdv','Ave+2*Sdv'
          write (*,6030) xbox(nbox),yave,ysdv,ymin,ymax,
     +      yave-2.0*ysdv,yave+2.0*ysdv
c
 6000   format (/1x,2a12,a6,5a11)
 6010   format (1x,1p,2e12.4,i6,5e11.3)
 6020   format (1x,8x,7a11)
 6030   format (1x,1p,8x,7e11.3)
c
          x1 = xnow+0.03*boxbin
          x2 = xnow+0.5*boxbin
          x3 = xnow+0.97*boxbin
c
          if (linter) then
            call xtx_2d_move (x1,y10)
            call xtx_2d_draw (x3,y10)
            call xtx_2d_move (x2,y10)
            call xtx_2d_draw (x2,y25)
            call xtx_2d_move (x1,y25)
            call xtx_2d_draw (x3,y25)
            call xtx_2d_draw (x3,y75)
            call xtx_2d_draw (x1,y75)
            call xtx_2d_draw (x1,y25)
            call xtx_2d_move (x1,y50)
            call xtx_2d_draw (x3,y50)
            call xtx_2d_move (x2,y75)
            call xtx_2d_draw (x2,y90)
            call xtx_2d_move (x1,y90)
            call xtx_2d_draw (x3,y90)
c
c ... + sign for average
c
            x6 = max (x1 + 0.05*boxbin, xbox(nbox) - 0.2*boxbin)
            x7 = min (x3 - 0.05*boxbin, xbox(nbox) + 0.2*boxbin)
            call xtx_2d_move (x6,yave)
            call xtx_2d_draw (x7,yave)
            call xtx_2d_move (xbox(nbox),yave-0.2*(y75-y25))
            call xtx_2d_draw (xbox(nbox),yave+0.2*(y75-y25))
c
c          call xtx_2d_symbol (1,x2-0.2*boxbin,x2+0.2*boxbin,
c     +      yave-0.2*(y75-y25),yave+0.2*(y75-y25))
c
            do i=1,n
              if (ybuf(indlo(j)-1+i) .lt. y10 .or.
     +            ybuf(indlo(j)-1+i) .gt. y90) then
                call xtx_2d_symbol (mtyp,
     +            xbuf(indlo(j)-1+i)-dx,xbuf(indlo(j)-1+i)+dx,
     +            ybuf(indlo(j)-1+i)-dy,ybuf(indlo(j)-1+i)+dy)
              end if
            end do
          end if
c
          if (dops) then
            call xps_move (x1,y10)
            call xps_draw (x3,y10)
            call xps_move (x2,y10)
            call xps_draw (x2,y25)
            call xps_move (x1,y25)
            call xps_draw (x3,y25)
            call xps_draw (x3,y75)
            call xps_draw (x1,y75)
            call xps_draw (x1,y25)
            call xps_move (x1,y50)
            call xps_draw (x3,y50)
            call xps_move (x2,y75)
            call xps_draw (x2,y90)
            call xps_move (x1,y90)
            call xps_draw (x3,y90)
c
c ... + sign for average
c
            x6 = max (x1 + 0.05*boxbin, xbox(nbox) - 0.2*boxbin)
            x7 = min (x3 - 0.05*boxbin, xbox(nbox) + 0.2*boxbin)
            call xps_move (x6,yave)
            call xps_draw (x7,yave)
            call xps_move (xbox(nbox),yave-0.2*(y75-y25))
            call xps_draw (xbox(nbox),yave+0.2*(y75-y25))
c
c          call xps_symbol (1,x2-0.2*boxbin,x2+0.2*boxbin,
c     +      yave-0.2*(y75-y25),yave+0.2*(y75-y25))
c
            do i=1,n
              if (ybuf(indlo(j)-1+i) .lt. y10 .or.
     +            ybuf(indlo(j)-1+i) .gt. y90) then
                call xps_symbol (mtyp,
     +            xbuf(indlo(j)-1+i)-dx,xbuf(indlo(j)-1+i)+dx,
     +            ybuf(indlo(j)-1+i)-dy,ybuf(indlo(j)-1+i)+dy)
              end if
            end do
          end if
c
        end if
c
  200   continue
c
      end do
c
      if (nbox .gt. 1) then
c
        write (*,*)
        call xystat (xbox,ybox,nbox,
     +    rmsd,shap,corr,rf1,rf2,x6,x7,x8,x9)
        write (line,'(a,1p,e12.4,0p,a,f8.3,a,i6,a)')
     +    'Stats box-average X versus Y: RMSD ',rmsd,' ... CC ',
     +    corr,' ... ',nbox,' bins'
        call pretty (line)
        call xps_legend (line)
        line = ' '//line
        call prompt (line)
c
        call fitlin (nbox,xbox,ybox,slope,icept)
c
      end if
c
      write (*,*)
c
      return
      end
