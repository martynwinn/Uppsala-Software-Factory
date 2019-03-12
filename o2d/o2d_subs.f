c
      logical function xtx_implemented ()
c
code ...
c
      xtx_implemented = .false.
c
      return
      end
c
c ... the rest are dummy entries
c
      subroutine xtxini
      entry xtx_initialise
      return
      end
c
      subroutine xtxlab (iwin,icol,irow)
      entry xtx_labels_draw (iwin,icol,irow)
      return
      end
c
      subroutine xtxslb (lab)
      entry xtx_labels_on_off (lab)
      return
      end
c
      subroutine xtxdfl (horstr,verstr,ortstr)
      entry xtx_labels_define (horstr,verstr,ortstr)
      return
      end
c
      subroutine xtxdim (iwmax,icmax,irmax)
      entry xtx_dimensions (iwmax,icmax,irmax)
      return
      end
c
      subroutine xtxoff (xoff,yoff)
      entry xtx_set_2d_offset (xoff,yoff)
      return
      end
c
      subroutine xtxsca (xsc,ysc)
      entry xtx_set_2d_scales (xsc,ysc)
      return
      end
c
      subroutine xtxopn (u_name,w_nr,idim,ierr)
      entry xtx_open_window (u_name,w_nr,idim,ierr)
      return
      end
c
      subroutine xtxero (iwin)
      entry xtx_objects_erase (iwin)
      return
      end
c
      subroutine xtxalo (nrobj,o_name,lnew,ierr)
      entry xtx_object_alloc (nrobj,o_name,lnew,ierr)
      return
      end
c
      subroutine xtxedo (nrobj,ierr)
      entry xtx_object_edit (nrobj,ierr)
      return
      end
c
      subroutine xtxclo (nrobj,ierr)
      entry xtx_object_close (nrobj,ierr)
      return
      end
c
      subroutine xtxcuo
      entry xtx_objects_cleanup
      return
      end
c
      subroutine xtxiso (string,nobj)
      entry xtx_object_identify (string,nobj)
      return
      end
c
      subroutine xtxndo (string)
      entry xtx_named_object_delete (string)
      return
      end
c
      subroutine xtxdlo (nobj,ierr)
      entry xtx_object_delete (nobj,ierr)
      return
      end
c
      subroutine xtxslo (string)
      entry xtx_object_sleep (string)
      return
      end
c
      subroutine xtxaco (string)
      entry xtx_object_activate (string)
      return
      end
c
      subroutine xtxcpo (string,iwin,icol,irow)
      entry xtx_object_copy (string,iwin,icol,irow)
      return
      end
c
      subroutine xtxlio
      entry xtx_objects_list
      return
      end
c
      subroutine xtxera (icolor)
      entry xtx_erase_window (icolor)
      return
      end
c
      subroutine xtxcls
      entry xtx_close_window
      return
      end
c
      subroutine xtxcla
      entry xtx_close_all
      return
      end
c
      subroutine xtxrdw (w_nr)
      entry xtx_redraw_window (w_nr)
      return
      end
c
      subroutine xtxfup
      entry xtx_redraw_force
      return
      end
c
      subroutine xtxsrm(i)
      entry xtx_set_redraw(i)
      return
      end
c
      subroutine xtxgrm(i)
      entry xtx_get_redraw(i)
      return
      end
c
      subroutine xtxupd
      entry xtx_update_windows
      return
      end
c
      subroutine xtxslt (iw_nr,icol,irow)
      entry xtx_select_window (iw_nr,icol,irow)
      return
      end
c
      subroutine xtxmov (x,y)
      entry xtx_2d_move (x,y)
      return
      end
c
      subroutine xtxdrw(x,y)
      entry xtx_2d_draw (x,y)
      return
      end
c
      subroutine xtxasp(ix,iy)
      entry xtx_aspect_ratio (ix,iy)
      return
      end
c
      subroutine xtxsvp (w_nr,numcol,numrow,icolor)
      entry xtx_set_viewports (w_nr,numcol,numrow,icolor)
      return
      end
c
      subroutine xtxswi (icol,irow,xxlo,xxhi,yylo,yyhi)
      entry xtx_set_2d_world (icol,irow,xxlo,xxhi,yylo,yyhi)
      return
      end
c
      subroutine xtxdvc (iopt,icol,irow,ic2,ir2)
      entry xtx_vpt_connect (iopt,icol,irow,ic2,ir2)
      return
      end
c
      subroutine xtxdas (idir,itype,inter,ratio)
      entry xtx_axis_scale (idir,itype,inter,ratio)
      return
      end
c
      subroutine xtxgwi (icol,irow,xlo,xhi,ylo,yhi)
      entry xtx_get_2d_world (icol,irow,xlo,xhi,ylo,yhi)
      return
      end
c
      subroutine xtxcmv (iwindo,numcol,numrow,xphys,yphys)
      entry xtx_vpt_calc (iwindo,numcol,numrow,xphys,yphys)
      return
      end
c
      subroutine xtxktc
      entry xtx_keep_text_cursor
      return
      end
c
      subroutine xtxstc
      entry xtx_set_text_cursor
      return
      end
c
      subroutine xtxnmc (lokay)
      entry xtx_next_menu_cursor (lokay)
      return
      end
c
      subroutine xtxpzc (ddx,ddy)
      entry xtx_pan_zoom (ddx,ddy)
      return
      end
c
      subroutine xtxgvc (ix,iy,ic,ir,wx,wy)
      entry xtx_get_2d_coordinates (ix,iy,ic,ir,wx,wy)
      return
      end
c
      subroutine xtxlcc (iw,ic,ir,wx,wy)
      entry xtx_last_2d_cursor (iw,ic,ir,wx,wy)
      return
      end
c
      subroutine xtxstr (x,y,mytext)
      entry xtx_2d_string (x,y,mytext)
      return
      end
c
      subroutine xtxcol (icolor)
      entry xtx_set_colour (icolor)
      return
      end
c
      subroutine xtxgco (icolor)
      entry xtx_get_colour (icolor)
      return
      end
c
      subroutine xtxaxc(icolor)
      entry xtx_get_axis_colour (icolor)
      return
      end
c
      subroutine xtxerc(icolor)
      entry xtx_get_erase_colour (icolor)
      return
      end
c
      subroutine xtxcrs (jchar,x,y)
      entry xtx_cursor (jchar,x,y)
      entry xtxzcr (jchar,x,y,iwin,icol,irow)
      entry xtx_multi_cursor (jchar,x,y,iwin,icol,irow)
      return
      end
c
      subroutine xtxkcr
      entry xtx_keep_cursor
      return
      end
c
      subroutine xtxpnt (x,y)
      entry xtx_2d_pixel (x,y)
      return
      end
c
      subroutine xtxine (x,y,npt)
      entry xtx_2d_polyline (x,y,npt)
      return
      end
c
      subroutine xtxlin (x1,y1,x2,y2)
      entry xtx_2d_line (x1,y1,x2,y2)
      return
      end
c
      subroutine xtxbox (x1,x2,y1,y2)
      entry xtx_2d_box (x1,x2,y1,y2)
      return
      end
c
      subroutine xtxcir (x1,x2,y1)
      entry xtx_2d_circle (x1,x2,y1)
      return
      end
c
      subroutine xtxell (x1,x2,y1,y2)
      entry xtx_2d_ellipse (x1,x2,y1,y2)
      return
      end
c
      subroutine xtxgel (x1,x2,y1,y2,i1,i2,i3)
      entry xtx_2d_gen_ell (x1,x2,y1,y2,i1,i2,i3)
      return
      end
c
      subroutine xtxmrk (x1,x2,y1,y2)
      entry xtx_2d_mark (x1,x2,y1,y2)
      return
      end
c
      subroutine xtxcro (x1,x2,y1,y2)
      entry xtx_2d_cross (x1,x2,y1,y2)
      return
      end
c
      subroutine xtxpls (x1,x2,y1,y2)
      entry xtx_2d_plus (x1,x2,y1,y2)
      return
      end
c
      subroutine xtxgrd (x1,x2,xin,y1,y2,yin)
      entry xtx_2d_grid (x1,x2,xin,y1,y2,yin)
      return
      end
c
      subroutine xtxsym (isym,xx,x2,yy,y2)
      entry xtx_2d_symbol (isym,xx,x2,yy,y2)
      return
      end
c
      subroutine xtxdot (xx,x2,yy,y2,zmark)
      entry xtx_2d_dotline (xx,x2,yy,y2,zmark)
      return
      end
c
      subroutine xtxnpp (text,nr)
      entry xtx_popup_new (text,nr)
      return
      end
c
      subroutine xtxapp (nr,text,isubr)
      entry xtx_popup_add (nr,text,isubr)
      return
      end
c
      subroutine xtxdpp (nr)
      entry xtx_popup_delete (nr)
      return
      end
c
      subroutine xtxpia (icolor,iw,icol1,irow1,worl1x,worl1y,
     +                   icol2,irow2,worl2x,worl2y)
      entry xtx_pick_2d_area (icolor,iw,icol1,irow1,worl1x,worl1y,
     +                   icol2,irow2,worl2x,worl2y)
      return
      end
c
      subroutine xtxldi
      entry xtx_list_dims
      return
      end
c
      subroutine xt3svp (w_nr,icolor)
      entry xtx_3D_set_viewports (w_nr,icolor)
      return
      end
c
      subroutine xt3cub (iwindo,iobj)
      entry xtx_3D_draw_cube (iwindo)
      return
      end
c
      subroutine xt3per (ifovy,asp,znear,zfar)
      entry xtx_3D_perspective (ifovy,asp,znear,zfar)
      return
      end
c
      subroutine xt3win (xmin,xmax,ymin,ymax,zmin,zmax)
      entry xtx_3D_world (xmin,xmax,ymin,ymax,zmin,zmax)
      return
      end
c
      subroutine xt3mov (x,y,z)
      entry xtx_3D_move (x,y,z)
      return
      end
c
      subroutine xt3drw (x,y,z)
      entry xtx_3D_draw (x,y,z)
      return
      end
c
      subroutine xt3str (x,y,z,mytext)
      entry xtx_3d_string (x,y,z,mytext)
      return
      end
c
      subroutine xt3lat (vx,vy,vz,px,py,pz,itwist)
      entry xtx_3d_lookat (vx,vy,vz,px,py,pz,itwist)
      return
      end
