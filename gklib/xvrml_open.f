c
c ... XVRML_SUBS.F - subroutines for handling VRML files
c
c ... Gerard J Kleywegt - July, 1997
c
c ... Update 20051205
c
c      xvrml_init ()
c      xvrml_open (iunit,bg1,bg2,bg3)
c      xvrml_close ()
c
c      xvrml_encode_rgb (r,g,b,irgb)
c      xvrml_decode_rgb (r,g,b,irgb)
c      xvrml_colour (rgb1,rgb2,rgb3)
c      xvrml_col_index (i,j)
c      xvrml_col_rgb (i,r,g,b)
c      xvrml_col_list ()
c      xvrml_col_name (col,r,g,b,icol)
c      xvrml_rgb_name (colnam,r,g,b)
c
c      xvrml_cpk (nat,xyz,rad)
c      xvrml_cpk_col (nat,xyz,rad,icol)
c
c      xvrml_wire (nat,xyz,maxbnd,nbonds,bndptr)
c      xvrml_wire_col (nat,xyz,maxbnd,nbonds,bndptr,icol)
c
c      xvrml_stick (nat,xyz,maxbnd,nbonds,bndptr,strad)
c      xvrml_stick_col (nat,xyz,maxbnd,nbonds,bndptr,strad,icol)
c
c      xvrml_polyline (nat,xyz)
c      xvrml_pointset (nat,xyz)
c
c      xvrml_wire_surf (nx,ny,z,x1,x2,y1,y2)
c      xvrml_face_surf (nx,ny,z,x1,x2,y1,y2)
c
c      xvrml_scale (sx,sy,sz)
c      xvrml_plus (nat,xyz,icol,size)
c      xvrml_get_cyl (x1,x2,a,b,c,d,e,xc)
c
c
      subroutine xvrml_open (iunit,bg1,bg2,bg3)
c
      include 'xvrml.incl'
c
      real bg1,bg2,bg3
c
      integer iunit,leng1
c
      character line*128
c
code ...
c
      ixvrml = 99
      if (iunit .gt. 0) ixvrml = iunit
      lxvopn = .true.
c
      write (ixvrml,'(a)',err=9999)
     +  '#VRML V1.0 ascii'
c
      call stamp (line)
      write (ixvrml,'(a1,1x,a)',err=9999)
     +  '#',line(1:leng1(line))
c
      write (ixvrml,'(3a)',err=9999) 
     +  'DEF SceneInfo Info { string "',verstr,'" }'
c
      write (ixvrml,'(a/3a)',err=9999) 
     +  'DEF Title Info { ',
     +  'string "',line(1:leng1(line)),'" }'
c
      write (ixvrml,6000,err=9999) bg1,bg2,bg3
 6000 format ('DEF BackgroundColor Info { string "',
     +  3f8.3,'" }')
c
      write (ixvrml,6010,err=9999)
 6010 format ('FontStyle { size 2 }')
c
      call prompt (' Opened VRML file')
c
      return
c
 9999 continue
      call errcon ('While writing VRML file')
      close (ixvrml)
      lxvopn = .false.
c
      return
      end
