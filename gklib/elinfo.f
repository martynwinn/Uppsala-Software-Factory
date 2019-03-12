c
c ==========================================================================
c
      subroutine elinfo (abbr,name,nr,mass,radius,lprint)
c
c ... ELINFO - given the 2-character element symbol ABBR
c     (right justified), this routine returns:
c     NAME   = full name of the element
c     NR     = atomic nr (= nr of electrons if uncharged = nr of protons)
c     MASS   = atomic mass (Da ? amu ?)
c     RADIUS = covalent bond radius (A)
c
c ... error if returned NR <= 0
c
c ... if ABBR = '?' or '??' then the complete list is printed
c
c ... Gerard Kleywegt @ 950721/970707
c
      integer maxelm
      parameter (maxelm=103)
c
      real mass,radius
      integer nr
      character abbr*2,name*(*)
c
      real masses(maxelm),radii(maxelm)
c
      character symbol(maxelm)*2,names(maxelm)*15,myab*2
c
      integer i,length
c
      logical lprint
c
      data symbol /
     +  ' H',                              'HE',
     +  'LI','BE',' B',' C',' N',' O',' F','NE',
     +  'NA','MG','AL','SI',' P',' S','CL','AR',
     +  ' K','CA','SC','TI',' V','CR','MN','FE','CO','NI','CU',
     +       'ZN','GA','GE','AS','SE','BR','KR',
     +  'RB','SR',' Y','ZR','NB','MO','TC','RU','RH','PD','AG',
     +       'CD','IN','SN','SB','TE',' I','XE',
     +  'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB',
     +       'DY','HO','ER','TM','YB','LU','HF','TA',' W','RE',
     +       'OS','IR','PT','AU','HG','TL','PB','BI','PO','AT',
     +                                     'RN',
     +  'FR','RA','AC','TH','PA',' U','NP','PU','AM','CM','BK',
     +       'CF','ES','FM','MD','NO','LR'/
c
      data radii /
     +  0.23,  -1.0,0.68,0.35,0.83,0.68,
     +  0.68,0.68,0.64,  -1.0,0.97,1.10,
     +  1.35,1.20,1.05,1.02,0.99,
     +    -1.0,1.33,0.99,1.44,1.47,
     +  1.33,1.35,1.35,1.34,1.33,1.50,
     +  1.52,1.45,1.22,1.17,1.21,1.22,
     +  1.21,  -1.0,1.47,1.12,1.78,
     +  1.56,1.48,1.47,1.35,
     +  1.40,1.45,1.50,1.59,1.69,
     +  1.63,1.46,1.46,1.47,1.40,  -1.0,
     +  1.67,1.34,1.87,1.83,1.82,
     +  1.81,1.80,1.80,1.99,
     +  1.79,1.76,1.75,1.74,1.73,
     +  1.72,1.94,1.72,1.57,1.43,
     +  1.37,1.35,1.37,1.32,1.50,1.50,
     +  1.70,1.55,1.54,1.54,1.68,
     +    -1.0,  -1.0,  -1.0,1.90,1.88,
     +  1.79,1.61,1.58,1.55,
     +  1.53,1.51,  -1.0,  -1.0,
     +    -1.0,  -1.0,  -1.0,  -1.0,
     +    -1.0,  -1.0/
c
      data names /
     +  'Hydrogen','Helium','Lithium','Beryllium','Boron','Carbon',
     +  'Nitrogen','Oxygen','Fluorine','Neon','Sodium','Magnesium',
     +  'Aluminium','Silicon','Phosphorous','Sulfur','Chlorine',
     +  'Argon','Potassium','Calcium','Scandium','Titanium',
     +  'Vanadium','Chromium','Manganese','Iron','Cobalt','Nickel',
     +  'Copper','Zinc','Gallium','Germanium','Arsenic','Selenium',
     +  'Bromine','Krypton','Rubidium','Strontium','Yttrium',
     +  'Zirconium','Niobium','Molybdenum','Technetium',
     +  'Ruthenium','Rhodium','Palladium','Silver','Cadmium',
     +  'Indium','Tin','Antimony','Tellurium','Iodine','Xenon',
     +  'Caesium','Barium','Lanthanium','Cerium','Praseodymium',
     +  'Neodymium','Promethium','Samarium','Europium',
     +  'Gadolinium','Terbium','Dysprosium','Holmium','Erbium',
     +  'Thulium','Ytterbium','Lutetium','Hafnium','Tantalum',
     +  'Tungsten','Rhenium','Osmium','Iridium','Platinum','Gold',
     +  'Mercury','Thallium','Lead','Bismuth','Polonium',
     +  'Astatine','Radon','Francium','Radium','Actinium',
     +  'Thorium','Protactinium','Uranium','Neptunium',
     +  'Plutonium','Americium','Curium','Berkelium',
     +  'Californium','Einsteinium','Fermium','Mendelevium',
     +  'Nobelium','Lawrencium'/
c
      data masses /
     +  1.008,4.003,6.941,9.012,10.81,12.011,14.007,15.999,18.998,
     +  20.18,22.99,24.305,26.982,28.086,30.974,32.066,35.453,
     +  39.948,39.098,40.478,44.956,47.88,50.942,51.996,54.938,
     +  55.847,58.933,58.69,63.546,65.39,69.723,72.61,74.922,
     +  78.96,79.904,83.8,85.468,87.62,88.906,91.224,92.906,95.94,
     +  98,101.07,102.906,106.42,107.868,112.411,114.82,118.71,
     +  121.75,127.6,126.905,131.29,132.905,137.327,138.906,
     +  140.115,140.908,144.24,145,150.36,151.965,157.25,158.925,
     +  162.5,164.93,167.26,168.934,173.04,174.967,178.49,180.948,
     +  183.85,186.207,190.2,192.22,195.08,196.967,200.59,204.383,
     +  207.2,208.98,209,210,222,223,226.025,227.028,232.038,
     +  231.026,238.029,237.048,244,243,247,247,251,252,257,258,
     +  259,260/
c
code ...
c
      nr = -1
      mass = -1.0
      radius = -1.0
      name = '???'
c
      myab = abbr
      call upcase (myab)
      if (length(myab) .eq. 0) then
        if (lprint)
     +    call errcon ('ELINFO - no element symbol supplied')
        return
      else if (length(myab) .eq. 1) then
        myab = ' '//myab(1:1)
      end if
c
      if (myab .eq. ' ?' .or. myab .eq. '??') then
        do i=1,maxelm
          write (*,6000) i,symbol(i),names(i),masses(i)
        end do
        return
      end if
c
 6000 format (' At. # ',i3,' ',a2,' = ',a15,' Mass = ',f8.3,
     +  ' Radius = ',f8.3)
c
      do i=1,maxelm
        if (myab .eq. symbol(i)) then
          nr = i
          mass = masses (i)
          radius = radii (i)
          name = names (i)
          return
        end if
      end do
c
      if (lprint) then
        call errcon ('ELINFO - element symbol not recognised')
        call textut (' Symbol :',myab)
      end if
c
      return
      end
