C *************************************************************************
C 									  *
C  		 Copyright (C) 1984, Silicon Graphics, Inc.		  *
C 									  *
C   These coded instructions, statements, and computer programs  contain  *
C   unpublished  proprietary  information of Silicon Graphics, Inc., and  *
C   are protected by Federal copyright law.  They  may  not be disclosed  *
C   to  third  parties  or copied or duplicated in any form, in whole or  *
C   in part, without the prior written consent of Silicon Graphics, Inc.  *
C 									  *
C ************************************************************************


C  Graphics Libary constants 





C  various hardware/software limits 
       integer*4   ATTRIB
       parameter ( ATTRIB = 10 ) 	
       integer*4   VPSTAC
       parameter ( VPSTAC = 8 ) 	
       integer*4   MATRIX
       parameter ( MATRIX = 32 ) 	
       integer*4   NAMEST
       parameter ( NAMEST = 1025 ) 	

C  special pre-defined tags 
       integer*4   STARTT
       parameter ( STARTT = -2 ) 	
       integer*4   ENDTAG
       parameter ( ENDTAG = -3 ) 	

C  names for colors in color map loaded by greset 
       integer*4   BLACK
       parameter ( BLACK = 0 ) 	
       integer*4   RED
       parameter ( RED 	= 1 ) 	
       integer*4   GREEN
       parameter ( GREEN = 2 ) 	
       integer*4   YELLOW
       parameter ( YELLOW = 3 ) 	
       integer*4   BLUE
       parameter ( BLUE = 4 ) 	
       integer*4   MAGENT
       parameter ( MAGENT = 5 ) 	
       integer*4   CYAN
       parameter ( CYAN = 6 ) 	
       integer*4   WHITE
       parameter ( WHITE = 7 ) 	

C  popup colors 
       integer*4   PUPCLR
       parameter ( PUPCLR = 0 ) 	
       integer*4   PUPCOL
       parameter ( PUPCOL = 1 ) 	
       integer*4   PUPBLK
       parameter ( PUPBLK = 2 ) 	
       integer*4   PUPWHT
       parameter ( PUPWHT = 3 ) 	

C  defines for drawmode and mswapbuffers
       integer*4   NORMDR
       parameter ( NORMDR = 16 ) 	
       integer*4   PUPDRW
       parameter ( PUPDRW = 32 ) 	
       integer*4   OVRDRW
       parameter ( OVRDRW = 64 ) 	
       integer*4   UNDRDR
       parameter ( UNDRDR = 128 ) 	
       integer*4   CURSDR
       parameter ( CURSDR = 256 ) 	
       integer*4   DUALDR
       parameter ( DUALDR = 512 ) 	
       integer*4   GANGDR
       parameter ( GANGDR = 512 ) 	
       integer*4   HYPER0
       parameter ( HYPER0 = 65536 ) 	
       integer*4   HYPER1
       parameter ( HYPER1 = 131072 ) 	
       integer*4   HYPER2
       parameter ( HYPER2 = 262144 ) 	

C  defines for defpattern 
       integer*4   PAT16
       parameter ( PAT16 = 16 ) 	
       integer*4   PAT32
       parameter ( PAT32 = 32 ) 	
       integer*4   PAT64
       parameter ( PAT64 = 64 ) 	


C  defines for readsource 
       integer*4   SRCNUL
       parameter ( SRCNUL = -1 ) 	
       integer*4   SRCAUT
       parameter ( SRCAUT = 0 ) 	
       integer*4   SRCFRO
       parameter ( SRCFRO = 1 ) 	
       integer*4   SRCFL
       parameter ( SRCFL = 1 ) 	
       integer*4   SRCBAC
       parameter ( SRCBAC = 2 ) 	
       integer*4   SRCBL
       parameter ( SRCBL = 2 ) 	
       integer*4   SRCZBU
       parameter ( SRCZBU = 3 ) 	
       integer*4   SRCPUP
       parameter ( SRCPUP = 4 ) 	
       integer*4   SRCOVE
       parameter ( SRCOVE = 5 ) 	
       integer*4   SRCUND
       parameter ( SRCUND = 6 ) 	
       integer*4   SRCFRA
       parameter ( SRCFRA = 7 ) 	
       integer*4   SRCFR
       parameter ( SRCFR = 8 ) 	
       integer*4   SRCBR
       parameter ( SRCBR = 9 ) 	
       integer*4   SIL1
       parameter ( SIL1 = 1000 ) 	
       integer*4   SIL2
       parameter ( SIL2 = 1001 ) 	
       integer*4   SIL3
       parameter ( SIL3 = 1002 ) 	
       integer*4   SIL4
       parameter ( SIL4 = 1003 ) 	
       integer*4   SIL5
       parameter ( SIL5 = 1004 ) 	
       integer*4   SIL6
       parameter ( SIL6 = 1005 ) 	
       integer*4   SIL7
       parameter ( SIL7 = 1006 ) 	
       integer*4   SIL8
       parameter ( SIL8 = 1007 ) 	
       integer*4   SIL9
       parameter ( SIL9 = 1008 ) 	
       integer*4   SIL10
       parameter ( SIL10 = 1009 ) 	
       integer*4   SIL11
       parameter ( SIL11 = 1010 ) 	
       integer*4   SIL12
       parameter ( SIL12 = 1011 ) 	
       integer*4   SIL13
       parameter ( SIL13 = 1012 ) 	
       integer*4   SIL14
       parameter ( SIL14 = 1013 ) 	
       integer*4   SIL15
       parameter ( SIL15 = 1014 ) 	
       integer*4   SIL16
       parameter ( SIL16 = 1015 ) 	

C  minmax 
       integer*4   MMENA
       parameter ( MMENA = 0 ) 	
       integer*4   MMDIS
       parameter ( MMDIS = 1 ) 	
       integer*4   MMINI
       parameter ( MMINI = 2 ) 	

C  histogram 
       integer*4   HISTEN
       parameter ( HISTEN = 0 ) 	
       integer*4   HISTDI
       parameter ( HISTDI = 1 ) 	
       integer*4   HISTCL
       parameter ( HISTCL = 2 ) 	

C  defines for readcomponent 
       integer*4   RCABGR
       parameter ( RCABGR = 4 ) 	
       integer*4   RCALPH
       parameter ( RCALPH = 0 ) 	
       integer*4   RCBLUE
       parameter ( RCBLUE = 1 ) 	
       integer*4   RCGREE
       parameter ( RCGREE = 2 ) 	
       integer*4   RCRED
       parameter ( RCRED = 3 ) 	

C  defines for blendfunction 
       integer*4   BFZERO
       parameter ( BFZERO = 0 ) 	
       integer*4   BFONE
       parameter ( BFONE = 1 ) 	
       integer*4   BFDC
       parameter ( BFDC = 2 ) 	
       integer*4   BFSC
       parameter ( BFSC = 2 ) 	
       integer*4   BFMDC
       parameter ( BFMDC = 3 ) 	
       integer*4   BFMSC
       parameter ( BFMSC = 3 ) 	
       integer*4   BFSA
       parameter ( BFSA = 4 ) 	
       integer*4   BFMSA
       parameter ( BFMSA = 5 ) 	
       integer*4   BFDA
       parameter ( BFDA = 6 ) 	
       integer*4   BFMDA
       parameter ( BFMDA = 7 ) 	
       integer*4   BFMINS
       parameter ( BFMINS = 8 ) 	
       integer*4   BFCA
       parameter ( BFCA = 9 ) 	
       integer*4   BFMCA
       parameter ( BFMCA = 10 ) 	
       integer*4   BFCC
       parameter ( BFCC = 11 ) 	
       integer*4   BFMCC
       parameter ( BFMCC = 12 ) 	
       integer*4   BFMIN
       parameter ( BFMIN = 13 ) 	
       integer*4   BFMAX
       parameter ( BFMAX = 14 ) 	

C  defines for afunction 
       integer*4   AFNEVE
       parameter ( AFNEVE = 0 ) 	
       integer*4   AFLESS
       parameter ( AFLESS = 1 ) 	
       integer*4   AFEQUA
       parameter ( AFEQUA = 2 ) 	
       integer*4   AFLEQU
       parameter ( AFLEQU = 3 ) 	
       integer*4   AFGREA
       parameter ( AFGREA = 4 ) 	
       integer*4   AFNOTE
       parameter ( AFNOTE = 5 ) 	
       integer*4   AFGEQU
       parameter ( AFGEQU = 6 ) 	
       integer*4   AFALWA
       parameter ( AFALWA = 7 ) 	

C  defines for zfunction 
       integer*4   ZFNEVE
       parameter ( ZFNEVE = 0 ) 	
       integer*4   ZFLESS
       parameter ( ZFLESS = 1 ) 	
       integer*4   ZFEQUA
       parameter ( ZFEQUA = 2 ) 	
       integer*4   ZFLEQU
       parameter ( ZFLEQU = 3 ) 	
       integer*4   ZFGREA
       parameter ( ZFGREA = 4 ) 	
       integer*4   ZFNOTE
       parameter ( ZFNOTE = 5 ) 	
       integer*4   ZFGEQU
       parameter ( ZFGEQU = 6 ) 	
       integer*4   ZFALWA
       parameter ( ZFALWA = 7 ) 	

C  defines for zsource 
       integer*4   ZSRCDE
       parameter ( ZSRCDE = 0 ) 	
       integer*4   ZSRCCO
       parameter ( ZSRCCO = 1 ) 	

C  defines for pntsmooth 
       integer*4   SMPOFF
       parameter ( SMPOFF = 0 ) 	
       integer*4   SMPON
       parameter ( SMPON = 1 ) 	
       integer*4   SMPSMO
       parameter ( SMPSMO = 2 ) 	

C  defines for linesmooth 
       integer*4   SMLOFF
       parameter ( SMLOFF = 0 ) 	
       integer*4   SMLON
       parameter ( SMLON = 1 ) 	
       integer*4   SMLSMO
       parameter ( SMLSMO = 2 ) 	
       integer*4   SMLEND
       parameter ( SMLEND = 4 ) 	

C  defines for polysmooth, NOTE, these must be decimal to allow
C  macros in STAPUFT pe ucode (modes.m) to build
C 
       integer*4   PYSMOF
       parameter ( PYSMOF = 0 ) 	
       integer*4   PYSMON
       parameter ( PYSMON = 1 ) 	
       integer*4   PYSMSH
       parameter ( PYSMSH = 2 ) 	

C  defines for msalpha 
       integer*4   MSAMAS
       parameter ( MSAMAS = 0 ) 	
       integer*4   MSAONE
       parameter ( MSAONE = 1 ) 	
       integer*4   MSAALP
       parameter ( MSAALP = 2 ) 	

C  defines for mssample 
       integer*4   MSSPOI
       parameter ( MSSPOI = 0 ) 	
       integer*4   MSSARE
       parameter ( MSSARE = 1 ) 	
       integer*4   MSSCEN
       parameter ( MSSCEN = 2 ) 	

C  defines for mspattern 
       integer*4   MSPDEF
       parameter ( MSPDEF = 0 ) 	
       integer*4   MSP2P0
       parameter ( MSP2P0 = 1 ) 	
       integer*4   MSP2P1
       parameter ( MSP2P1 = 2 ) 	
       integer*4   MSP4P0
       parameter ( MSP4P0 = 3 ) 	
       integer*4   MSP4P1
       parameter ( MSP4P1 = 4 ) 	
       integer*4   MSP4P2
       parameter ( MSP4P2 = 5 ) 	
       integer*4   MSP4P3
       parameter ( MSP4P3 = 6 ) 	

C  dither modes	
       integer*4   DTOFF
       parameter ( DTOFF = 0 ) 	
       integer*4   DTON
       parameter ( DTON = 1 ) 	

C  defines for setpup 
       integer*4   PUPNON
       parameter ( PUPNON = 0 ) 	
       integer*4   PUPGRE
       parameter ( PUPGRE = 1 ) 	
       integer*4   PUPBOX
       parameter ( PUPBOX = 2 ) 	
       integer*4   PUPCHE
       parameter ( PUPCHE = 4 ) 	

C  defines for glcompat 
       integer*4   GLCOLD
       parameter ( GLCOLD = 0 ) 	
       integer*4   GLCZRA
       parameter ( GLCZRA = 1 ) 	
       integer*4   GLCMQU
       parameter ( GLCMQU = 2 ) 	
       integer*4   GLCSOF
       parameter ( GLCSOF = 3 ) 	
       integer*4   GLCMAN
       parameter ( GLCMAN = 4 ) 	
       integer*4   GLCSLO
       parameter ( GLCSLO = 5 ) 	
       integer*4   GLCINP
       parameter ( GLCINP = 6 ) 	
       integer*4   GLCNOB
       parameter ( GLCNOB = 7 ) 	
       integer*4   GLCSET
       parameter ( GLCSET = 8 ) 	
       integer*4   GLCGET
       parameter ( GLCGET = 9 ) 	
       integer*4   GLCVSY
       parameter ( GLCVSY = 10 ) 	
       integer*4   GLCFOR
       parameter ( GLCFOR = 11 ) 	
       integer*4   GLCTAB
       parameter ( GLCTAB = 12 ) 	
       integer*4   GLCNOD
       parameter ( GLCNOD = 13 ) 	
       integer*4   GLCNOR
       parameter ( GLCNOR = 14 ) 	
       integer*4   GLCFLU
       parameter ( GLCFLU = 15 ) 	
       integer*4   GLCILC
       parameter ( GLCILC = 1008 ) 	
       integer*4   GILRCM
       parameter ( GILRCM = 0 ) 	
       integer*4   GILRZM
       parameter ( GILRZM = 1 ) 	
       integer*4   GILCVL
       parameter ( GILCVL = 2 ) 	
       integer*4   GILLRP
       parameter ( GILLRP = 3 ) 	
       integer*4   GLCILH
       parameter ( GLCILH = 4 ) 	
       integer*4   GILGHG
       parameter ( GILGHG = 5 ) 	
       integer*4   GILMMX
       parameter ( GILMMX = 6 ) 	
       integer*4   GILGMM
       parameter ( GILGMM = 7 ) 	
       integer*4   GILPXM
       parameter ( GILPXM = 8 ) 	
       integer*4   GILPXT
       parameter ( GILPXT = 9 ) 	
       integer*4   GLCILI
       parameter ( GLCILI = 10 ) 	
       integer*4   GILDRW
       parameter ( GILDRW = 11 ) 	
       integer*4   GILSTL
       parameter ( GILSTL = 12 ) 	
       integer*4   GLCILL
       parameter ( GLCILL = 13 ) 	
       integer*4   GLCILR
       parameter ( GLCILR = 14 ) 	
       integer*4   GLCILS
       parameter ( GLCILS = 15 ) 	
       integer*4   GLCILM
       parameter ( GLCILM = 16 ) 	
       integer*4   GLCILF
       parameter ( GLCILF = 17 ) 	
       integer*4   GILTXL
       parameter ( GILTXL = 18 ) 	
       integer*4   GLCILB
       parameter ( GLCILB = 19 ) 	

C  value for GLC_MQUEUERATE which is compatible with 3.3 
       integer*4   GLCCOM
       parameter ( GLCCOM = 15 ) 	
C  value for GLC_MQUEUERATE which causes maximum compression 
       integer*4   GLCCMP
       parameter ( GLCCMP = 1 ) 	

C  defines for curstype 
       integer*4   C16X1
       parameter ( C16X1 = 0 ) 	
       integer*4   C16X2
       parameter ( C16X2 = 1 ) 	
       integer*4   C32X1
       parameter ( C32X1 = 2 ) 	
       integer*4   C32X2
       parameter ( C32X2 = 3 ) 	
       integer*4   CCROSS
       parameter ( CCROSS = 4 ) 	

C  defines for shademodel 
       integer*4   FLAT
       parameter ( FLAT = 0 ) 	
       integer*4   GOURAU
       parameter ( GOURAU = 1 ) 	

C  defines for logicop 
       integer*4   LOZERO
       parameter ( LOZERO = 0 ) 	
       integer*4   LOAND
       parameter ( LOAND = 1 ) 	
       integer*4   LOANDR
       parameter ( LOANDR = 2 ) 	
       integer*4   LOSRC
       parameter ( LOSRC = 3 ) 	
       integer*4   LOANDI
       parameter ( LOANDI = 4 ) 	
       integer*4   LODST
       parameter ( LODST = 5 ) 	
       integer*4   LOXOR
       parameter ( LOXOR = 6 ) 	
       integer*4   LOOR
       parameter ( LOOR = 7 ) 	
       integer*4   LONOR
       parameter ( LONOR = 8 ) 	
       integer*4   LOXNOR
       parameter ( LOXNOR = 9 ) 	
       integer*4   LONDST
       parameter ( LONDST = 10 ) 	
       integer*4   LOORR
       parameter ( LOORR = 11 ) 	
       integer*4   LONSRC
       parameter ( LONSRC = 12 ) 	
       integer*4   LOORI
       parameter ( LOORI = 13 ) 	
       integer*4   LONAND
       parameter ( LONAND = 14 ) 	
       integer*4   LOONE
       parameter ( LOONE = 15 ) 	

C  define for scrnselect 
       integer*4   INFOCU
       parameter ( INFOCU = -2 ) 	

C  defines for stencil 
       integer*4   STKEEP
       parameter ( STKEEP = 0 ) 	
       integer*4   STZERO
       parameter ( STZERO = 1 ) 	
       integer*4   STREPL
       parameter ( STREPL = 2 ) 	
       integer*4   STINCR
       parameter ( STINCR = 3 ) 	
       integer*4   STDECR
       parameter ( STDECR = 4 ) 	
       integer*4   STINVE
       parameter ( STINVE = 5 ) 	
       integer*4   SFNEVE
       parameter ( SFNEVE = 0 ) 	
       integer*4   SFLESS
       parameter ( SFLESS = 1 ) 	
       integer*4   SFEQUA
       parameter ( SFEQUA = 2 ) 	
       integer*4   SFLEQU
       parameter ( SFLEQU = 3 ) 	
       integer*4   SFGREA
       parameter ( SFGREA = 4 ) 	
       integer*4   SFNOTE
       parameter ( SFNOTE = 5 ) 	
       integer*4   SFGEQU
       parameter ( SFGEQU = 6 ) 	
       integer*4   SFALWA
       parameter ( SFALWA = 7 ) 	

C  defines for scrsubdivide 
       integer*4   SSOFF
       parameter ( SSOFF = 0 ) 	
       integer*4   SSDEPT
       parameter ( SSDEPT = 1 ) 	

C  defines for polymode 
       integer*4   PYMFIL
       parameter ( PYMFIL = 1 ) 	
       integer*4   PYMPOI
       parameter ( PYMPOI = 2 ) 	
       integer*4   PYMLIN
       parameter ( PYMLIN = 3 ) 	
       integer*4   PYMHOL
       parameter ( PYMHOL = 4 ) 	
       integer*4   PYMLNF
       parameter ( PYMLNF = 5 ) 	

C  defines for fogvertex 
       integer*4   FGOFF
       parameter ( FGOFF = 0 ) 	
       integer*4   FGON
       parameter ( FGON = 1 ) 	
       integer*4   FGDEFI
       parameter ( FGDEFI = 2 ) 	
       integer*4   FGVE
       parameter ( FGVE = 2 ) 	
       integer*4   FGVL
       parameter ( FGVL = 3 ) 	
       integer*4   FGPE
       parameter ( FGPE = 4 ) 	
       integer*4   FGPL
       parameter ( FGPL = 5 ) 	
       integer*4   FGVE2
       parameter ( FGVE2 = 6 ) 	
       integer*4   FGPE2
       parameter ( FGPE2 = 7 ) 	

C  defines for pixmode 
       integer*4   PMSHIF
       parameter ( PMSHIF = 0 ) 	
       integer*4   PMEXPA
       parameter ( PMEXPA = 1 ) 	
       integer*4   PMC0
       parameter ( PMC0 = 2 ) 	
       integer*4   PMC1
       parameter ( PMC1 = 3 ) 	
       integer*4   PMADD2
       parameter ( PMADD2 = 4 ) 	
       integer*4   PMSIZE
       parameter ( PMSIZE = 5 ) 	
       integer*4   PMOFFS
       parameter ( PMOFFS = 6 ) 	
       integer*4   PMSTRI
       parameter ( PMSTRI = 7 ) 	
       integer*4   PMTTOB
       parameter ( PMTTOB = 8 ) 	
       integer*4   PMRTOL
       parameter ( PMRTOL = 9 ) 	
       integer*4   PMZDAT
       parameter ( PMZDAT = 10 ) 	
       integer*4   PMWARP
       parameter ( PMWARP = 11 ) 	
       integer*4   PMRDX
       parameter ( PMRDX = 12 ) 	
       integer*4   PMRDY
       parameter ( PMRDY = 13 ) 	
       integer*4   PMCDX
       parameter ( PMCDX = 14 ) 	
       integer*4   PMCDY
       parameter ( PMCDY = 15 ) 	
       integer*4   PMXSTA
       parameter ( PMXSTA = 16 ) 	
       integer*4   PMYSTA
       parameter ( PMYSTA = 17 ) 	
       integer*4   PMVO1
       parameter ( PMVO1 = 1000 ) 	

C  defines for VENICE float path pixels/texels 
       integer*4   PMIFMT
       parameter ( PMIFMT = 11 ) 	
       integer*4   PMITYP
       parameter ( PMITYP = 12 ) 	
       integer*4   PMOFMT
       parameter ( PMOFMT = 13 ) 	
       integer*4   PMOTYP
       parameter ( PMOTYP = 14 ) 	
       integer*4   PMABGR
       parameter ( PMABGR = 0 ) 	
       integer*4   PMBGR
       parameter ( PMBGR = 1 ) 	
       integer*4   PMRGBA
       parameter ( PMRGBA = 2 ) 	
       integer*4   PMRGB
       parameter ( PMRGB = 3 ) 	
       integer*4   PMLUMI
       parameter ( PMLUMI = 4 ) 	
       integer*4   PMLUMA
       parameter ( PMLUMA = 5 ) 	
       integer*4   PMALPH
       parameter ( PMALPH = 6 ) 	
       integer*4   PMCOLO
       parameter ( PMCOLO = 7 ) 	
       integer*4   PMSTEN
       parameter ( PMSTEN = 8 ) 	
       integer*4   PMDEPT
       parameter ( PMDEPT = 9 ) 	
       integer*4   PMRED
       parameter ( PMRED = 10 ) 	
       integer*4   PMGREE
       parameter ( PMGREE = 11 ) 	
       integer*4   PMBLUE
       parameter ( PMBLUE = 12 ) 	


       integer*4   PMBITM
       parameter ( PMBITM = 0 ) 	
       integer*4   PMBYTE
       parameter ( PMBYTE = 1 ) 	
       integer*4   PMUNSI
       parameter ( PMUNSI = 2 ) 	
       integer*4   PMSH12
       parameter ( PMSH12 = 3 ) 	
       integer*4   PMUS12
       parameter ( PMUS12 = 4 ) 	
       integer*4   PMSHOR
       parameter ( PMSHOR = 5 ) 	
       integer*4   PMUSHT
       parameter ( PMUSHT = 6 ) 	
       integer*4   PMINT
       parameter ( PMINT = 7 ) 	
       integer*4   PMUINT
       parameter ( PMUINT = 8 ) 	
       integer*4   PMFLOA
       parameter ( PMFLOA = 9 ) 	

C  defines for convolve 
       integer*4   CVOFF
       parameter ( CVOFF = 0 ) 	
       integer*4   CVGENE
       parameter ( CVGENE = 1 ) 	
       integer*4   CVSEPA
       parameter ( CVSEPA = 2 ) 	
       integer*4   CVREDU
       parameter ( CVREDU = 1 ) 	

C  defines for pixelmap 
       integer*4   MAPITI
       parameter ( MAPITI = 0 ) 	
       integer*4   MAPITR
       parameter ( MAPITR = 1 ) 	
       integer*4   MAPITG
       parameter ( MAPITG = 2 ) 	
       integer*4   MAPITB
       parameter ( MAPITB = 3 ) 	
       integer*4   MAPITA
       parameter ( MAPITA = 4 ) 	
       integer*4   MAPRTR
       parameter ( MAPRTR = 5 ) 	
       integer*4   MAPGTG
       parameter ( MAPGTG = 6 ) 	
       integer*4   MAPBTB
       parameter ( MAPBTB = 7 ) 	
       integer*4   MAPATA
       parameter ( MAPATA = 8 ) 	

C  defines for pixeltransfer 
       integer*4   PMMPCO
       parameter ( PMMPCO = 0 ) 	
       integer*4   PTMAPD
       parameter ( PTMAPD = 0 ) 	
       integer*4   PMMPCL
       parameter ( PMMPCL = 1 ) 	
       integer*4   PTMAPE
       parameter ( PTMAPE = 2 ) 	

       integer*4   PTREDS
       parameter ( PTREDS = 1 ) 	
       integer*4   PTREDB
       parameter ( PTREDB = 2 ) 	
       integer*4   PTGRNS
       parameter ( PTGRNS = 3 ) 	
       integer*4   PTGRNB
       parameter ( PTGRNB = 4 ) 	
       integer*4   PTBLUS
       parameter ( PTBLUS = 5 ) 	
       integer*4   PTBLUB
       parameter ( PTBLUB = 6 ) 	
       integer*4   PTALPS
       parameter ( PTALPS = 7 ) 	
       integer*4   PTALPB
       parameter ( PTALPB = 8 ) 	

C  defines for nmode 
       integer*4   NAUTO
       parameter ( NAUTO = 0 ) 	
       integer*4   NNORMA
       parameter ( NNORMA = 1 ) 	

C  defines for acbuf 
       integer*4   ACCLR
       parameter ( ACCLR = 0 ) 	
       integer*4   ACACCU
       parameter ( ACACCU = 1 ) 	
       integer*4   ACCLRA
       parameter ( ACCLRA = 2 ) 	
       integer*4   ACRETU
       parameter ( ACRETU = 3 ) 	
       integer*4   ACMULT
       parameter ( ACMULT = 4 ) 	
       integer*4   ACADD
       parameter ( ACADD = 5 ) 	

C  defines for clipplane 
       integer*4   CPOFF
       parameter ( CPOFF = 0 ) 	
       integer*4   CPON
       parameter ( CPON = 1 ) 	
       integer*4   CPDEFI
       parameter ( CPDEFI = 2 ) 	

C  defines for scrbox 
       integer*4   SBRESE
       parameter ( SBRESE = 0 ) 	
       integer*4   SBTRAC
       parameter ( SBTRAC = 1 ) 	
       integer*4   SBHOLD
       parameter ( SBHOLD = 2 ) 	


C  defines for readdisplay 
       integer*4   RDFREZ
       parameter ( RDFREZ = 1 ) 	
       integer*4   RDALP1
       parameter ( RDALP1 = 2 ) 	
       integer*4   RDIGUN
       parameter ( RDIGUN = 4 ) 	
       integer*4   RDIGOV
       parameter ( RDIGOV = 8 ) 	
       integer*4   RDIGPU
       parameter ( RDIGPU = 16 ) 	
       integer*4   RDOFFS
       parameter ( RDOFFS = 32 ) 	

C  defines for getgconfig 
       integer*4   GCCMOD
       parameter ( GCCMOD = 0 ) 	
       integer*4   GCRED
       parameter ( GCRED = 1 ) 	
       integer*4   GCGREE
       parameter ( GCGREE = 2 ) 	
       integer*4   GCBLUE
       parameter ( GCBLUE = 3 ) 	
       integer*4   GCALPH
       parameter ( GCALPH = 4 ) 	
       integer*4   GCZBUF
       parameter ( GCZBUF = 5 ) 	
       integer*4   GCZMIN
       parameter ( GCZMIN = 6 ) 	
       integer*4   GCZMAX
       parameter ( GCZMAX = 7 ) 	
       integer*4   GCSTEN
       parameter ( GCSTEN = 8 ) 	
       integer*4   GCACBU
       parameter ( GCACBU = 9 ) 	
       integer*4   GCMS
       parameter ( GCMS = 10 ) 	
       integer*4   GCMSZ
       parameter ( GCMSZ = 11 ) 	
       integer*4   GCMZMI
       parameter ( GCMZMI = 12 ) 	
       integer*4   GCMZMA
       parameter ( GCMZMA = 13 ) 	
       integer*4   GCMSS
       parameter ( GCMSS = 14 ) 	
       integer*4   GCSTER
       parameter ( GCSTER = 15 ) 	
       integer*4   GCDOUB
       parameter ( GCDOUB = 16 ) 	


C  
C  START defines for getgdesc 
C 

       integer*4   GDXPMA
       parameter ( GDXPMA = 0 ) 	
       integer*4   GDYPMA
       parameter ( GDYPMA = 1 ) 	
       integer*4   GDXMMA
       parameter ( GDXMMA = 2 ) 	
       integer*4   GDYMMA
       parameter ( GDYMMA = 3 ) 	
       integer*4   GDZMIN
       parameter ( GDZMIN = 4 ) 	
       integer*4   GDZMAX
       parameter ( GDZMAX = 5 ) 	
       integer*4   GDBNSR
       parameter ( GDBNSR = 6 ) 	
       integer*4   GDBNSG
       parameter ( GDBNSG = 7 ) 	
       integer*4   GDBNSB
       parameter ( GDBNSB = 8 ) 	
       integer*4   GDBNDR
       parameter ( GDBNDR = 9 ) 	
       integer*4   GDBNDG
       parameter ( GDBNDG = 10 ) 	
       integer*4   GDBNDB
       parameter ( GDBNDB = 11 ) 	
       integer*4   GDBNSC
       parameter ( GDBNSC = 12 ) 	
       integer*4   GDBNDC
       parameter ( GDBNDC = 13 ) 	
       integer*4   GDBNSM
       parameter ( GDBNSM = 14 ) 	
       integer*4   GDBNDM
       parameter ( GDBNDM = 15 ) 	
       integer*4   GDBNZB
       parameter ( GDBNZB = 16 ) 	
       integer*4   GDBOSC
       parameter ( GDBOSC = 17 ) 	
       integer*4   GDBUSC
       parameter ( GDBUSC = 18 ) 	
       integer*4   GDBPSC
       parameter ( GDBPSC = 19 ) 	
       integer*4   GDBNSA
       parameter ( GDBNSA = 21 ) 	
       integer*4   GDBNDA
       parameter ( GDBNDA = 22 ) 	
       integer*4   GDBCUR
       parameter ( GDBCUR = 23 ) 	
       integer*4   GDOVER
       parameter ( GDOVER = 24 ) 	
       integer*4   GDBLEN
       parameter ( GDBLEN = 25 ) 	
       integer*4   GDCIFR
       parameter ( GDCIFR = 26 ) 	
       integer*4   GDXHCI
       parameter ( GDXHCI = 27 ) 	
       integer*4   GDDITH
       parameter ( GDDITH = 28 ) 	
       integer*4   GDLSCM
       parameter ( GDLSCM = 30 ) 	
       integer*4   GDLSRG
       parameter ( GDLSRG = 31 ) 	
       integer*4   GDLOGI
       parameter ( GDLOGI = 33 ) 	
       integer*4   GDNSCR
       parameter ( GDNSCR = 35 ) 	
       integer*4   GDNURB
       parameter ( GDNURB = 36 ) 	
       integer*4   GDNBLI
       parameter ( GDNBLI = 37 ) 	
       integer*4   GDNVPO
       parameter ( GDNVPO = 39 ) 	
       integer*4   GDPATS
       parameter ( GDPATS = 40 ) 	
       integer*4   GDPSCM
       parameter ( GDPSCM = 41 ) 	
       integer*4   GDPSRG
       parameter ( GDPSRG = 42 ) 	
       integer*4   GDPUPO
       parameter ( GDPUPO = 43 ) 	
       integer*4   GDRSRC
       parameter ( GDRSRC = 44 ) 	
       integer*4   GDRSZB
       parameter ( GDRSZB = 48 ) 	
       integer*4   GDSTER
       parameter ( GDSTER = 50 ) 	
       integer*4   GDSBLI
       parameter ( GDSBLI = 51 ) 	
       integer*4   GDSBPN
       parameter ( GDSBPN = 52 ) 	
       integer*4   GDSBPO
       parameter ( GDSBPO = 53 ) 	
       integer*4   GDTROD
       parameter ( GDTROD = 54 ) 	
       integer*4   GDWSYS
       parameter ( GDWSYS = 55 ) 	
       integer*4   GDZDGM
       parameter ( GDZDGM = 57 ) 	
       integer*4   GDZDPX
       parameter ( GDZDPX = 58 ) 	
       integer*4   GDSTYP
       parameter ( GDSTYP = 61 ) 	
       integer*4   GDTEXT
       parameter ( GDTEXT = 62 ) 	
       integer*4   GDNMMA
       parameter ( GDNMMA = 63 ) 	
       integer*4   GDFRAM
       parameter ( GDFRAM = 64 ) 	
       integer*4   GDTIME
       parameter ( GDTIME = 66 ) 	
       integer*4   GDDBBO
       parameter ( GDDBBO = 67 ) 	
       integer*4   GDAFUN
       parameter ( GDAFUN = 68 ) 	
       integer*4   GDALPH
       parameter ( GDALPH = 69 ) 	
       integer*4   GDBIAC
       parameter ( GDBIAC = 70 ) 	
       integer*4   GDBIAH
       parameter ( GDBIAH = 71 ) 	
       integer*4   GDBIST
       parameter ( GDBIST = 72 ) 	
       integer*4   GDCLIP
       parameter ( GDCLIP = 73 ) 	
       integer*4   GDFOGV
       parameter ( GDFOGV = 74 ) 	
       integer*4   GDLITW
       parameter ( GDLITW = 76 ) 	
       integer*4   GDPOMO
       parameter ( GDPOMO = 77 ) 	
       integer*4   GDPOSM
       parameter ( GDPOSM = 78 ) 	
       integer*4   GDSCRB
       parameter ( GDSCRB = 79 ) 	
       integer*4   GDTXTR
       parameter ( GDTXTR = 80 ) 	
       integer*4   GDFOGP
       parameter ( GDFOGP = 81 ) 	
       integer*4   GDTXPS
       parameter ( GDTXPS = 82 ) 	
       integer*4   GDMXPS
       parameter ( GDMXPS = 83 ) 	
       integer*4   GDMULT
       parameter ( GDMULT = 84 ) 	
       integer*4   GDTX3D
       parameter ( GDTX3D = 85 ) 	
       integer*4   GDTXLU
       parameter ( GDTXLU = 86 ) 	
       integer*4   GDTXSH
       parameter ( GDTXSH = 87 ) 	
       integer*4   GDTXDE
       parameter ( GDTXDE = 88 ) 	
       integer*4   GDSTIW
       parameter ( GDSTIW = 89 ) 	
       integer*4   GDBLCO
       parameter ( GDBLCO = 90 ) 	
       integer*4   GDLISP
       parameter ( GDLISP = 91 ) 	
       integer*4   GDLIA2
       parameter ( GDLIA2 = 92 ) 	
       integer*4   GDAFMO
       parameter ( GDAFMO = 93 ) 	
       integer*4   GDTXMS
       parameter ( GDTXMS = 94 ) 	
       integer*4   GDIGLO
       parameter ( GDIGLO = 95 ) 	
       integer*4   GDISNA
       parameter ( GDISNA = 96 ) 	


C  return value for inquiries when there is no limit 
       integer*4   GDNOLI
       parameter ( GDNOLI = -2 ) 	

C  return values for GD_WSYS 
       integer*4   GDWSNO
       parameter ( GDWSNO = 0 ) 	
       integer*4   GDWS4S
       parameter ( GDWS4S = 1 ) 	

C  return values for GD_SCRNTYPE 
       integer*4   GDSTWM
       parameter ( GDSTWM = 0 ) 	
       integer*4   GDSTNW
       parameter ( GDSTNW = 1 ) 	

C  
C  END defines for getgdesc 
C 


C  
C  START NURBS interface definitions 
C 

C  NURBS Rendering Properties 
       integer*4   NPIXEL
       parameter ( NPIXEL = 1 ) 	
       integer*4   NCULLI
       parameter ( NCULLI = 2 ) 	
       integer*4   NDISPL
       parameter ( NDISPL = 3 ) 	
       integer*4   NERROR
       parameter ( NERROR = 4 ) 	
       integer*4   NSUBDI
       parameter ( NSUBDI = 5 ) 	
       integer*4   NSSTEP
       parameter ( NSSTEP = 6 ) 	
       integer*4   NTSTEP
       parameter ( NTSTEP = 7 ) 	
       integer*4   NTILES
       parameter ( NTILES = 8 ) 	
       integer*4   NTMP1
       parameter ( NTMP1 = 9 ) 	
       integer*4   NTMP2
       parameter ( NTMP2 = 10 ) 	
       integer*4   NTMP3
       parameter ( NTMP3 = 11 ) 	
       integer*4   NTMP4
       parameter ( NTMP4 = 12 ) 	
       integer*4   NTMP5
       parameter ( NTMP5 = 13 ) 	
       integer*4   NTMP6
       parameter ( NTMP6 = 14 ) 	

       real*4      NFILL
       parameter ( NFILL = 1.0 ) 	
       real*4      NOUTPO
       parameter ( NOUTPO = 2.0 ) 	
       real*4      NOUTPA
       parameter ( NOUTPA = 5.0 ) 	
       real*4      NISOLI
       parameter ( NISOLI = 12.0 ) 	

C ---------------------------------------------------------------------------
C  FLAGS FOR NURBS SURFACES AND CURVES			
C  WARNING: Any changes to these flags should be checked against the 
C  decoding macros in nurbs.h.
C 
C  Bit: 876 5432 10 
C      |ttt|nnnn|rr|   :  rr - 2 bits = 1 if rational coordinate exists
C 		       : nnn - 4 bits for number of coordinates
C 		       : ttt - 3 bits for type of data (color, position, etc.)
C 	
C 
C  NURBS data type
C  N_T_ST	 	0	 parametric space data
C  N_T_XYZ		1	 model space data
C  N_T_TEX		2	 texture coordinate data
C  N_T_RGBA		3	 color data
C 
C  Number of coordinates per datum
C  N_COORD2	 	2	 2 coords
C  N_COORD3		3	 3 coords
C  N_COORD4		4	 4 coords
C  N_COORD5		5	 5 coords
C 
C  rational or non-rational data and position in memory 
C  N_NONRATIONAL	0	 non-rational data
C  N_RATIONAL		1	 rational data with rat coord after rest
C 
C 	
C ---------------------------------------------------------------------------
C 
       integer*4   NST
       parameter ( NST 	= 8 ) 	
       integer*4   NSTW
       parameter ( NSTW = 13 ) 	
       integer*4   NXYZ
       parameter ( NXYZ = 76 ) 	
       integer*4   NXYZW
       parameter ( NXYZW = 81 ) 	
       integer*4   NTEX
       parameter ( NTEX = 136 ) 	
       integer*4   NTEXW
       parameter ( NTEXW = 141 ) 	
       integer*4   NRGBA
       parameter ( NRGBA = 208 ) 	
       integer*4   NRGBAW
       parameter ( NRGBAW = 213 ) 	

C  New versions of above constants 

       integer*4   NP2D
       parameter ( NP2D = 8 ) 	
       integer*4   NP2DR
       parameter ( NP2DR = 13 ) 	
       integer*4   NV3D
       parameter ( NV3D = 76 ) 	
       integer*4   NV3DR
       parameter ( NV3DR = 81 ) 	
       integer*4   NT2D
       parameter ( NT2D = 136 ) 	
       integer*4   NT2DR
       parameter ( NT2DR = 141 ) 	
       integer*4   NC4D
       parameter ( NC4D = 208 ) 	
       integer*4   NC4DR
       parameter ( NC4DR = 213 ) 	

C  
C  END NURBS interface definitions 
C 


C  
C  START lighting model defines 
C 

       real*4      LMNULL
       parameter ( LMNULL = 0.0 ) 	

C  MATRIX modes	
       integer*4   MSINGL
       parameter ( MSINGL = 0 ) 	
       integer*4   MPROJE
       parameter ( MPROJE = 1 ) 	
       integer*4   MVIEWI
       parameter ( MVIEWI = 2 ) 	
       integer*4   MTEXTU
       parameter ( MTEXTU = 3 ) 	

C  LIGHT constants 
       integer*4   MAXLIG
       parameter ( MAXLIG = 8 ) 	
       integer*4   MAXRES
       parameter ( MAXRES = 4 ) 	

C  MATERIAL properties 
       integer*4   DEFMAT
       parameter ( DEFMAT = 0 ) 	
       integer*4   EMISSI
       parameter ( EMISSI = 1 ) 	
       integer*4   AMBIEN
       parameter ( AMBIEN = 2 ) 	
       integer*4   DIFFUS
       parameter ( DIFFUS = 3 ) 	
       integer*4   SPECUL
       parameter ( SPECUL = 4 ) 	
       integer*4   SHININ
       parameter ( SHININ = 5 ) 	
       integer*4   COLORI
       parameter ( COLORI = 6 ) 	
       integer*4   ALPHA
       parameter ( ALPHA = 7 ) 	

C  LIGHT properties 
       integer*4   DEFLIG
       parameter ( DEFLIG = 100 ) 	
       integer*4   LCOLOR
       parameter ( LCOLOR = 101 ) 	
       integer*4   POSITI
       parameter ( POSITI = 102 ) 	
       integer*4   SPOTDI
       parameter ( SPOTDI = 103 ) 	
       integer*4   SPOTLI
       parameter ( SPOTLI = 104 ) 	


C  LIGHTINGMODEL properties 
       integer*4   DEFLMO
       parameter ( DEFLMO = 200 ) 	
       integer*4   LOCALV
       parameter ( LOCALV = 201 ) 	
       integer*4   ATTENU
       parameter ( ATTENU = 202 ) 	
       integer*4   ATTEN2
       parameter ( ATTEN2 = 203 ) 	
       integer*4   TWOSID
       parameter ( TWOSID = 204 ) 	


C  TARGET constants 
       integer*4   MATERI
       parameter ( MATERI = 1000 ) 	
       integer*4   BACKMA
       parameter ( BACKMA = 1001 ) 	
       integer*4   LIGHT0
       parameter ( LIGHT0 = 1100 ) 	
       integer*4   LIGHT1
       parameter ( LIGHT1 = 1101 ) 	
       integer*4   LIGHT2
       parameter ( LIGHT2 = 1102 ) 	
       integer*4   LIGHT3
       parameter ( LIGHT3 = 1103 ) 	
       integer*4   LIGHT4
       parameter ( LIGHT4 = 1104 ) 	
       integer*4   LIGHT5
       parameter ( LIGHT5 = 1105 ) 	
       integer*4   LIGHT6
       parameter ( LIGHT6 = 1106 ) 	
       integer*4   LIGHT7
       parameter ( LIGHT7 = 1107 ) 	
       integer*4   LMODEL
       parameter ( LMODEL = 1200 ) 	

C  lmcolor modes 
       integer*4   LMCCOL
       parameter ( LMCCOL = 0 ) 	
       integer*4   LMCEMI
       parameter ( LMCEMI = 1 ) 	
       integer*4   LMCAMB
       parameter ( LMCAMB = 2 ) 	
       integer*4   LMCDIF
       parameter ( LMCDIF = 3 ) 	
       integer*4   LMCSPE
       parameter ( LMCSPE = 4 ) 	
       integer*4   LMCAD
       parameter ( LMCAD = 5 ) 	
       integer*4   LMCNUL
       parameter ( LMCNUL = 6 ) 	

C 
C  END lighting model defines 
C 


C  
C  START texturing defines 
C 

C  texdef param token values 
       integer*4   TXMINF
       parameter ( TXMINF = 256 ) 	
       integer*4   TXMAGF
       parameter ( TXMAGF = 512 ) 	
       integer*4   TXMAFC
       parameter ( TXMAFC = 2560 ) 	
       integer*4   TXMAFA
       parameter ( TXMAFA = 2816 ) 	
       integer*4   TXWRAP
       parameter ( TXWRAP = 768 ) 	
       integer*4   TXWRPS
       parameter ( TXWRPS = 784 ) 	
       integer*4   TXWRPT
       parameter ( TXWRPT = 800 ) 	
       integer*4   TXWRPR
       parameter ( TXWRPR = 816 ) 	
       integer*4   TXTILE
       parameter ( TXTILE = 1024 ) 	
       integer*4   TXBORD
       parameter ( TXBORD = 1280 ) 	
       integer*4   TXDETA
       parameter ( TXDETA = 3584 ) 	
       integer*4   TXFAST
       parameter ( TXFAST = 3840 ) 	
       integer*4   TXNOCO
       parameter ( TXNOCO = TXFAST ) 	
       integer*4   TXFRAM
       parameter ( TXFRAM = 4096 ) 	
       integer*4   TXSUBT
       parameter ( TXSUBT = 8192 ) 	
       integer*4   TXFORM
       parameter ( TXFORM = 8192 ) 	
       integer*4   TXFSTC
       parameter ( TXFSTC = 8208 ) 	
       integer*4   TXCPNT
       parameter ( TXCPNT = 3072 ) 	
       integer*4   TXCCLP
       parameter ( TXCCLP = 3328 ) 	
       integer*4   TXNULL
       parameter ( TXNULL = 0 ) 	

C  texture external formats 
       integer*4   TXEXTE
       parameter ( TXEXTE = 1792 ) 	
       integer*4   TXPK08
       parameter ( TXPK08 = 1808 ) 	
       integer*4   TXPK16
       parameter ( TXPK16 = 1824 ) 	
       integer*4   TXNOPA
       parameter ( TXNOPA = 1840 ) 	
       integer*4   TXPIXM
       parameter ( TXPIXM = 12288 ) 	
       integer*4   TXBCF
       parameter ( TXBCF = 16384 ) 	

       integer*4   TXINTE
       parameter ( TXINTE = 1536 ) 	
C  choices for TX_INTERNAL_FORMAT 
       integer*4   TXI12
       parameter ( TXI12 = 1552 ) 	
       integer*4   TXI12A
       parameter ( TXI12A = 1552 ) 	
       integer*4   TXI8
       parameter ( TXI8 = 1568 ) 	
       integer*4   TXIA8
       parameter ( TXIA8 = 1568 ) 	
       integer*4   TXRGB5
       parameter ( TXRGB5 = 1584 ) 	
       integer*4   TXRGBA
       parameter ( TXRGBA = 1600 ) 	
       integer*4   TXIA12
       parameter ( TXIA12 = 1616 ) 	
       integer*4   RGBA8
       parameter ( RGBA8 = 1632 ) 	
       integer*4   TXRGB8
       parameter ( TXRGB8 = 1632 ) 	
       integer*4   TXRGB1
       parameter ( TXRGB1 = 1664 ) 	
       integer*4   RGBA12
       parameter ( RGBA12 = 1648 ) 	
       integer*4   TXI16
       parameter ( TXI16 = 1680 ) 	


       integer*4   TXMIPM
       parameter ( TXMIPM = 2304 ) 	

C  texture filter choices 
       integer*4   TXPOIN
       parameter ( TXPOIN = 272 ) 	
       integer*4   TXBILI
       parameter ( TXBILI = 544 ) 	
       integer*4   TXMM
       parameter ( TXMM = 288 ) 	
       integer*4   TXMMP
       parameter ( TXMMP = 289 ) 	
       integer*4   TXMML
       parameter ( TXMML = 290 ) 	
       integer*4   TXMMBL
       parameter ( TXMMBL = 291 ) 	
       integer*4   TXMMTL
       parameter ( TXMMTL = 292 ) 	
       integer*4   TXBICU
       parameter ( TXBICU = 560 ) 	
       integer*4   TXSHAR
       parameter ( TXSHAR = 576 ) 	
       integer*4   TXMODU
       parameter ( TXMODU = 592 ) 	
       integer*4   TXADDD
       parameter ( TXADDD = 608 ) 	
       integer*4   TXTRIL
       parameter ( TXTRIL = 624 ) 	
       integer*4   TXMMQL
       parameter ( TXMMQL = 640 ) 	
       integer*4   TXBCGE
       parameter ( TXBCGE = 656 ) 	
       integer*4   TXBCLE
       parameter ( TXBCLE = 672 ) 	
       integer*4   TXBLGE
       parameter ( TXBLGE = 688 ) 	
       integer*4   TXBLLE
       parameter ( TXBLLE = 704 ) 	


C  texture wrapping access choices 
       integer*4   TXREPE
       parameter ( TXREPE = 769 ) 	
       integer*4   TXCLAM
       parameter ( TXCLAM = 770 ) 	
       integer*4   TXSELE
       parameter ( TXSELE = 771 ) 	

C  texture targets 
       integer*4   TXTEXT
       parameter ( TXTEXT = 0 ) 	
       integer*4   TXTEXD
       parameter ( TXTEXD = 1 ) 	
       integer*4   TXTEXI
       parameter ( TXTEXI = 2 ) 	

C  texture environment definitions 
       integer*4   TVMODU
       parameter ( TVMODU = 257 ) 	
       integer*4   TVBLEN
       parameter ( TVBLEN = 258 ) 	
       integer*4   TVDECA
       parameter ( TVDECA = 259 ) 	
       integer*4   TVCOLO
       parameter ( TVCOLO = 512 ) 	
       integer*4   TVSHAD
       parameter ( TVSHAD = 260 ) 	
       integer*4   TVALPH
       parameter ( TVALPH = 261 ) 	
       integer*4   TVCOMP
       parameter ( TVCOMP = 768 ) 	
       integer*4   TVIGR
       parameter ( TVIGR = 784 ) 	
       integer*4   TVIGG
       parameter ( TVIGG = 800 ) 	
       integer*4   TVIGB
       parameter ( TVIGB = 816 ) 	
       integer*4   TVIGA
       parameter ( TVIGA = 832 ) 	
       integer*4   TVIARG
       parameter ( TVIARG = 848 ) 	
       integer*4   TVIABA
       parameter ( TVIABA = 864 ) 	
       integer*4   TVIGI
       parameter ( TVIGI = 880 ) 	
       integer*4   TVNULL
       parameter ( TVNULL = 0 ) 	

C  texture lookup table definitions 
       integer*4   TLNULL
       parameter ( TLNULL = 0 ) 	
       integer*4   TLTLT0
       parameter ( TLTLT0 = 0 ) 	

C  texture environment targets 
       integer*4   TVENV0
       parameter ( TVENV0 = 0 ) 	

C  defines for texgen 
       integer*4   TXS
       parameter ( TXS 	= 0 ) 	
       integer*4   TXT
       parameter ( TXT 	= 1 ) 	
       integer*4   TXR
       parameter ( TXR 	= 2 ) 	
       integer*4   TXQ
       parameter ( TXQ 	= 3 ) 	
       integer*4   TGOFF
       parameter ( TGOFF = 0 ) 	
       integer*4   TGON
       parameter ( TGON = 1 ) 	
       integer*4   TGCONT
       parameter ( TGCONT = 2 ) 	
       integer*4   TGLINE
       parameter ( TGLINE = 3 ) 	
       integer*4   TGSPHE
       parameter ( TGSPHE = 4 ) 	
       integer*4   TGREFR
       parameter ( TGREFR = 5 ) 	

C 
C  END texturing defines 
C 

C  
C  START multi-byte character types for lcharstr and lstrwidth
C 

       integer*4   STRB
       parameter ( STRB = 0 ) 	
       integer*4   STR2B
       parameter ( STR2B = 1 ) 	
       integer*4   STR3B
       parameter ( STR3B = 16 ) 	
       integer*4   STR4B
       parameter ( STR4B = 17 ) 	
       integer*4   STR16
       parameter ( STR16 = 256 ) 	
       integer*4   STR32
       parameter ( STR32 = 257 ) 	

C  
C  END mbcharstr defines
C 

C  
C  START Distributed Graphics Library defines 
C 

       integer*4   DGLSIN
       parameter ( DGLSIN = 0 ) 	
       integer*4   DGLLOC
       parameter ( DGLLOC = 1 ) 	
       integer*4   DGLTSO
       parameter ( DGLTSO = 2 ) 	
       integer*4   DGL4DD
       parameter ( DGL4DD = 3 ) 	

C  
C  END Distributed Graphics Library defines 
C 


C  
C  START obsolete defines - included for backwards compatibility 
C 


       integer*4   PUPCUR
       parameter ( PUPCUR = PUPCOL ) 	

       integer*4   FATAL
       parameter ( FATAL = 1 ) 	
       integer*4   WARNIN
       parameter ( WARNIN = 2 ) 	
       integer*4   ASKCON
       parameter ( ASKCON = 3 ) 	
       integer*4   ASKRES
       parameter ( ASKRES = 4 ) 	

C  high-resolution monitor 
       integer*4   XMAXSC
       parameter ( XMAXSC = 1279 ) 	
       integer*4   YMAXSC
       parameter ( YMAXSC = 1023 ) 	

C  medium-resolution monitor 
       integer*4   XMAXME
       parameter ( XMAXME = 1023 ) 	
       integer*4   YMAXME
       parameter ( YMAXME = 767 ) 	

C  RS-170 
       integer*4   XMAX17
       parameter ( XMAX17 = 645 ) 	
       integer*4   YMAX17
       parameter ( YMAX17 = 484 ) 	

C  PAL 
       integer*4   XMAXPA
       parameter ( XMAXPA = 779 ) 	
       integer*4   YMAXPA
       parameter ( YMAXPA = 574 ) 	

C  
C  END obsolete defines
C 


C
C  START GL Function Declarations
C

       integer*4    blkqre
       integer*4    dglope
       integer*4    dopup
       integer*4    endfee
       integer*4    endpic
       integer*4    endsel
       integer*4    genobj
       integer*4    gentag
       integer*4    getbac
       integer*4    getbuf
       logical      getbut
       logical      getcmm
       integer*4    getcol
       logical      getdcm
       integer*4    getdes
       integer*4    getdis
       integer*4    getdra
       integer*4    getfon
       integer*4    getgco
       integer*4    getgde
       integer*4    gethei
       integer*4    gethit
C  gethit: not recommended
       logical      getlsb
C  getlsb: not recommended
       integer*4    getlsr
       integer*4    getlst
       integer*4    getlwi
       integer*4    getmap
       integer*4    getmmo
       integer*4    getmon
       logical      getmul
       integer*4    getope
       integer*4    getoth
C  getoth: obsolete
       integer*4    getpat
       integer*4    getpla
       logical      getres
C  getres: not recommended
       integer*4    getsha
C  getsha: obsolete
       integer*4    getsm
       integer*4    getval
       integer*4    getvid
       integer*4    getwri
       integer*4    getwsc
       logical      getzbu
       integer*4    gversi
       integer*4    ilbuff
       logical      ismex
C  ismex: obsolete
       logical      isobj
       logical      isqueu
       logical      istag
       integer*4    istexl
       integer*4    lrectr
       integer*4    lstwid
       integer*4    newpup
       integer*4    qgetfd
       integer*4    qread
       integer*4    qtest
       integer*4    readdi
       integer*4    readpi
       integer*4    readRG
       integer*4    rectre
       integer*4    scrnat
       integer*4    scrnse
       integer*4    strwid
       integer*4    swinop
       integer*4    winatt
C  winatt: obsolete
       integer*4    windep
       integer*4    winget
       integer*4    winope

C
C  END GL Function Declarations
C

