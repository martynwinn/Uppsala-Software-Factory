c
c ==========================================================================
c
c ... Latest update @ 2007-08-29
c
      subroutine gkrefs (prog,nr)
c
c ... subroutine to print literature references
c
      implicit none
c
      integer nr,ll,mm,length,index,lp1,lp2,li1,li2
c
      logical lrave,lvoid,ldeja,lspaz,lsbin
c
      character prog*(*),mypr*80
      character pubpre*100,pubpos*100,iucpre*100,iucpos*100
c
code ...
c
      nr = 0
      ll = length(prog)
      if (ll .lt. 1) return
c
c ... PubMed URL format
c
      pubpre = '[http://www.ncbi.nlm.nih.gov/entrez/' //
     +  'query.fcgi?db=PubMed&cmd=Retrieve&list_uids='
      pubpos = '&dopt=Citation]'
c
c ... IUCr URL format
c
      iucpre = '[http://scripts.iucr.org/cgi-bin/paper?'
      iucpos = ']'
c
      lp1 = length(pubpre)
      lp2 = length(pubpos)
      li1 = length(iucpre)
      li2 = length(iucpos)
c
      mypr = '|'//prog(1:ll)//'|'
      mm = ll+2
c
      lvoid = (index('|VOIDOO|MAPROP|FLOOD|',mypr(1:mm)).gt.0)
      ldeja = (index('|DEJAVU|PRO1|PRO2|',mypr(1:mm)).gt.0 .or.
     +         index('|POST|DEJANA|LSQMAN|',mypr(1:mm)).gt.0 .or.
     +         index('|GETSSE|',mypr(1:mm)).gt.0)
      lrave = (index('|IMP|AVE|COMAP|COMDEM|',mypr(1:mm)).gt.0 .or.
     +         index('|MAVE|MAPFIX|NCS6D|',mypr(1:mm)).gt.0 .or.
     +         index('|MAPMAN|MAMA|DATAMAN|',mypr(1:mm)).gt.0 .or.
     +         index('|CRAVE|ESSENS|MAPPAGE|',mypr(1:mm)).gt.0 .or.
     +         index('|PROF|SOLEX|SPANCSI|',mypr(1:mm)).gt.0 .or.
     +         index('|SITE2RT|COMA|MASKIT|',mypr(1:mm)).gt.0 .or.
     +         index('|SSENCS|FINDNCS|',mypr(1:mm)).gt.0)
      lspaz = (index('|SPASM|MKSPAZ|RIGOR|',mypr(1:mm)).gt.0 .or.
     +         index('|AUTOMOTIF|NHANCE|',mypr(1:mm)).gt.0 .or.
     +         index('|SAVANT|MAKRIG|DEJANA|',mypr(1:mm)).gt.0 .or.
     +         index('|SAVANA|SPANA|MKSPAN|',mypr(1:mm)).gt.0)
      lsbin = (index('|STRUPRO|ZPROF|PRF2MSEQ|',mypr(1:mm)).gt.0 .or.
     +         index('|MSEQPRO|STRUPAT|MSEQ2ALSQ|',mypr(1:mm)).gt.0)
c
c ... 1986 ..............................................................
c
      if (index('|MAPMAN|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'T.A. Jones & S. Thirup (1986). Using known substructures in',
     +'protein model building and crystallography. EMBO J 5, 819-822.',
     +pubpre(1:lp1)//'3709525'//pubpos(1:lp2),
     +' ',' ',' ')
      end if
c
c ... 1992 ..............................................................
c
      if (lrave) then
        call priref (nr,
     +'T.A. Jones (1992). A, yaap, asap, @#*?  A set of averaging',
     +'programs. In "Molecular Replacement", edited by E.J. Dodson,',
     +'S. Gover and W. Wolf. SERC Daresbury Laboratory, Warrington,',
     +'pp. 91-105.',
     +' ',' ')
      end if
c
c ... 1993 ..............................................................
c
      if (index('|MAMA|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1993).  Masks Made Easy.',
     +'CCP4/ESF-EACBM Newsletter on Protein Crystallography 28,',
     +'May 1993, pp. 56-59.',
     +'[http://xray.bmc.uu.se/usf/factory_1.html]',
     +' ',' ')
      end if
c
      if (lvoid) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1993). Biomacromolecular',
     +'Speleology. CCP4/ESF-EACBM Newsletter on Protein',
     +'Crystallography 29, November 1993, pp. 26-28.',
     +'[http://xray.bmc.uu.se/usf/factory_2.html]',
     +' ',' ')
      end if
c
c ... 1994 ..............................................................
c
      if (lvoid) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1994).  Detection, delineation,',
     +'measurement and display of cavities in macromolecular',
     +'structures. Acta Cryst D50, 178-185.',
     +iucpre(1:li1)//'gr0263'//iucpos(1:li2),
     +' ',' ')
      end if
c
      if (lrave .or. ldeja) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1994).  Halloween ... Masks and',
     +'Bones. In "From First Map to Final Model", edited by',
     +'S. Bailey, R. Hubbard and D. Waller.  SERC Daresbury',
     +'Laboratory, Warrington, pp. 59-66.',
     +'[http://xray.bmc.uu.se/gerard/papers/halloween.html]',
     +' ')
      end if
c
      if (index('|OOPS|OOPS2|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1994). OOPS-a-daisy.',
     +'CCP4/ESF-EACBM Newsletter on Protein Crystallography 30,',
     +'June 1994, pp. 20-24.',
     +'[http://xray.bmc.uu.se/usf/factory_3.html]',
     +' ',' ')
      end if
c
      if (index('|LSQMAN|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1994). A super position.',
     +'CCP4/ESF-EACBM Newsletter on Protein Crystallography 31,',
     +'November 1994, pp. 9-14.',
     +'[http://xray.bmc.uu.se/usf/factory_4.html]',
     +' ',' ')
      end if
c
c ... 1995 ..............................................................
c
      if (index('|LSQMAN|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1995). Where freedom is given,',
     +'liberties are taken. Structure 3, 535-540.',
     +pubpre(1:lp1)//'8590014'//pubpos(1:lp2),
     +' ',' ',' ')
      end if
c
      if (index('|MOLEMAN2|XPLO2D|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt (1995). Dictionaries for Heteros.',
     +'CCP4/ESF-EACBM Newsletter on Protein Crystallography 31,',
     +'June 1995, pp. 45-50.',
     +'[http://xray.bmc.uu.se/usf/factory_5.html]',
     +' ',' ')
      end if
c
      if (index('|ESSENS|SPASM|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'M. Harel, G.J. Kleywegt, R.B.G. Ravelli, I. Silman & J.L.',
     +'Sussman (1995). Crystal structure of an acetylcholinesterase-',
     +'fasciculin complex: interaction of a three-fingered toxin',
     +'from snake venom with its target. Structure 3, 1355-1366.',
     +pubpre(1:lp1)//'8747462'//pubpos(1:lp2),
     +' ')
      end if
c
c ... 1996 ..............................................................
c
      if (index('|MAPMAN|DATAMAN|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1996). xdlMAPMAN and xdlDATAMAN -',
     +'programs for reformatting, analysis and manipulation of',
     +'biomacromolecular electron-density maps and reflection data',
     +'sets. Acta Cryst D52, 826-828.',
     +iucpre(1:li1)//'gr0468'//iucpos(1:li2),
     +' ')
      end if
c
      if (index('|OOPS|OOPS2|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1996). Efficient rebuilding of',
     +'protein structures. Acta Cryst D52, 829-832.',
     +iucpre(1:li1)//'gr0469'//iucpos(1:li2),
     +' ',' ',' ')
      end if
c
      if (index('|LSQMAN|QDB|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt (1996). Use of non-crystallographic symmetry',
     +'in protein structure refinement. Acta Cryst D52, 842-857.',
     +iucpre(1:li1)//'gr0471'//iucpos(1:li2),
     +' ',' ',' ')
      end if
c
      if (index('|MOLEMAN2|XPLO2D|LSQMAN|',mypr(1:mm)).gt.0 .or.
     +    index('|SEAMAN|MAPMAN|PACMAN|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt (1996). Making the most of your search model.',
     +'CCP4/ESF-EACBM Newsletter on Protein Crystallography 32,',
     +'June 1996, pp. 32-36.',
     +'[http://xray.bmc.uu.se/usf/factory_6.html]',
     +' ',' ')
      end if
c
      if (index('|DATAMAN|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt & A.T. Brunger (1996). Checking your',
     +'imagination: applications of the free R value.',
     +'Structure 4, 897-904.',
     +pubpre(1:lp1)//'8805582'//pubpos(1:lp2),
     +' ',' ')
      end if
c
      if (index('|MOLEMAN2|OOPS|OOPS2|LSQMAN|',mypr(1:mm))
     +    .gt. 0) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1996). Phi/Psi-chology:',
     +'Ramachandran revisited. Structure 4, 1395-1400.',
     +pubpre(1:lp1)//'8994966'//pubpos(1:lp2),
     +' ',' ',' ')
      end if
c
c ... 1997 ..............................................................
c
      if (index('|ESSENS|SOLEX|DEJAVU|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1997). Taking the fun out of map',
     +'interpretation. CCP4/ESF-EACBM Newsletter on Protein',
     +'Crystallography 33, January 1997, pp. 19-21.',
     +'[http://xray.bmc.uu.se/usf/factory_7.html]',
     +' ',' ')
      end if
c
      if (index('|ESSENS|SOLEX|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1997).  Template convolution to',
     +'enhance or detect structural features in macromolecular',
     +'electron-density maps. Acta Cryst D53, 179-185.',
     +iucpre(1:li1)//'gr0647'//iucpos(1:li2),
     +' ',' ')
      end if
c
      if (index('|SOD|O2D|ODBMAN|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt (1997). Les amis d"O. CCP4/ESF-EACBM',
     +'Newsletter on Protein Crystallography 34, September 1997,',
     +'pp. 5-8.',
     +'[http://xray.bmc.uu.se/usf/factory_8.html]',
     +' ',' ')
      end if
c
      if (index('|MOLEMAN2|XPLO2D|OOPS|OOPS2|',mypr(1:mm))
     +    .gt. 0) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1997). Model-building and',
     +'refinement practice. Methods in Enzymology 277, 208-230.',
     +'[http://xray.bmc.uu.se/gerard/gmrp/gmrp.html]',
     +' ',' ',' ')
      end if
c
      if (ldeja .or. (index('|LSQMAN|',mypr(1:mm)).gt.0)) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1997). Detecting folding motifs',
     +'and similarities in protein structures. Methods in',
     +'Enzymology 277, 525-545.',
     +' ',' ',' ')
      end if
c
      if (index('|MOLEMAN2|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt (1997). Validation of protein models from CA',
     +'coordinates alone. J Mol Biol 273, 371-376.',
     +pubpre(1:lp1)//'9344745'//pubpos(1:lp2),
     +' ',' ',' ')
      end if
c
      if (lrave) then
        call priref (nr,
     +'G.J. Kleywegt & R.J. Read (1997). Not your average density.',
     +'Structure 5, 1557-1569.',
     +pubpre(1:lp1)//'9438862'//pubpos(1:lp2),
     +' ',' ',' ')
      end if
c
c ... 1998 ..............................................................
c
      if (lspaz) then
        call priref (nr,
     +'G.J. Kleywegt (1998). Deja-vu all over again. CCP4 Newsletter',
     +'on Protein Crystallography 35, July 1998, pp. 10-12.',
     +'[http://xray.bmc.uu.se/usf/factory_9.html]',
     +' ',' ',' ')
      end if
c
      if (lspaz .or. lsbin .or.
     +    (index('|HETZE|',mypr(1:mm)).gt.0)) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1998). Databases in protein',
     +'crystallography. Acta Cryst D54, 1119-1131.',
     +'[http://xray.bmc.uu.se/gerard/papers/databases.html]',
     +pubpre(1:lp1)//'10089488'//pubpos(1:lp2),
     +iucpre(1:li1)//'ba0001'//iucpos(1:li2),
     +' ')
      end if
c
c ... 1999 ..............................................................
c
      if (lspaz) then
        call priref (nr,
     +'G.J. Kleywegt (1999). Recognition of spatial motifs in',
     +'protein structures. J Mol Biol 285, 1887-1897.',
     +pubpre(1:lp1)//'9917419'//pubpos(1:lp2),
     +  ' ',' ',' ')
      end if
c
      if (index('|MAMA|COMA|MASKIT|',mypr(1:mm)).gt.0 .or.
     +    index('|MAPMAN|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (1999). Software for handling',
     +'macromolecular envelopes. Acta Cryst D55, 941-944.',
     +pubpre(1:lp1)//'10089342'//pubpos(1:lp2),
     +iucpre(1:li1)//'se0256'//iucpos(1:li2),
     +' ',' ')
      end if
c
      if (index('|LSQMAN|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'T.A. Jones & G.J. Kleywegt (1999). CASP3 comparative',
     +'modelling evaluation.',
     +'Proteins: Struct. Funct. Genet. Suppl. 3, 30-46.',
     +pubpre(1:lp1)//'10526350'//pubpos(1:lp2),
     +'[http://xray.bmc.uu.se/casp3]',
     +' ')
      end if
c
      if (index('|MAMA|COMA|MASKIT|',mypr(1:mm)).gt.0 .or.
     +    index('|MAPMAN|LSQMAN|MOLEMAN2|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'G.J. Kleywegt (1999). Experimental assessment of',
     +'differences between related protein crystal structures.',
     +'Acta Cryst. D55, 1878-1857.',
     +pubpre(1:lp1)//'10531486'//pubpos(1:lp2),
     +iucpre(1:li1)//'se0283'//iucpos(1:li2),
     +' ')
      end if
c
c ... 2000
c
      if (index('|MOLEMAN2|HETZE|OOPS|OOPS2|',mypr(1:mm))
     +    .gt. 0) then
        call priref (nr,
     +'G.J. Kleywegt (2000). Validation of protein crystal',
     +'structures. Acta Cryst. D56, 249-265 (Topical Review).',
     +pubpre(1:lp1)//'10713511'//pubpos(1:lp2),
     +iucpre(1:li1)//'gr0949'//iucpos(1:li2),
     +' ',' ')
      end if
c
      if (index('|SEAMAN|LSQMAN|',mypr(1:mm))
     +    .gt. 0) then
        call priref (nr,
     +'Y.W. Chen, E.J. Dodson & G.J. Kleywegt (2000). Does NMR',
     +'mean "Not for Molecular Replacement" ? Using NMR-based',
     +'search models to solve protein crystal structures.',
     +'Structure 8, R213-R220.',
     +pubpre(1:lp1)//'11080645'//pubpos(1:lp2),
     +' ')
      end if
c
c ... 2001 ..............................................................
c
      if (lrave) then
        call priref (nr,
     +'R.J. Read & G.J. Kleywegt (2001). Density modification:',
     +'theory and practice. In: "Methods in Macromolecular',
     +'Crystallography" (D Turk & L Johnson, Eds.), IOS Press,',
     +'Amsterdam, pp. 123-135.',
     +' ',' ')
      end if
c
      if (index('|MOLEMAN2|HETZE|OOPS|OOPS2|',mypr(1:mm))
     +    .gt. 0) then
        call priref (nr,
     +'G.J. Kleywegt (2001). Validation of protein crystal structures.',
     +'In: "International Tables for Crystallography, Volume F.',
     +'Crystallography of Biological Macromolecules" (Rossmann, M.G.',
     +'& Arnold, E., Editors). Chapter 21.1, pp. 497-506, 526-528.',
     +'Dordrecht: Kluwer Academic Publishers, The Netherlands.',
     +' ')
      end if
c
c ... 2002 ..............................................................
c
      if (lspaz .or. ldeja) then
        call priref (nr,
     +'D. Madsen & G.J. Kleywegt (2002). Interactive motif and',
     +'fold recognition in protein structures. J. Appl. Cryst.',
     +'35, 137-139.',
     +iucpre(1:li1)//'wt0007'//iucpos(1:li2),
     +' ',' ')
      end if
c
      if (index('|MOLEMAN2|OOPS|OOPS2|',mypr(1:mm))
     +    .gt. 0) then
        call priref (nr,
     +'G.J. Kleywegt & T.A. Jones (2002). Homo Crystallographicus -',
     +'Quo Vadis ? Structure 10, 465-472.',
     +pubpre(1:lp1)//'11937051'//pubpos(1:lp2),
     +' ',' ',' ')
      end if
c
c ... 2003 ..............................................................
c
      if (index('|XPLO2D|',mypr(1:mm))
     +    .gt. 0) then
        call priref (nr,
     +'G.J. Kleywegt, K. Henrick, E.J. Dodson & D.M.F. van Aalten.',
     +'Pound-wise but penny-foolish: How well do micromolecules',
     +'fare in macromolecular refinement ?',
     +'Structure, 11, 1051-1059 (2003).',
     +pubpre(1:lp1)//'12962624'//pubpos(1:lp2),
     +' ')
      end if
c
c ... 2004 ..............................................................
c
      if (ldeja) then
        call priref (nr,
     +'M. Novotny, D. Madsen & G.J. Kleywegt (2004). An evaluation',
     +'of protein-fold-comparison servers. Proteins, 54, 260-270',
     +'(2004).',
     +pubpre(1:lp1)//'14696188'//pubpos(1:lp2),
     + ' ',' ')
      end if
c
      if (index('|MOLEMAN2|MAPMAN|DATAMAN|',mypr(1:mm)) .gt. 0 .or.
     +    index('|STAT2O|FILPDB|ELAL|',mypr(1:mm)) .gt. 0) then
        call priref (nr,
     +'G.J. Kleywegt, M.R. Harris, J. Zou, T.C. Taylor, ',
     +'A. Wahlby & T.A. Jones (2004). The Uppsala Electron-',
     +'Density Server. Acta Cryst. D60, 2240-2249.',
     +pubpre(1:lp1)//'15572777'//pubpos(1:lp2),
     +iucpre(1:li1)//'ba5061'//iucpos(1:li2),
     + ' ')
      end if
c
      if (ldeja .or. (index('|LSQMAN|',mypr(1:mm)).gt.0)) then
        call priref (nr,
     +'M.L. Sierk & G.J. Kleywegt (2004). Deja vu all over again:',
     +'finding and analyzing protein structure similarities.',
     +'Structure 12, 2103-2111.',
     +pubpre(1:lp1)//'15576025'//pubpos(1:lp2),
     + ' ',' ')
      end if
c
c ... 2005 ..............................................................
c
      if (index('|DATAMAN|',mypr(1:mm)).gt.0) then
        call priref (nr,
     +'E. Jakobsson, J. Nilsson, U. Kallstrmm, D. Ogg & G.J.',
     +'G.J. Kleywegt (2005). Crystallization of a truncated',
     +'soluble human semicarbazide-sensitive amine oxidase.',
     +'Acta Cryst. F61, 274-278.',
     +pubpre(1:lp1)//'16511016'//pubpos(1:lp2),
     + ' ')
      end if
c
      if (lspaz) then
        call priref (nr,
     +'M. Novotny & G.J. Kleywegt (2005). A survey of',
     +'left-handed helices in protein structures.',
     +'J. Mol. Biol. 347, 231-241.',
     +pubpre(1:lp1)//'15740737'//pubpos(1:lp2),
     + ' ',' ')
      end if
c
c ... 2006 ..............................................................
c
      if (index('|MOLEMAN2|OOPS|OOPS2|',mypr(1:mm)) .gt. 0) then
        call priref (nr,
     +'G.J. Kleywegt (2006). Quality control and validation.',
     +'Methods in Molecular Biology, 364, 255-272.',
     +pubpre(1:lp1)//'17172770'//pubpos(1:lp2),
     +' ',' ',' ')
      end if
c
c ... 2007 ..............................................................
c
      if (index('|LIGCOM|XPLO2D|HETZE|',mypr(1:mm)) .gt. 0) then
        call priref (nr,
     +'G.J. Kleywegt (2007). Crystallographic refinement of ligand',
     +'complexes. Acta Cryst. D63 (CCP4 Proceedings), 270-274.',
     +pubpre(1:lp1)//'17164531'//pubpos(1:lp2),
     +' ',' ',' ')
      end if
c
      if (index('|LIGCOM|',mypr(1:mm)) .gt. 0) then
        call priref (nr,
     +'G.J. Kleywegt & M.R. Harris (2007). ValLigURL: a server for',
     +'ligand-structure comparison and validation. Acta Cryst. D63,',
     +'935-938.',
     +pubpre(1:lp1)//'17642521'//pubpos(1:lp2),
     + ' ',' ')
      end if
c
      if (index('|DATAMAN|',mypr(1:mm)) .gt. 0) then
        call priref (nr,
     +'G.J. Kleywegt (2007). Separating model optimization and model',
     +'validation in statistical cross-validation as applied to',
     +'crystallography. Acta Cryst. D63, 939-940.',
     +pubpre(1:lp1)//'17704561'//pubpos(1:lp2),
     + ' ',' ')
      end if




c
c ... to be published, in the press, etc. ...............................
c




c
c ... unpublished program
c
      if (nr .le. 0) then
        call priref (nr,
     +'G.J. Kleywegt (1992-2006).',
     +'Uppsala University, Uppsala, Sweden.',
     +'Unpublished program.',
     +' ',' ',' ')
      end if
c
c ... Int Tables chapter covers most software
c
c234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c        1         2         3         4         5         6         7         8         9        10
c
      call priref (nr,
     +'Kleywegt, G.J., Zou, J.Y., Kjeldgaard, M. & Jones, T.A. (2001).',
     +'Around O. In: "International Tables for Crystallography, Vol. F.'
     +,'Crystallography of Biological Macromolecules" (Rossmann, M.G.',
     +'& Arnold, E., Editors). Chapter 17.1, pp. 353-356, 366-367.',
     +'Dordrecht: Kluwer Academic Publishers, The Netherlands.',
     +' ')
c
c ... point to Web site
c
      write (*,6000)
 6000 format (/
     +  ' ==> For manuals and up-to-date references, visit:'/
     +  ' ==>     http://xray.bmc.uu.se/usf'/
     +  ' ==> For reprints, visit:'/
     +  ' ==>     http://xray.bmc.uu.se/gerard'/
     +  ' ==> For downloading up-to-date versions, visit:'/
     +  ' ==>     ftp://xray.bmc.uu.se/pub/gerard')
c
      return
      end
