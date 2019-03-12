c
c ... nuc_subs.f - subroutines for dealing with nucleic acids
c
c ... gerard j kleywegt @ 2000-12-12
c
c ......................................................................
c
      subroutine nuctyp (restyp,itype)
c
c ... map various PDB three-letter codes to the 5 types of nucleotides
c
c ... gerard j kleywegt
c     @ 2000-12-12
c     @ 2001-01-08,05-07,07-26,08-08,10-09,11-12,12-20
c     @ 2002-01-07,02-01,07-03,07-09,08-23,10-01
c     @ 2003-01-27,03-10,04-23,06-05,11-03
c     @ 2004-02-25,03-02,10-04
c     @ 2005-02-01,12-05
c     @ 2006-02-23,04-18,08-02
c     @ 2007-05-07,09-13,11-23
c     @ 2008-06-19,11-27
c
c --------------------------------------------------------------
c
c get list of undefined types from: tail -100 ~/progs/spasm/spana.log
c
C USE:
C
C OLD: alias what 'set x=\!^ ; curl -v "http://ligand-depot.rutgers.edu/pub/"$x/$x".cif.html"  |& head -100 | grep -e "^_chem_comp\." -e "\;" | head -10'
c
C NEW: alias what 'set x=\!^ ; set y=`echo $x | cut -c1-1` ; curl -v "http://ligand-expo.rcsb.org/reports/"$y/$x/$x".cif" |& head -100 | grep -e "^_chem_comp\." -e "\;" | head -10'
c
C then:
c
C what C34
C _chem_comp.name                      'N4-METHYL-2'-DEOXY-CYTIDINE-5'-MONOPHOSPHATE' 
C _chem_comp.type                      'DNA LINKING' 
C
C or:
c
c \rm qqq ; touch qqq
c foreach nuc ( 2BA GGT F2A)
c   (echo "" ; echo $nuc ; what $nuc) >> qqq
c end
c
C put results in text file, search for NON-POLYMER then mon_nstd_parent_comp_id
c and classify as ACGTU - all remaining 'DNA LINKING' and 'RNA LINKING' become
c class 6 (rest are other polymer types, e.g. amino acids or sugars)
c
c or check at:
c
c   http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/vctr/indhet.pl?name=XAL
c
c NOT nucleotides (as part of RNA or DNA strand): 
c 
c
c Not decided yet (or nucleotide unlike AGCTU !): 
C
c  1 = A: TCY S4A
c  2 = G: SDG S4G
c  3 = C: 4PC 10C C4S N5M
c  4 = T: 5PY T5S
c  5 = U: BVP U2N
C  6 = X: 4MF N5I XAR FFD
c  NON: TCB 145 BNG CDF BVD MTA LBT 3AD ADV 3AT HSS AZZ N5P GFH 2OA RFZ TM2 WSA N5C 2BA GGT F2A
c
c Residue types to check: 
c
c Nr of unknown nucleotide types : (         49) 
c Unknown nucleotide types : ( 
c         
c
c
c --------------------------------------------------------------
c
c      data typ3lc /'  A','  G','  C','  T','  U','XXX'/
c                      1     2     3     4     5     6
c
      implicit none
c
      integer itype
c
      character restyp*3
c
code ...
c
      itype = 1
      if (index('|  A| DA| +A|ADE|',restyp) .gt. 0) return
c
      itype = 2
      if (index('|  G| DG| +G|GUA|',restyp) .gt. 0) return
c
      itype = 3
      if (index('|  C| DC| +C|CYT|',restyp) .gt. 0) return
c
      itype = 4
      if (index('|  T| DT| +T|THY|',restyp) .gt. 0) return
c
      itype = 5
      if (index('|  U| +U|URI|',restyp) .gt. 0) return
c
c ... special cases
C
C      if (index('|||||||',restyp) .gt. 0) return
c
c ... AAAAAAAAA
c
      itype = 1
      if (index('|  R|  Y|  E|1MA|SRA|ADP|',restyp) .gt. 0) return
      if (index('|AMP|A23|2AR|AP2|ABM|6HA|',restyp) .gt. 0) return
      if (index('|MA6|12A|AAM|ABR|ABS|ACP|',restyp) .gt. 0) return
      if (index('|RMP|SMP|ADN|6MC|A3A|6MT|',restyp) .gt. 0) return
      if (index('| DA|FA2|2BU|MTU|5AA|A2M|',restyp) .gt. 0) return
      if (index('|A2L|MIA|ADX|MAD|6MA|FA5|',restyp) .gt. 0) return
      if (index('|AD2|AVC|2DA|5CA|7DA|AP7|',restyp) .gt. 0) return
      if (index('| AS|MA7|XAL|XAD|RIA|A44|',restyp) .gt. 0) return
      if (index('|A40|A38|A47|A43|TCY|S4A|',restyp) .gt. 0) return
      if (index('|A8N|2FE|||||',restyp) .gt. 0) return
c
c ... GGGGGGGGG
c
      itype = 2
      if (index('|LGP|GSR|  P|  X|GSS|M2G|',restyp) .gt. 0) return
      if (index('|OMG|7MG|1MG|IGU|2GP|GDP|',restyp) .gt. 0) return
      if (index('|GMP|PGP|5GP|2MG|GN7|SGP|',restyp) .gt. 0) return
      if (index('|FAG|6HG|DFG|8OG|3GP|YYG|',restyp) .gt. 0) return
      if (index('| IG|S6G|DGT|GMS|G2S|BGM|',restyp) .gt. 0) return
      if (index('|G25|8MG| LG| DG|PG7|G2L|',restyp) .gt. 0) return
      if (index('|G4P|DGP|QUO|G7M|GAO|PGN|',restyp) .gt. 0) return
      if (index('|5CG|GH3|MRG|DDG|6OG|7GU|',restyp) .gt. 0) return
      if (index('|8FG|PPW|GNE| GS|XUG|XGL|',restyp) .gt. 0) return
      if (index('|XGU|KAG|2EG|G42|G32|BZG|',restyp) .gt. 0) return
      if (index('|G36|G47|G31|G49|G38|G48|',restyp) .gt. 0) return
      if (index('|GFL|GDR|MG1|SDG|S4G|FMG|',restyp) .gt. 0) return
c
c ... CCCCCCCCC
c
      itype = 3
      if (index('|5MC|MCY|OMC|EDC|IMC|CDP|',restyp) .gt. 0) return
      if (index('|CMP|CCC|55C| CH|6HC|DFC|',restyp) .gt. 0) return
      if (index('|C2P|C3P|CAR| IC|DCZ|GCK|',restyp) .gt. 0) return
      if (index('|DNR|CMR|C2S|5CM|C25|5IC|',restyp) .gt. 0) return
      if (index('| LC|I5C|CBR| DC|TC1|5FC|',restyp) .gt. 0) return
      if (index('|C2L|CSL|A5M|M5M|S4C|5NC|',restyp) .gt. 0) return
      if (index('|DOC|DCT|5PC|TPC|4SC| SC|',restyp) .gt. 0) return
      if (index('|XCL|XCT|C34|C38|C31|C45|',restyp) .gt. 0) return
      if (index('|C49|C36|C46|C42|C43|CFL|',restyp) .gt. 0) return
      if (index('|CBV|1SC|4PC|10C|C4S|N5M|',restyp) .gt. 0) return
c
c ... TTTTTTTTT
c
      itype = 4
      if (index('|THX|5AT|TDP|TMP|6CT|T23|',restyp) .gt. 0) return
      if (index('|TAF|TLC|ATD|6HT|TSP|DRT|',restyp) .gt. 0) return
      if (index('|ADT|SMT|MTR|T2S|TLB|DPB|',restyp) .gt. 0) return
      if (index('|NMT|NMS|2BT|TTM|2NT|2OT|',restyp) .gt. 0) return
      if (index('|EIT|TFE|P2T|2AT|2GT|BOE|',restyp) .gt. 0) return
      if (index('|T3P|D4M|S2M|5IT|TCP|MMT|',restyp) .gt. 0) return
      if (index('|ATL| RT|SPT|TGP|64T|ATM|',restyp) .gt. 0) return
      if (index('|2DT|D3T| TS|CTG|XTL|XTH|',restyp) .gt. 0) return
      if (index('|T41|T49|T32|T48|T38|T39|',restyp) .gt. 0) return
      if (index('|TA3|5PY|T5S||||',restyp) .gt. 0) return
c
c ... UUUUUUUUU
c
      itype = 5
      if (index('|PSU|H2U|5MU|5IU|UDP|UMP|',restyp) .gt. 0) return
      if (index('|U8U|5MD|U5P|U3P|UM3|ADU|',restyp) .gt. 0) return
      if (index('|DHU|70U|2MU|127|125|4SU|',restyp) .gt. 0) return
      if (index('|126|BRU|HEU|UMS|SSU|U2P|',restyp) .gt. 0) return
      if (index('|U25|3ME|GMU|5BU|LHU|DRM|',restyp) .gt. 0) return
      if (index('|OMU|UR3|MNU|U2L|ZDU|FMU|',restyp) .gt. 0) return
      if (index('|P2U|UAR|FHU| DU|PDU|5HU|',restyp) .gt. 0) return
      if (index('| IU|UD5|SUR|S4U|2ST|U34|',restyp) .gt. 0) return
      if (index('|U31|U36|UFR|UCL|5FU|US1|',restyp) .gt. 0) return
      if (index('|BVP|U2N|UFT||||',restyp) .gt. 0) return
c
c ... non-standard nucleotides that (might) occur in DNA/RNA chains
c
      itype = 6
      if (index('| YG| D3|DFT|NP3|MBZ|PYP|',restyp) .gt. 0) return
      if (index('|PYY|3DR|T6A|PPU|MEP|TLN|',restyp) .gt. 0) return
      if (index('|  N|D1P|DPY|DRP|M1G|FMP|',restyp) .gt. 0) return
      if (index('|IRP|P5P|MTR|6MI|2DM|HOL|',restyp) .gt. 0) return
      if (index('|HOB|1AP| DI|LCG|2BD|DXD|',restyp) .gt. 0) return
      if (index('|DXN|HDP|P1P|YRR|PRN|ONE|',restyp) .gt. 0) return
      if (index('|2MA|ASU|EDA|PST|  Z|DDX|',restyp) .gt. 0) return
      if (index('|DRZ|GOM|PQ1|AET|TTD|E1X|',restyp) .gt. 0) return
      if (index('|XCY|PR5|N6G|NF2|T2T|XTY|',restyp) .gt. 0) return
      if (index('|XCS|XGA|XAE|T4S|NCX|PBT|',restyp) .gt. 0) return
      if (index('|6IA|RTP|X4A|B1P|CM0|6MZ|',restyp) .gt. 0) return
      if (index('|NYM|4MF|N5I|XAR|FFD|2FI|',restyp) .gt. 0) return
      if (index('|IRN||||||',restyp) .gt. 0) return
c
c ... unknown
c
      itype = -1
c
      return
      end
