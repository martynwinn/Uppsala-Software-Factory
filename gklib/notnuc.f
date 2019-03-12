c
c
c
      logical function notnuc (restyp)
c
c ... NOTNUC - return .TRUE. if residue type RESTYP is known NOT
c       to be a nucleotide that can be part of a DNA/RNA chain
c
      implicit none
c
      character restyp*3
c
code ...
c
      notnuc = .true.
c
      if (index('|ATP|GTP|CTP|UTP|TTP|PRH|',restyp) .gt. 0) return
      if (index('|TIA|IMU|ZEB|IMP|XMP|NCN|',restyp) .gt. 0) return
      if (index('|2AM|UVC|OMP|TNT|DM1|DM2|',restyp) .gt. 0) return
      if (index('|DM3|DM4|DM5|DM6|DM7|DM8|',restyp) .gt. 0) return
      if (index('|DM9|DMM|  I|CRH|PET|NGM|',restyp) .gt. 0) return
      if (index('|PAP|ATR|BDA|ABP|DAP|MAR|',restyp) .gt. 0) return
      if (index('|BRN|PNT|CPH|BAP|D24|FMN|',restyp) .gt. 0) return
      if (index('|RBF|FAD|PIN|ANP| PU|D18|',restyp) .gt. 0) return
      if (index('|D34|D35|FMC|RIB|PUA|CPG|',restyp) .gt. 0) return
      if (index('|DCG|MNG|PSO|BGF|CMD|D19|',restyp) .gt. 0) return
      if (index('|NCS|APR|BRO|DAG|NOD|DXA|',restyp) .gt. 0) return
      if (index('|PIQ|NGU|DTP|DCP|AGS|NGP|',restyp) .gt. 0) return
      if (index('|1RB|7RP|7RA|A2P|137|G3D|',restyp) .gt. 0) return
      if (index('|ABF|OAD|THM|BLS|CB2|RP5|',restyp) .gt. 0) return
      if (index('|PUY|2DF|MDR|MTP|GNP|NRI|',restyp) .gt. 0) return
      if (index('|ROB|R1P|PRT|5UD|AR4|3D1|',restyp) .gt. 0) return
      if (index('|SNI|CG2|AIS|25A|XR2|PAX|',restyp) .gt. 0) return
      if (index('|OR5|LOF|UA3|UMF|SR1|SAH|',restyp) .gt. 0) return
      if (index('|AAB|HAM|M7G|SFG|MGT|CPR|',restyp) .gt. 0) return
      if (index('|FYA|FUP|HCI|AMO|BOG|GPL|',restyp) .gt. 0) return
      if (index('|SAM|A3P|THP|LAD|NPF|PTP|',restyp) .gt. 0) return
      if (index('|QSI|TSB|APC|MRC|SSA|NEA|',restyp) .gt. 0) return
      if (index('|VAA|GAP|GMC|LKC|LCA|LCC|',restyp) .gt. 0) return
      if (index('|LMS|PSD|ATG|PED|TTE|TYM|',restyp) .gt. 0) return
      if (index('|MGP|ADI|DG3|ILA|HPD|SUC|',restyp) .gt. 0) return
      if (index('|AV2|ABT|IDP|OIP|CH1|U3H|',restyp) .gt. 0) return
      if (index('|RVP|UBB|GNG|P5A|A5A|2AD|',restyp) .gt. 0) return
      if (index('|DUR|CDM|NHE|MSP|MOD|HDF|',restyp) .gt. 0) return
      if (index('|DAD|3DA|DDY|AHX|HCX|BGL|',restyp) .gt. 0) return
      if (index('|FOX|A3S|UFP|2PR|MXA|YSA|',restyp) .gt. 0) return
      if (index('|HXC|LCH|VMS|SRP|2VA|HC4|',restyp) .gt. 0) return
      if (index('|SAP|8HG|G2P|D5M|SA8|FDA|',restyp) .gt. 0) return
      if (index('|PYO|A1P|TYA|MG7|UNK|AZT|',restyp) .gt. 0) return
      if (index('|WWF|O2C|THV|THW|TZD|GCQ|',restyp) .gt. 0) return
      if (index('|FLQ|4AD|NSS|DUP|AIR|ICR|',restyp) .gt. 0) return
      if (index('|TPP|GSU|N3E|2AU|3PD|U33|',restyp) .gt. 0) return
      if (index('|ABG|23T|APK|K05|GFF|GGH|',restyp) .gt. 0) return
      if (index('|C5P|IRF| TT|DUT|CXR|CGR|',restyp) .gt. 0) return
      if (index('|NGD|AQP|BNR|NAE|NAQ|NDC|',restyp) .gt. 0) return
      if (index('|PPY|MAL|GDX|RBZ|PMD|44D|',restyp) .gt. 0) return
      if (index('|DCM|AEN|UPG|UGA|BIZ|UFG|',restyp) .gt. 0) return
      if (index('|YMP|UD1|TIZ|YLY|ANZ|BPF|',restyp) .gt. 0) return
      if (index('|UMA|LLP|TCB|145|BNG|CDF|',restyp) .gt. 0) return
      if (index('|BVD|MTA|LBT|3AD|ADV|3AT|',restyp) .gt. 0) return
      if (index('|HSS|AZZ|N5P|GFH|2OA|RFZ|',restyp) .gt. 0) return
      if (index('|TM2|WSA|N5C|2BA|GGT|F2A|',restyp) .gt. 0) return
      if (index('|AD3|2FA|ARJ|RDD|Q22|3DH|',restyp) .gt. 0) return
      if (index('|RF5|RGT|EEM|CTC|AIF||',restyp) .gt. 0) return
c
c      if (index('|||||||',restyp) .gt. 0) return
c
      notnuc = .false.
      return
c
c  999 continue
c      notnuc = .true.
c
      return
      end
