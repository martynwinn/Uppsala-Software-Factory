REMA
REMA *************************************************************
REMA
REMA maprop.polar - "polarity" library used by program MAPROP
REMA              - code: polar +1, hydrophobic -1, rest 0
REMA
REMA v. 0.1 - Gerard Kleijwegt @ 940317
REMA
REMA
REMA use MAPROP, MAPMAN and O as follows:
REMA power = 1, constant = 1, cuton = 1, cutoff = ~7
REMA use combination method SUM to effectively calculate
REMA the 'polarity' or 'hydrophobicity' at each point
REMA mappage the NEWEZD map with MAPMAN; contour at e.g. +2 in O
REMA with colour blue to show polar patches
REMA multiply map with -1 in MAPMAN; contour at +2 in O
REMA with colour green to show hydrophobic patches
REMA
REMA *** elements *** TWO characters; RIGHT justified
REMA
ELEM ' N' 1.0
ELEM ' C' -1.0
ELEM ' O' 1.0
ELEM ' S' 0.0
ELEM ' H' 0.0
ELEM ' P' 0.0
REMA
ELEM 'NA' 1.0
ELEM 'MG' 1.0
ELEM 'CL' 1.0
ELEM 'CA' 1.0
ELEM 'MN' 1.0
ELEM 'FE' 1.0
ELEM 'ZN' 1.0
ELEM 'CD' 1.0
REMA
REMA *************************************************************
REMA
REMA *** atom types *** FOUR characters; element is first TWO
REMA
ATOM ' N  '  1.0
ATOM ' C  '  0.0
ATOM ' O  '  1.0
ATOM ' CA ' -1.0
REMA
REMA *************************************************************
REMA
REMA *** residue/atom combinations *** EIGHT chracaters: RRR*AAAA
REMA first three=residue, four=anything, five-eight=atom name
REMA
SPAT 'ARG* CD '  0.0
SPAT 'ARG* CZ '  0.0
REMA
SPAT 'ASN* CG '  0.0
REMA
SPAT 'ASP* CG '  0.0
REMA
SPAT 'GLN* CD '  0.0
REMA
SPAT 'GLU* CG '  0.0
SPAT 'GLU* CD '  0.0
REMA
SPAT 'HIS* CG '  0.0
SPAT 'HIS* CD2'  0.0
SPAT 'HIS* CE1'  0.0
REMA
SPAT 'LYS* CE '  0.0
REMA
SPAT 'SER* CB '  0.0
REMA
SPAT 'THR* CB '  0.0
REMA
SPAT 'TYR* CZ '  0.0
REMA
REMA *************************************************************
REMA
REMA *** allowed residue types ***
REMA
RESI 'ALA'
RESI 'ARG'
RESI 'ASN'
RESI 'ASP'
RESI 'CYS'
RESI 'GLN'
RESI 'GLU'
RESI 'GLY'
RESI 'HIS'
RESI 'ILE'
RESI 'LEU'
RESI 'LYS'
RESI 'MET'
RESI 'PHE'
RESI 'PRO'
RESI 'SER'
RESI 'THR'
RESI 'TRP'
RESI 'TYR'
RESI 'VAL'
RESI 'CPR'
RESI 'CYH'
RESI 'PYR'
RESI 'PCA'
REMA
REMA *************************************************************
REMA
END 
