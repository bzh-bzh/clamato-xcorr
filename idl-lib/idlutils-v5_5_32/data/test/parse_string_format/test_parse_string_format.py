#!/usr/bin/env python

aa_string='AA'
bb_string='QWERTYUIIOP'

aa_int=-10
bb_int=10
cc_int=12345000
dd_int=-12324550

aa_float=-10.0234000
bb_float=10.00000
cc_float=12345000.230000
dd_float=-12324550.00001

fp= open('test.txt')
ofp= open('test-python.txt','w')
for line in fp.readlines():
    line=line.rstrip('\n')
    out= line.format(aa_string=aa_string, bb_string=bb_string,
                     aa_int=aa_int, bb_int=bb_int, cc_int=cc_int, dd_int=dd_int,
                     aa_float=aa_float, bb_float=bb_float, cc_float=cc_float,
                     dd_float=dd_float )
    ofp.write(line+' '+out+'\n')

fp.close()
ofp.close()
