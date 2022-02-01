pro test_parse_string_format

aa_string='AA'
bb_string='QWERTYUIIOP'

aa_int=-10
bb_int=10
cc_int=12345000
dd_int=-12324550

aa_float=-10.0234000d
bb_float=10.00000d
cc_float=12345000.230000d
dd_float=-12324550.00001d

openr, unit, 'test.txt', /get_lun
openw, unit2, 'test-idl.txt', /get_lun

while(NOT eof(unit)) do begin
    pattern=''
    readf, unit, pattern
    printf, unit2, pattern+' '+ $
      parse_string_format(pattern, aa_string=aa_string, $
                          bb_string=bb_string, $
                          aa_int=aa_int, bb_int=bb_int, $
                          cc_int=cc_int, dd_int=dd_int, $
                          aa_float=aa_float, bb_float=bb_float, $
                          cc_float=cc_float, dd_float=dd_float)
endwhile

free_lun, unit
free_lun, unit2

end
