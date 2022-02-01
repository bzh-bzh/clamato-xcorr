pro test_types

readcol, 'test_types.txt', f='(a)', types

for i=0L, n_elements(types)-1L do begin 
    path1= sdss_name(types[i], 1336, 4, 100, rerun=301, filter='g')
    path2= sdss_filename(types[i], run=1336, rerun=301, camcol=4, field=100, filter='g', type='sky') 
    if(path1 ne path2) then begin
        splog, types[i]
        splog, ' '+path1
        splog, ' '+path2
    endif
endfor 

end
