;+
; NAME:
;   write_ecsv
;
; PURPOSE:
;   Write ECSV file
;
; CALLING SEQUENCE:
;   write_ecsv, filename, pdata, [ description=, unit=, extname= ]
;
; INPUTS:
;   filename   - Output file name
;   pdata      - Data structure array
;
; OPTIONAL INPUTS:
;   description- String array with description of each column
;   unit       - String array with units of each column
;   extname    - If set, then use this for the structure name
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Not all data types are supported.
;
; REVISION HISTORY:
;   04-Nov-2016  Written by D. Schlegel, Berkeley Lab
;-
;------------------------------------------------------------------------------
pro write_ecsv, filename, pdata, description=description1, unit=unit1, $
 extname=extname1

   tname = ['string', 'int16', 'int16', 'int32', 'int64', 'float32', 'float64', 'uint16']
   idlname = ['STRING', 'BYTE', 'INT', 'LONG', 'LONG64', 'FLOAT', 'DOUBLE', 'UINT']

   tags = tag_names(pdata)
   ntag = n_elements(tags)
   if (keyword_set(description1)) then description = description1 $
    else description = strarr(ntag)
   if (n_elements(description) NE ntag) then $
    message, 'Wrong number of elements in DESCRIPTION!'
   if (keyword_set(unit1)) then unit = unit1 $
    else unit = strarr(ntag)
   if (n_elements(unit) NE ntag) then $
    message, 'Wrong number of elements in UNIT!'
   if (keyword_set(extname1)) then extname = extname1 $
    else extname = tag_names(pdata, /structure_name)
   if (NOT keyword_set(extname)) then extname = 'STRUCT1'

   get_lun, olun
   openw, olun, filename

   printf, olun, '# %ECSV 0.9'
   printf, olun, '# ---'
   printf, olun, '# datatype:'
   for i=0, ntag-1 do begin
      type1 = tname[where(idlname EQ size(pdata[0].(i), /tname))]
      if (keyword_set(description[i])) then $
       descrip1 = ', description: '+description[i] $
      else descrip1 = ''
      if (keyword_set(unit[i])) then $
       unit1 = ', unit: '+unit1 $
      else unit1 = ''
      printf, olun, '# - {name: '+tags[i]+unit1+', datatype: '+type1+descrip1+'}'
   endfor
   printf, olun, '# meta: !!omap'
   printf, olun, '# - {EXTNAME: '+extname+'}'
   for i=0, ntag-1 do printf, olun, tags[i], format='(a," ",$)'
   printf, olun, ''
   struct_print, pdata, lun=olun, /no_head

   close, olun
   free_lun, olun

   return
end
;------------------------------------------------------------------------------
