;-----------------------------------------------------------------------
;+
; NAME:
;   mrdfits_many
;
; PURPOSE:
;   Read in and append many FITS binary tables to a single structure.
;
; CALLING SEQUENCE:
;   res = mrdfits_many(files, [ exten, _EXTRA= ])
;
; INPUTS:
;   files     - File names, or expression for matching file names
;
; OPTIONAL INPUTS:
;   exten     - Extension number; default to 1
;   _EXTRA    - Keywords for MRDFITS
;
; OUTPUTS:
;   res       - Output structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This function will be fast if the structures are identical and
;   can therefore be copied in memory.  If they are not identical, then
;   the first non-empty file defines the structure, and matching column names
;   from other structures are copied using COPY_STRUCT_INX.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   copy_struct_inx
;   mrdfits()
;
; REVISION HISTORY:
;   14-Jun-2017  Written by David Schlegel, Berkeley Lab
;-
;------------------------------------------------------------------------------
; Return 1 if the two structures are the same and can therefore be appended
; Each should be a scalar structure
function mrdfits_many_same, res1, res2

   ntag1 = n_tags(res1)
   ntag2 = n_tags(res2)
   if (ntag1 NE ntag2) then return, 0B

   tag1 = tag_names(res1)
   tag2 = tag_names(res2)
   if (total(tag1 NE tag2) GT 0) then return, 0B

   for i=0, ntag1-1 do $
    if (size(res1.(i),/type) NE size(res2.(i),/type)) then return, 0B

   for i=0, ntag1-1 do $
    if (n_elements(res1.(i)) NE n_elements(res2.(i))) then return, 0B

   return, 1B
end
;------------------------------------------------------------------------------
function mrdfits_many, files, exten1, _EXTRA=EXTRA

   ; Need at least 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - res = mrdfits_many(files)'
      return, 0
   endif
   if (n_elements(exten1) GT 0) then exten = exten1[0] $
    else exten = 1

   ;----------
   ; Find the file names

   if (n_elements(files) EQ 1) then $
    allfiles = file_search(files, count=nfile) $
   else $
    allfiles = files
   nfile = n_elements(allfiles)

   ; Count the total number of elements in all files
   nper = lonarr(nfile)
   for i=0L, nfile-1L do begin
print,i,string(13b),format='(i6,a,$)'
      hdr1 = headfits(allfiles[i], exten=exten)
      if (size(hdr1,/tname) EQ 'STRING') then $
       nper[i] = sxpar(hdr1, 'NAXIS2')
      if (nper[i] GT 0 AND keyword_set(res1) EQ 0) then $
       res1 = mrdfits(allfiles[i], exten, range=[0,0], _EXTRA=EXTRA, /silent)
   endfor
   ntot = total(nper)

   ; Create the blank output structure
   if (keyword_set(res1) EQ 0) then return, 0
   struct_assign, {junk:0}, res1
   res = replicate(res1, ntot)

   ; Now loop over and read the data
   ntot = 0L
   for i=0L, nfile-1L do begin
print,i,string(13b),format='(i6,a,$)'
      if (nper[i] GT 0) then begin
         res1 = mrdfits(allfiles[i], exten, _EXTRA=EXTRA, /silent)
         if (n_elements(res1) NE nper[i]) then $
          message, 'File length changed between access times '+allfiles[i]
         ; Directly copy in the structure contents if identical, otherwise
         ; do the slower copy of only matching elements
         if (mrdfits_many_same(res[0],res1[0])) then $
          res[ntot:ntot+nper[i]-1] = res1 $
         else $
          copy_struct_inx, res1, res, index_from=lindgen(nper[i]), $
           index_to=ntot+lindgen(nper[i])
         ntot += nper[i]
      endif
   endfor
print

   return, res
end 
;------------------------------------------------------------------------------
