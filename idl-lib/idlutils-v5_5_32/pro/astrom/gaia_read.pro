;+
; NAME:
;       gaia_read
; PURPOSE:
;       Read GAIA stars near a given RA, dec
;
; INPUTS:
;   racen     - RA of region center (J2000)    [degrees]
;   deccen    - dec of region center (J2000)   [degrees]
;   rad       - radius of region               [degrees]
;
; OPTIONAL INPUTS:
;   dr        - data release name (default 'dr1')
;
; OUTPUTS:
;   result    - data structure defined by FITS file
;
; COMMENTS:    
;   Assumes data in $GAIA_DATA/[dr]/gaia_source/fits_sorted
;   produced by gaia_setup.py. 
;
; REVISION HISTORY
;   2016-Oct-20   MRB, NYU
;
;------------------------------------------------------------------------  
function gaia_read, racen, deccen, rad, dr=dr

  common com_gaia_read, gaia_index

  if(NOT keyword_set(dr)) then $
     dr = 'dr1'

  fitspath = concat_dir(getenv('GAIA_DATA'), dr+'/gaia_source/fits_sorted')
  if fitspath EQ 'slice' then begin 
     splog, 'You must set your $GAIA_DATA environment variable!'
     stop
  endif

  if(n_tags(gaia_index) eq 0) then begin
     gaia_index = mrdfits(fitspath+'/gaia-index.fits',1)
  endif

  spherematch, gaia_index.racen, gaia_index.deccen, racen, deccen, $
               rad + 1., m1, m2, d12, maxmatch=0

  uu = uniq(m1,sort(m1))
  m1 = m1[uu]

  if(m1[0] eq -1) then $
     message, 'GAIA data should exist everywhere.'
  
  data = 0
  for i = 0L, n_elements(m1)-1L do begin
     indx = m1[i]
     ira = strtrim(string(gaia_index[indx].ira),2)
     idec = strtrim(string(gaia_index[indx].idec),2)
     gaia_file = fitspath+'/gaia-'+dr+'-'+ira+'-'+idec+'.fits'
     tmp_data = mrdfits(gaia_file, 1, /silent)
     if(n_tags(tmp_data) gt 0) then begin
        spherematch, tmp_data.ra, tmp_data.dec, racen, deccen, rad, $
                     cm1, cm2, d12, maxmatch=0
        if(cm1[0] ne -1) then begin
           if(n_tags(data) gt 0) then begin
              data = [data, tmp_data[cm1]] 
           endif else begin
              data = tmp_data[cm1]
           endelse
        endif
     endif
  endfor

  return, data
end
