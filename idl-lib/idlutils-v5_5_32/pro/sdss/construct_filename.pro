;+
; NAME:
;   construct_filename
; PURPOSE:
;   Return an SDSS file name given its type and identifiers
; CALLING SEQUENCE:
;   filename= construct_filename(config, filetype [, function_base=, keywords])
; INPUTS:
;   config - configuration file object (from mg_read_config)
;   filetype - name of file type
; OPTIONAL INPUTS:
;   function_base - base for function calls (default '')
;   keywords - named keywords with values necessary to define 
;              file location on disk
; COMMENTS:
;   Usually called by a wrapper (e.g. sdss_filename)
; REVISION HISTORY:
;   2015-06-02: written, MRB, NYU
;-
function construct_filename, config, filetype, function_base= function_base, $
                             _EXTRA=keywords

if(n_elements(function_base) eq 0) then $
  function_base=''

filename= config->get(filetype, section='paths')

;; handle python format variables
filename= parse_string_format(filename, _EXTRA=keywords)

;; handle environmental variables
wordchars='\$([A-Za-z0-9_]*)'
indx= stregex(filename, wordchars, length=length)
while(indx ne -1) do begin
    mid=strmid(filename, indx+1, length-1)
    value= getenv(mid)
    filename= strmid(filename, 0, indx)+value+strmid(filename, indx+length)
    indx= stregex(filename, wordchars, length=length)
endwhile

;; execute functions
wordchars='\%([A-Za-z0-9_]*)'
indx= stregex(filename, wordchars, length=length)
while(indx ne -1) do begin
    mid=strmid(filename, indx+1, length-1)
    value= call_function(function_base+mid, _EXTRA=keywords)
    filename= strmid(filename, 0, indx)+value+strmid(filename, indx+length)
    indx= stregex(filename, wordchars, length=length)
endwhile

return, filename

end

