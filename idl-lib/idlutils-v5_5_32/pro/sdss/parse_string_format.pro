;+
; NAME:
;   parse_string_format
; PURPOSE:
;   Parse string similarly to python's string format method
; CALLING SEQUENCE:
;   out= parse_string_format(pyformat [, keywords])
; INPUTS:
;   pyformat - string with python3-esque formatting 
; COMMENTS:
;   Primary use case is for construction of SDSS path names using
;    standard python-oriented sdss_paths.ini file. There are a number
;    of limitations of this implementation and its general use 
;    in other situations is warned against.
;   Limitations:
;    - It uses only named, not numbered, keywords 
;    - It handles ONLY types s,d,f in the python convention
;    - It ignores !r and !s directives
;    - It fails if presented with # or , directives
;    - It will fail for floats with > 100 digits
;    - It does not dereference arrays or structures (e.g. {array[4]}
;      or {object.attribute} don't work).
; REVISION HISTORY:
;   2015-06-02: written, MRB, NYU
;-
function py2idl_format, in_pyformat, value

pyformat= in_pyformat

;; ignore !r or !s
if(strmid(pyformat, 0, 1) eq '!') then $
  pyformat= strmid(pyformat, 2)

;; check colon
colon=strmid(pyformat, 0, 1)
if(colon ne ':') then $
  message, 'Illegal format: '+in_pyformat
pyformat= strmid(pyformat, 1)

;; spawn error if # or , is used
if(strmatch(pyformat, '#') ne 0) then $
  message, '"#" not supported'
if(strmatch(pyformat, ',') ne 0) then $
  message, '"," not supported'

;; all types allowed by python
types='bcdeEfFgGnosxX%'

;; parse the formating string
subs= stregex(pyformat, $
              '^([^><=^]?)([><=\^]?)([\+\-\ ]?)(0?)([0-9]*)(\.?[0-9]*)(['+types+']?)$', $
              /sub, /extr)
fill=subs[1]
justify=subs[2]
sign=subs[3]
zero=subs[4]
width=subs[5]
precision=subs[6]
type=subs[7]

;; if zero precedes width parameter then set fill and justify
;; appropriately, unless they are already set
if(zero ne '') then begin
    if(fill eq '') then $
      fill='0'
    if(justify eq '') then $
      justify='='
endif

if(justify eq '') then $
  justify='>'

;; fill with spaces by default
if(fill eq '') then fill=' '

;; make sure width is an integer
if(width eq '') then $
  width=0 $
else $
  width=fix(width)

;; make sure precision is an integer
if(precision eq '') then begin
    precision=-1
endif else begin
    precision=fix(strmid(precision,1))
endelse

;; determine type
if(type eq '') then begin
    type_name= size(value, /tname)
    if(type_name eq 'BYTE' OR $
       strmatch(type_name, '*INT') ne 0 OR $
       strmatch(type_name,'*LONG*') ne 0) then $
      type='d'
    if(type_name eq 'FLOAT' OR $
       type_name eq 'DOUBLE') then $
      type='f'
    if(type_name eq 'STRING') then $
      type='s'
    if(type eq '') then $
      message, 'Type not supported: '+type_name
endif

;; spawn error if not s, d, or g type
if(type ne 's' and type ne 'd' and type ne 'f') then $
  message, 'Only types s, d, and f supported'

;; string type -- deal with justification
if(type eq 's') then begin
    sign_string=''
    if(strlen(value) lt precision OR $
       precision eq -1) then $
      value_string=value $
    else $
      value_string=strmid(value, 0, precision)
endif else begin

    ;; handle sign options
    sign_string=''
    if(sign eq ' ') then begin
        sign_string=string(7b)
        if(justify eq '=') then $
          width= width-1
    endif
    if(sign eq '+') then begin
        sign_string='+'
        if(justify eq '=') then $
          width= width-1
    endif
    if(value lt 0) then begin
        sign_string='-'
        if((sign eq '-' or sign eq '') and $
           justify eq '=') then $
          width= width-1
    endif

    ;; case of integers
    if(type eq 'd') then begin
        value_string= string(abs(value), format='(i0)')
        if(justify ne '=') then begin
            value_string= sign_string+value_string
            sign_string=''
        endif 
    endif

    ;; case of floating point
    if(type eq 'f') then begin
        if(precision eq -1) then $
          precision=6
        if(precision eq 0) then $
          value_string= string(long(abs(value)), format='(i0)') $
        else $
          value_string= strtrim(string(abs(value), $
                                       format='(f100.'+ $
                                       strtrim(string(precision),2)+')'),2)
        if(justify ne '=') then begin
            value_string= sign_string+value_string
            sign_string=''
        endif 
    endif
endelse

;; now make sure it is the right width and justification
if(width gt strlen(value_string)) then begin
    jsign=''
    if(justify eq '^') then $
      message, 'Centered strings not supported'
    if(justify eq '<') then $
      jsign='-'
    if(justify eq '=') then begin
        jsign=''
    endif
    value_string= $
      string(value_string, format='(a'+jsign+strtrim(string(width),2)+')')
    value_string= repstr(value_string, ' ', fill)
endif

;; Add the sign at front
value_string= sign_string+value_string
value_string= repstr(value_string, string(7b), ' ')

return, value_string

end
;
function parse_string_format, in_str, _EXTRA=keywords

str=in_str

word='[A-Za-z0-9_><().!%: =+-]'
match= '(\{'+word+'*\})'

indx=stregex(str, match, length=length)
while(indx ne -1) do begin
    ;; get format string
    mid= strmid(str, indx+1, length-2)

    ;; divide format strong into keyword and any formatting commands
    words= strsplit(mid, ':!', /extr)
    key=words[0]

    ;; gather the value from the keyword list
    itag= tag_indx(keywords, key)
    if(itag eq -1) then begin
        message, 'Required input "'+key+'"'
    endif
    value= keywords.(itag)

    ;; translate format
    format=''
    if(n_elements(words) gt 1) then begin
        iwords= strsplit(mid, ':!')
        pyformat=strmid(mid, iwords[1]-1)
        strvalue=py2idl_format(pyformat, value)
    endif else begin
        strvalue= strtrim(string(value),2)
    endelse
    
    ;; put the value
    str= strmid(str,0, indx)+strvalue+strmid(str, indx+length)
    indx=stregex(str, match, length=length)
endwhile 

return, str

end
