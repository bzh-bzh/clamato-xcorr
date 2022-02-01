;+
; NAME:
;   sdss_flagname
;
; PURPOSE:
;   Return bitmask labels corresponding to bit numbers.
;
; CALLING SEQUENCE:
;   label = sdss_flagname(flagprefix, flagvalue, [ /concat, /silent ] )
;
; INPUTS:
;   flagprefix - Flag name (scalar string).
;   flagvalue  - Signed long with any number of its bits set.
;
; OPTIONAL KEYWORDS:
;   concat     - If set, then concatenate all of the output labels in
;                LABEL into a single whitespace-separated string.
;   silent     - If set, then don't print a warning when there is no bit label
;                corresponding to one of the bit values.
;
; OUTPUTS:
;   label      - String name(s) corresponding to each non-zero bit in FLAGVALUE.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This function is the inverse of SDSS_FLAGVAL().
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   sdss_maskbits
;   splog
;
; DATA FILES:
;   $IDLUTILS_DIR/data/sdss/sdssMaskbits.par
;
; REVISION HISTORY:
;   01-Apr-2002 Written by D. Schlegel, Princeton.
;   19-Aug-2008 Modified by A. Kim.  sdssMaskbits.par now has accompanying
;               information on the associated datatype for the bits.
;               Code modified to read in this information, check the validity
;               of the par file, and return values of the correct type.
;   15-Jan-2014 Simplify by the use of sdss_maskbits()
;-
;------------------------------------------------------------------------------
FUNCTION sdss_flagname, flagprefix, flagvalue, concat=concat, silent=silent
    ;
    ; Declare a common block so that the mask names are remembered between calls.
    ;
    COMMON com_maskbits, maskbits, masktype, maskalias
    ;
    ; Usage
    ;
    IF (N_PARAMS() NE 2 OR N_ELEMENTS(flagprefix) NE 1) THEN $
        MESSAGE, 'Syntax - label = sdss_flagname(flagprefix, flagvalue, [ /concat ] )'
    ;
    ; Ensure the common block is set
    ;
    sdss_maskbits
    ;
    ; Check for aliases
    ;
    aliasmatch = WHERE(STRUPCASE(flagprefix[0]) EQ STRUPCASE(maskalias.alias),ct)
    IF (ct GT 0) THEN prefix = maskalias[aliasmatch[0]].flag $
    ELSE prefix = flagprefix[0]
    ;
    ; Find the match for each non-zero bit.
    ;
    indx = WHERE(djs_int2bin(flagvalue), nret)
    IF (indx[0] EQ -1) THEN BEGIN
        retval = ''
    ENDIF ELSE BEGIN
        retval = STRARR(nret)
        FOR iret=0, nret-1 DO BEGIN
            j = WHERE((STRUPCASE(prefix) EQ STRUPCASE(maskbits.flag)) $
                AND (indx[iret] EQ maskbits.bit))
            IF (j[0] NE -1) THEN retval[iret] = maskbits[j].label $
            ELSE IF ~KEYWORD_SET(silent) THEN $
            splog, 'MESSAGE: Unknown bit ', indx[iret], $
                ' for flag ' + STRUPCASE(prefix)
        ENDFOR
    ENDELSE
    ;
    ; If /CONCAT is set, then concatenate all of the output strings
    ; into a single string separated only by whitespace.
    ;
    IF KEYWORD_SET(concat) THEN retval = STRJOIN(retval,' ')
    RETURN, retval
END
