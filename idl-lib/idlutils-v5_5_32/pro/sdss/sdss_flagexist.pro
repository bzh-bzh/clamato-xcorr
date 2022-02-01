;+
; NAME:
;   sdss_flagexist
;
; PURPOSE:
;   Returns whether a flag exists
;
; CALLING SEQUENCE:
;   exist= sdss_flagexist(flagprefix, label, flagexist=, whichexist=)
;
; INPUTS:
;   flagprefix - Flag name (scalar string).
;   label      - String name(s) corresponding to each non-zero bit in FLAGVALUE.
;
; OUTPUTS:
;   exist - 1 if label exists (or all labels exist) for this flag (0 otherwise)
;   flagexist - 1 if this flag exists (0 otherwise)
;   whichexist - byte array showing which individual labels exist.
;
; PROCEDURES CALLED:
;   sdss_maskbits
;
; DATA FILES:
;   $IDLUTILS_DIR/data/sdss/sdssMaskbits.par
;
; REVISION HISTORY:
;   2010-07-02: made by modifying sdss_flagval, MRB, NYU
;-
;------------------------------------------------------------------------------
FUNCTION sdss_flagexist, flagprefix, inlabel, flagexist=flagexist, whichexist=whichexist
    ;
    ; Declare a common block so that the mask names are remembered between calls.
    ;
    COMMON com_maskbits, maskbits, masktype, maskalias
    ;
    ; Usage
    ;
    IF (N_PARAMS() NE 2 OR N_ELEMENTS(flagprefix) NE 1) THEN $
        MESSAGE, 'Syntax - exist = sdss_flagexist(flagprefix, label [, flagexist=])'
    ;
    ; Generate a list of all non-blank labels as a string array
    ;
    alllabel = STRSPLIT(inlabel[0], /EXTRACT)
    FOR i=1, N_ELEMENTS(inlabel)-1 DO $
        alllabel = [alllabel, STRSPLIT(inlabel[i], /EXTRACT)]
    ilabel = WHERE(alllabel NE '', nlabel)
    IF (nlabel EQ 0) THEN RETURN, 0B
    alllabel = alllabel[ilabel]
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
    ; Check the flag
    ;
    imatch=WHERE(STRUPCASE(prefix) EQ STRUPCASE(maskbits.flag), ct)
    IF (ct NE 0) THEN flagexist=1B ELSE BEGIN
        flagexist=0B
        RETURN, 0B
    ENDELSE
    ;
    ; Check the labels
    ;
    whichexist = BYTARR(nlabel)
    FOR ilabel=0, nlabel-1 DO BEGIN
        imatch = WHERE((STRUPCASE(prefix) EQ STRUPCASE(maskbits.flag)) $
            AND (STRUPCASE(alllabel[ilabel]) EQ STRUPCASE(maskbits.label)), ct)
        whichexist[ilabel] = (ct EQ 1)
    ENDFOR
    exist = (TOTAL(whichexist) EQ nlabel)
    RETURN, exist
END
