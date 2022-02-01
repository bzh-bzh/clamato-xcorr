;+
; NAME:
;   print_maskbits
;
; PURPOSE:
;   Print the maskbits table for debugging purposes.
;
; CALLING SEQUENCE:
;   print_maskbits, [flagprefix], [/alias]
;
; INPUTS:
;   flagprefix - Flag name (scalar string).  If not set, everything will
;                be printed.
;
; OPTIONAL INPUTS:
;   alias      - If set, print any aliases associated with the flagprefix.
;
; PROCEDURES CALLED:
;   sdss_maskbits
;
; DATA FILES:
;   $IDLUTILS_DIR/data/sdss/sdssMaskbits.par
;
; REVISION HISTORY:
;   2014-03-27: written, BAW, NYU
;-
PRO print_maskbits, flagprefix, alias=alias
    ;
    ; Declare a common block so that the mask names are remembered between calls.
    ;
    COMMON com_maskbits, maskbits, masktype, maskalias
    ;
    ; Ensure the common block is set
    ;
    sdss_maskbits
    ;
    ;
    ;
    table = ['']
    IF KEYWORD_SET(flagprefix) THEN BEGIN
        ;
        ; Check for aliases
        ;
        aliasmatch = WHERE(STRUPCASE(flagprefix[0]) EQ STRUPCASE(maskalias.alias),ct)
        IF (ct GT 0) THEN BEGIN
            prefix = maskalias[aliasmatch[0]].flag
        ENDIF ELSE prefix = flagprefix[0]
        w = WHERE(maskbits.flag EQ prefix,ct)
        maxdesc = 0
        foo = maskbits[w]
        table = [table,'LABEL                          | BIT | DESCRIPTION']
        table = [table,'-------------------------------+-----+-']
        FOR i = 0, ct-1 DO BEGIN
            desc = STRLEN(foo[i].description)
            IF desc GT maxdesc THEN maxdesc = desc
            table=[table,STRING(foo[i].label, foo[i].bit, foo[i].description,$
                FORMAT='(A-30," | ",I3," | ",A-)')]
        ENDFOR
        FOR j = 0, maxdesc-1 DO table[2] = table[2] + '-'
        IF KEYWORD_SET(alias) THEN BEGIN
            IF prefix NE flagprefix[0] THEN $
                table = [table, '', STRUPCASE(flagprefix[0]) + ' is an alias for ' $
                    + STRUPCASE(prefix) + '.'] $
            ELSE BEGIN
                w = WHERE(prefix EQ maskalias.flag,ct)
                IF ct GT 0 THEN $
                    table = [table, '', STRUPCASE(prefix) + ' has these aliases: ' $
                        + STRJOIN(maskalias[w].alias,', ') + '.']
            ENDELSE
        ENDIF
    ENDIF ELSE BEGIN
        maxdesc = 0
        table = [table,'FLAG            | LABEL                          | BIT | DESCRIPTION']
        table = [table,'----------------+--------------------------------+-----+-']
        FOR k = 0, N_ELEMENTS(maskbits)-1 DO BEGIN
            desc = STRLEN(maskbits[k].description)
            IF desc GT maxdesc THEN maxdesc = desc
            table=[table,STRING(maskbits[k].flag, maskbits[k].label, maskbits[k].bit, maskbits[k].description,$
                FORMAT='(A-20," | ",A-30," | ",I3," | ",A-)')]
        ENDFOR
        FOR j = 0, maxdesc-1 DO table[2] = table[2] + '-'
    ENDELSE
    PRINT, STRJOIN(table,STRING(10B))
END
