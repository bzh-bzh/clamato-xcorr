;+
; NAME:
;   sdss_maskbits
;
; PURPOSE:
;   Read the maskbits file and set a common block
;
; CALLING SEQUENCE:
;   sdss_maskbits, [maskfile=]
;
; OPTIONAL KEYWORD INPUTS:
;   maskfile - Used to override the built-in maskbits file.
;
; PROCEDURES CALLED:
;   yanny_free
;   yanny_read
;
; DATA FILES:
;   $IDLUTILS_DIR/data/sdss/sdssMaskbits.par
;
; REVISION HISTORY:
;   2014-01-15 Written, BAW, NYU
;-
;------------------------------------------------------------------------------
PRO sdss_maskbits, maskfile=maskfile
    ;
    ; Declare a common block so that the mask names are remembered between calls.
    ;
    COMMON com_maskbits, maskbits, masktype, maskalias
    ;
    ; Allowed mask sizes
    ;
    legaltypes=[16,32,64]
    IF ~KEYWORD_SET(maskbits) THEN BEGIN
        IF ~KEYWORD_SET(maskfile) THEN $
            maskfile = FILEPATH('sdssMaskbits.par', $
            ROOT_DIR=GETENV('IDLUTILS_DIR'), $
            SUBDIRECTORY=['data','sdss'])
        IF ~KEYWORD_SET(maskfile) THEN $
            MESSAGE, 'File with mask bits not found'
        yanny_read, maskfile, pdat, stnames=stnames
        bitsindex=(WHERE(stnames EQ 'MASKBITS'))[0]
        typeindex=(WHERE(stnames EQ 'MASKTYPE'))[0]
        aliasindex=(WHERE(stnames EQ 'MASKALIAS'))[0]
        ;
        ; logic to see of the bits fit within the type
        ;
        ; get all the unique flag names
        ;
        allflags=(*pdat[typeindex]).flag
        alltypeflags=(*pdat[bitsindex]).flag
        allalias=(*pdat[aliasindex]).flag
        ;
        ; check the legality of each of the flags
        ;
        FOR i=0,(SIZE(allflags,/DIM))[0]-1 DO BEGIN
            ;
            ; check to see if this is a legal type
            ;
            itype=(*pdat[typeindex])[i].datatype
            legaltypeindex=WHERE(itype EQ legaltypes,num)
            IF (num EQ 0) THEN $
                MESSAGE, 'Illegal datatype used'
            wflag=WHERE(alltypeflags EQ allflags[i],num)
            ;
            ; find the maximum bit for this flag name
            ;
            IF (num GT 0) THEN BEGIN
                bits=(*pdat[bitsindex])[wflag].bit
                wbad = WHERE(bits GE legaltypes[legaltypeindex], badnum)
                IF (badnum NE 0) THEN $
                    MESSAGE, 'Illegal bit used'
            ENDIF
        ENDFOR
        maskbits = *pdat[bitsindex]
        masktype = *pdat[typeindex]
        maskalias = *pdat[aliasindex]
        yanny_free, pdat
    ENDIF
    RETURN
END
