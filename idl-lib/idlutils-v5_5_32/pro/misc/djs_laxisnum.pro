;-----------------------------------------------------------------------
;+
; NAME:
;   djs_laxisnum
;
; PURPOSE:
;   Return a longword integer array with the specified dimensions.
;   Each element of the array is set equal to its index number in
;   the specified axis.
;
; CALLING SEQUENCE:
;   result = djs_laxisnum( dimens, [ iaxis=iaxis ] )
;
; INPUT:
;   dimens:     Vector of the dimensions for the result.
;               Only up to 3 dimensions can be specified.
;   iaxis:      Axis number to use for indexing RESULT.  The first dimension
;               is axis number 0, the second 1, etc.  Default to 0.
;
; OUTPUTS:
;   result:     Output array
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   15-Jun-2001  Written by D. Schlegel, Princeton
;   2015-06-24 Code audit, BAW.
;-
;-----------------------------------------------------------------------
FUNCTION djs_laxisnum, dimens, iaxis=iaxis
    ;
    ; Need one parameter
    ;
    IF N_PARAMS() LT 1 THEN BEGIN
        PRINT, 'Syntax - result = djs_laxisnum( dimens, [iaxis= ] )'
        RETURN, -1
    ENDIF

    IF ~KEYWORD_SET(iaxis) THEN iaxis = 0

    ndimen = N_ELEMENTS(dimens)
    naxis = LONG(dimens) ; convert to type LONG

    IF iaxis GE ndimen THEN BEGIN
        PRINT, 'Invalid axis selection!'
        RETURN, -1
    ENDIF

    result = MAKE_ARRAY(DIMENSION=naxis, /LONG)

    CASE ndimen OF
        1 : result[*] = 0L
        2 : BEGIN
            CASE iaxis OF
                0 : FOR ii=0L, naxis[0]-1L DO result[ii,*] = ii
                1 : FOR ii=0L, naxis[1]-1L DO result[*,ii] = ii
            ENDCASE
        END
        3 : BEGIN
            CASE iaxis OF
                0 : FOR ii=0L, naxis[0]-1L DO result[ii,*,*] = ii
                1 : FOR ii=0L, naxis[1]-1L DO result[*,ii,*] = ii
                2 : FOR ii=0L, naxis[2]-1L DO result[*,*,ii] = ii
            ENDCASE
        END
        ELSE : BEGIN
            PRINT, ndimen, ' dimensions not supported!'
            result = -1
        END
    ENDCASE

    RETURN, result
END
;-----------------------------------------------------------------------
