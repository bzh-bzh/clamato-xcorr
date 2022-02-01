;-----------------------------------------------------------------------
;+
; NAME:
;   djs_laxisgen
;
; PURPOSE:
;   Return a longword integer array with the specified dimensions.
;   Each element of the array is set equal to its index number along
;   the dimension IAXIS.
;
; CALLING SEQUENCE:
;   result = djs_laxisgen( dimens, [ iaxis=iaxis ] )
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
; NOTES:
;   For any number of dimensions larger than 1, this routine simply calls
;   djs_laxisnum().
;
; PROCEDURES CALLED:
;   djs_laxisnum()
;
; REVISION HISTORY:
;   Written by D. Schlegel, 7 Oct 1997, Durham
;   Modified 12 May 1998 to pass one vector with all dimensions.
;   Removed code that was redunant with djs_laxisnum(), 2015-06-24, BAW.
;-
;-----------------------------------------------------------------------
FUNCTION djs_laxisgen, dimens, iaxis=iaxis
    ;
    ; Need one parameter
    ;
    IF N_PARAMS() LT 1 THEN BEGIN
        PRINT, 'Syntax - result = djs_laxisgen( dimens, [iaxis=iaxis] )'
        RETURN, -1
    ENDIF

    IF ~KEYWORD_SET(iaxis) THEN iaxis = 0

    ndimen = N_ELEMENTS(dimens)
    naxis = LONG(dimens) ; convert to type LONG

    IF ndimen EQ 1 THEN RETURN, LINDGEN(naxis[0]) $
    ELSE RETURN, djs_laxisnum(dimens,iaxis=iaxis)
END
;-----------------------------------------------------------------------
