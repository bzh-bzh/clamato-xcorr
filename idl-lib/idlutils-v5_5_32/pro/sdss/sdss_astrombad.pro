;+
; NAME:
;   sdss_astrombad
; PURPOSE:
;   For a list of RUN, CAMCOL, FIELD, return whether each field has bad astrometry
; CALLING SEQUENCE:
;   bad= sdss_astrombad(run, camcol, field)
; INPUTS:
;   run, camcol, field - [N] field identifiers
;   rerun - which rerun
; OUTPUTS:
;   bad - 0 for good, 1 for bad
; COMMENTS:
;   Reads data from:
;    $PHOTOLOG_DIR/opfiles/opBadFields.par
;   Note that if there is a problem with one camcol, we assume a
;   problem with all camcols.
; REVISION HISTORY:
;   10-Oct-2012  morphed by M. Blanton, NYU
;------------------------------------------------------------------------------
FUNCTION sdss_astrombad, run, camcol, field
    ;
    ; Declare a common block so that the fields are remembered between calls.
    ;
    COMMON com_astrombad, opbadfields
    ;
    ; Check inputs
    ;
    IF (N_ELEMENTS(run) EQ 0) THEN $
        MESSAGE, 'RUN, CAMCOL, FIELD need at least one element'
    IF ((N_ELEMENTS(run) NE N_ELEMENTS(camcol)) || $
        (N_ELEMENTS(camcol) NE N_ELEMENTS(field))) THEN $
        MESSAGE, 'RUN, CAMCOL, FIELD need the same number of elements'
    IF ~KEYWORD_SET(GETENV('PHOTOLOG_DIR')) THEN $
        MESSAGE, '$PHOTOLOG_DIR not set (photolog product not set up)'
    IF ~KEYWORD_SET(GETENV('BOSS_PHOTOOBJ')) THEN $
        MESSAGE, 'Environmental variable BOSS_PHOTOOBJ must be set'
    ;
    ; Read the file
    ;
    IF (N_TAGS(opbadfields) EQ 0) THEN BEGIN
        opbadfieldsfile= GETENV('PHOTOLOG_DIR')+'/opfiles/opBadfields.par'
        IF ~FILE_TEST(opbadfieldsfile) THEN $
            MESSAGE, 'File required: '+opbadfieldsfile
        opbadfields= yanny_readone(opbadfieldsfile)
    ENDIF
    ;
    ; Find the bad fields
    ;
    iastrom= WHERE((opbadfields.problem EQ 'astrom') OR $
        (opbadfields.problem EQ 'rotator'), nastrom)
    bad= BYTARR(N_ELEMENTS(run))
    FOR i=0L, nastrom-1L DO BEGIN
        irun= WHERE((run EQ opbadfields[iastrom[i]].run) AND $
           (field GE opbadfields[iastrom[i]].firstfield) AND $
           (field LE opbadfields[iastrom[i]].lastfield), nrun)
        IF (nrun GT 0) THEN bad[irun]=1B
    ENDFOR
    RETURN, bad
END
