;+
; NAME:
;   wise_ascii2fits
;
; PURPOSE:
;   Reformat the WISE catalog release as FITS binary table files
;
; CALLING SEQUENCE:
;   wise_ascii2fits
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The ASCII files in $WISE_DIR/ascii are converted to FITS binary table
;   files in $WISE_DIR/fits, with the same naming convention
;
; EXAMPLES:
;
; BUGS:
;    NULL entries for integer values are not written as NULL in the FITS
;    files, as I don't know how to do that with MWRFITS.
;    NULL entries for floating-point values are written as NaN.
;
; DATA FILES:
;    The input files are:
;      $WISE_DIR/wise*-cat-schema.txt
;      $WISE_DIR/ascii/wise*-cat-part??.gz
;    The output files are:
;      $WISE_DIR/fits/wise*-cat-part??.fits
;      $WISE_DIR/fits/wise*-cat-part??-radec.fits
;    where the latter file contains the RA,Dec coordinates per object only.
;
; PROCEDURES CALLED:
;   wise_ascii2fits_part
;
; REVISION HISTORY:
;   2016-02-18   Move individual file processing to wise_ascii2fits_part.
;   18-Apr-2012  Written by D. Schlegel, LBL
;-
;
;
;
pro wise_ascii2fits
    ;
    ; Check inputs.
    ;
    wiseDir = GETENV('WISE_DIR')
    IF ~KEYWORD_SET(wiseDir) THEN BEGIN
        PRINT, 'WISE_DIR is not defined!'
        RETURN
    ENDIF
    indir = wiseDir + '/ascii'
    outdir = wiseDir + '/fits'
    IF FILE_SEARCH(indir, /TEST_DIRECTORY) EQ '' THEN BEGIN
        PRINT, 'Input directory does not exist: '+indir
        RETURN
    ENDIF
    IF FILE_SEARCH(outdir, /TEST_DIRECTORY) EQ '' THEN BEGIN
        PRINT, 'Output directory does not exist: '+outdir
        RETURN
    ENDIF
    files = FILE_SEARCH(indir, 'wise-*-cat-part??.gz', COUNT=nfile)
    IF (nfile EQ 0) THEN BEGIN
        PRINT, 'No input files found!'
        RETURN
    ENDIF
    ;
    ; Loop over files.
    ;
    FOR ifile=0, nfile-1 DO BEGIN
        parts = STREGEX(fileandpath(files[ifile]), $
            'wise-(allsky|allwise)-cat-part([0-9][0-9])\.gz', $
            /SUBEXPR, /EXTRACT)
        catalog = parts[1]
        part = FIX(parts[2])
        wise_ascii2fits_part, part=part
    ENDFOR
    RETURN
END
