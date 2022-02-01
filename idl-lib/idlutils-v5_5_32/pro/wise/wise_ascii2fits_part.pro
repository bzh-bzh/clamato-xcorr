;+
; NAME:
;   wise_ascii2fits_part
;
; PURPOSE:
;   Reformat the WISE catalog release as FITS binary table files.  This
;   procedure only processes one ASCII file.  The procedure wise_ascii2fits
;   calls this procedure to process all the files.  This procedure can also
;   be run in parallel on a batch system.
;
; CALLING SEQUENCE:
;   wise_ascii2fits_part, part=part, [nrowchunk=nrowchunk]
;
; INPUTS:
;   part - (integer) The part (particular ASCII file) of the WISE catalog
;          to process.
;
; OPTIONAL INPUTS:
;   nrowchunk - (long) Number of lines to process at one go.
;
; OUTPUTS:
;   FITS files (see DATA FILES below).
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
;
; REVISION HISTORY:
;   2016-02-12: Written by J. Brownstein, U of Utah.
;-
;
;
;
FUNCTION wisefmt_nextline, ilun
    sline = ''
    READF, ilun, sline
    RETURN, sline
END
;
;
;
PRO wisefmt1, filename, outfile, nrowchunk=nrowchunk, schemafile=schemafile
    ;
    ; Read schema.
    ;
    readcol, schemafile, fname, ffmt, format='(a,a)'
    ntag = N_ELEMENTS(fname)
    IF ~KEYWORD_SET(fname) THEN $
        MESSAGE, 'Invalid data in schema file '+schemafile
    ;
    ; Process schema, creating a structure to hold the data.
    ;
    alldat = 0
    FOR i=0L, ntag-1L DO BEGIN
        IF STRMATCH(ffmt[i],'char*') THEN BEGIN
            thislen = LONG(STREGEX(ffmt[i],'[0-9]+',/EXTRACT))
            thisval = STRING('',FORMAT='(a'+STRTRIM(thislen,2)+')')
        ENDIF ELSE IF STRMATCH(ffmt[i], 'smallint') THEN thisval = 0 $
        ELSE IF STRMATCH(ffmt[i], 'integer') THEN thisval = 0L $
        ELSE IF STRMATCH(ffmt[i], 'int8') THEN thisval = 0LL $
        ELSE IF STRMATCH(ffmt[i], 'serial8') THEN thisval = 0LL $
        ELSE IF STRMATCH(ffmt[i], 'smallfloat') THEN thisval = 0. $
        ELSE IF STRMATCH(ffmt[i], 'decimal*') THEN thisval = 0.D $
        ELSE MESSAGE,'Could not determine the type of '+ffmt[i]
        IF (i EQ 0) THEN alldat = CREATE_STRUCT(fname[i], thisval) $
        ELSE alldat = create_struct(alldat, fname[i], thisval)
    ENDFOR
    ;
    ; We've already established the file exists and that it has a .gz
    ; extension by the time we reach this point.
    ; Count the number of lines in the file, then open the file for
    ; reading one line at a time.
    ;
    shortname = fileandpath(filename[0])
    SPAWN, 'gunzip -c '+filename+'|wc -l', maxlen
    maxlen = LONG(maxlen[0]) > 1L
    OPENR, ilun, filename, ERROR=err, /GET_LUN, /COMPRESS
    IF (err NE 0) THEN BEGIN
        IF KEYWORD_SET(ilun) THEN BEGIN
            CLOSE, ilun
            FREE_LUN, ilun
        ENDIF
        PRINT, 'ERROR opening ' + filename + $
            '! Error message was ' + !ERROR_STATE.MSG
        RETURN
    ENDIF
    ;
    ; Loop over all lines in the file
    ;
    alldat = REPLICATE(alldat, nrowchunk)
    sline = ''
    nullval = -0./0. ; NaN
    nchunk = CEIL((DOUBLE(maxlen)-0.5) / nrowchunk) ; -0.5 to avoid rounding err
    FOR ichunk=0L, nchunk-1L DO BEGIN
        PRINT, 'Parsing file ', filename, ' chunk ', ichunk, ' of ', nchunk
        IF (ichunk EQ nchunk-1L) THEN nthis = maxlen - ichunk*nrowchunk $
        ELSE nthis = nrowchunk
        FOR i=0L, nthis-1L DO BEGIN
            sline = wisefmt_nextline(ilun)
            words = str_sep(sline, '|')
            FOR j=0L, ntag-1L DO BEGIN
                IF SIZE(alldat[i].(j),/TNAME) EQ 'STRING' THEN BEGIN
                    ;
                    ; Force the string to be the same length as the blank
                    ; string in the structure.
                    ;
                    thisfmt = '(A'+STRTRIM(STRLEN(alldat[i].(j)),2)+')'
                    alldat[i].(j) = STRING(words[j]+STRING('',FORMAT=thisfmt), $
                        FORMAT=thisfmt)
                ENDIF ELSE BEGIN
                    ;
                    ; If the length of words[j] is 0, then this should be
                    ; stored as a NULL???  Apparently not, no one cares.
                    ;
                    alldat[i].(j) = words[j]
                ENDELSE
                ;
                ; Replace with NULL where appropriate.
                ;
                IF STRLEN(words[j]) EQ 0 THEN BEGIN
                    CASE SIZE(alldat[i].(j),/TNAME) OF
                        'FLOAT': alldat[i].(j) = nullval
                        'DOUBLE': alldat[i].(j) = nullval
                        ELSE:
                    ENDCASE
                ENDIF
            ENDFOR
        ENDFOR
        mwrfits_chunks, alldat[0L:nthis-1L], outfile, append=(ichunk GT 0), $
            create=(ichunk EQ 0)
    ENDFOR
    RETURN
END
;
;
;
PRO wise_ascii2fits_part, part=part, nrowchunk=nrowchunk
    ;
    ; number of rows to parse in memory at once
    ;
    IF ~KEYWORD_SET(nrowchunk) THEN nrowchunk = 100000L
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
    schemafile = wiseDir + '/wise-cat-schema.txt'
    IF FILE_SEARCH(indir, /TEST_DIRECTORY) EQ '' THEN BEGIN
        PRINT, 'Input directory does not exist: '+indir
        RETURN
    ENDIF
    IF FILE_SEARCH(outdir, /TEST_DIRECTORY) EQ '' THEN BEGIN
        PRINT, 'Output directory does not exist: '+outdir
        RETURN
    ENDIF
    IF ~KEYWORD_SET(part) THEN BEGIN
        PRINT, 'Part number required!'
        RETURN
    ENDIF
    ;
    ; Find the file.
    ;
    infile = 'wise-*-cat-part'+STRING(FORMAT='(i2.2)', part)+'.gz'
    files = FILE_SEARCH(indir, infile, COUNT=nfile)
    IF (nfile EQ 0) THEN BEGIN
        PRINT, 'No input files found!'
        RETURN
    ENDIF
    IF (nfile GT 1) THEN BEGIN
        PRINT, 'Multiple matching part files found!'
        RETURN
    ENDIF
    ;
    ; Process the file.
    ;
    PRINT, 'Working on file: '+files[0]
    fitsfile = repstr(fileandpath(files[0]), '.gz', '.fits')
    radecfile = repstr(fileandpath(files[0]), '.gz', '-radec.fits')
    outfile1 = FILE_SEARCH(outdir, fitsfile, COUNT=ct1)
    outfile2 = FILE_SEARCH(outdir, radecfile, COUNT=ct2)
    IF (ct1 EQ 0) THEN BEGIN
        wisefmt1, files[0], djs_filepath(fitsfile, root_dir=outdir), $
            nrowchunk=nrowchunk, schemafile=schemafile
    ENDIF ELSE PRINT, 'Found an existing copy of catalog file: ' + outfile1[0]
    IF (ct2 EQ 0) THEN BEGIN
        alldat = hogg_mrdfits(djs_filepath(fitsfile, root_dir=outdir), 1, $
            columns=['ra','dec'], nrowchunk=nrowchunk)
        mwrfits, alldat, djs_filepath(radecfile, root_dir=outdir), /create
    ENDIF ELSE PRINT, 'Found an existing copy of radec file: ' + outfile2[0]
    RETURN
END
