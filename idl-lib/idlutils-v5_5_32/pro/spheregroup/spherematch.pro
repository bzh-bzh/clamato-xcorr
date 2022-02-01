;+
; NAME:
;   spherematch
;
; PURPOSE:
;   Take two sets of ra/dec coords and efficiently match them. It
;   returns all matches between the sets of objects, and the distances
;   between the objects. The matches are returned sorted by increasing
;   distance. A parameter "maxmatch" can be set to the number of matches
;   returned for each object in either list. Thus, maxmatch=1 (the default)
;   returns the closest possible set of matches. maxmatch=0 means to
;   return all matches
;
; CALLING SEQUENCE:
;   spherematch, ra1, dec1, ra2, dec2, matchlength, match1, match2, $
;                distance12, [maxmatch=maxmatch]
;
; INPUTS:
;   ra1         - ra coordinates in degrees (N-dimensional array)
;   dec1        - dec coordinates in degrees (N-dimensional array)
;   ra2         - ra coordinates in degrees (N-dimensional array)
;   dec2        - dec coordinates in degrees (N-dimensional array)
;   matchlength - distance which defines a match (degrees)
;
; OPTIONAL INPUTS:
;   maxmatch    - Return only maxmatch matches for each object, at
;                 most. Defaults to maxmatch=1 (only the closest
;                 match for each object). maxmatch=0 returns all
;                 matches.
;   estnmatch   - Estimate of the TOTAL number of matches.  If this is
;                 absent or wrong, the C code is called twice,
;                 doubling execution time!
;
; OPTIONAL KEYWORDS:
;   /verbose    - Be verbose about warnings
;
; OUTPUTS:
;   match1     - List of indices of matches in list 1; -1 if no matches
;   match2     - List of indices of matches in list 2; -1 if no matches
;   distance12 - Distance of matches; 0 if no matches
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The code breaks the survey region into chunks of size
;   4*matchlength. Matches are then performed by considering only
;   objects in a given chunk and neighboring chunks. This makes the
;   code fast.
;
;   The matches are returned sorted by distance.
;
;   If you have a big list and a small list, call with the
;   BIG LIST FIRST!!!
;   i.e.
;
;   spherematch, BIGra, BIGdec, SMALLra, SMALLdec, matchlength, $
;                      matchBIG, matchSMALL, distance12
;
;   This method is inherently asymmetric.  Calling in this order will
;   exploit the asymmetry to reduce memory usage and execution time.
;
; EXAMPLES:
;
; BUGS:
;   Behavior at poles not well tested.
;
; PROCEDURES CALLED:
;   gcirc
;   idlutils_so_ext()
;   Dynamic link to spherematch.c
;
; REVISION HISTORY:
;   20-Jul-2001  Written by Mike Blanton, Fermiland
;   01-Mar-2006  estnmatch keyword added - D. Finkbeiner, Princeton
;          estnmatch allows the caller to estimate the number of
;          matches, so the wrapper can allocate memory for results before
;          calling the C code.  If the estimate is absent or wrong,
;          the code is called a second time (as before).
;   2014-09-09  Avoid calling C code if both sets of coordinates only contain
;          one object.
;   2015-10-21  Be quiet about warnings regarding efficiency if not 
;          explicitly told to be /verbose 
;-
;------------------------------------------------------------------------------
PRO spherematch, ra1, dec1, ra2, dec2, matchlength, match1, match2, $
                 distance12, maxmatch=maxmatch, chunksize=chunksize, $
                 estnmatch=estnmatch, verbose=verbose
    ;
    ; Set default return values
    ;
    match1 = -1L
    match2 = -1L
    distance12 = 0L
    ;
    ; Need at least 3 parameters
    ;
    IF (N_PARAMS() LT 7) THEN $
        MESSAGE, 'Syntax - spherematch, ra1, dec1, ra2, dec2, matchlength, ' + $
            'match1, match2, distance12, [maxmatch=]'
    IF (N_ELEMENTS(maxmatch) EQ 0) THEN maxmatch=1L ELSE $
        IF (maxmatch LT 0L) THEN MESSAGE, 'Illegal maxmatch value: '+maxmatch
    IF ~KEYWORD_SET(chunksize) THEN chunksize=MAX([4.*matchlength,0.1])
    ;
    ; Check array sizes.
    ;
    npoints1 = N_ELEMENTS(ra1)
    IF (npoints1 LE 0L) THEN $
        MESSAGE, 'Need array with > 0 elements.'
    IF (npoints1 NE N_ELEMENTS(dec1)) THEN $
        MESSAGE, 'ra1 and dec1 must have same length.'
    npoints2 = N_ELEMENTS(ra2)
    IF (npoints2 LE 0L) THEN $
        MESSAGE, 'Need array with > 0 elements.'
    IF (npoints2 NE N_ELEMENTS(dec2)) THEN $
        MESSAGE, 'ra2 and dec2 must have same length.'
    IF (matchlength LE 0L) THEN $
        MESSAGE, 'Need matchlength > 0'
    IF (npoints2 GT npoints1 and KEYWORD_SET(verbose) ne 0) THEN $
        MESSAGE, 'spherematch works best when ra1, dec1 contain more points than ra2, dec2.', $
            /INFORMATIONAL
    ;
    ; Check coordinate ranges.
    ;
    ibadra1 = WHERE(ra1 LT 0. OR ra1 GT 360., nbadra1)
    IF (nbadra1 GT 0) THEN $
        MESSAGE, 'spherematch does not accept RA outside 0 to 360 (RA1).'
    ibadra2 = WHERE(ra2 LT 0. OR ra2 GT 360., nbadra2)
    IF (nbadra2 GT 0) THEN $
        MESSAGE, 'spherematch does not accept RA outside 0 to 360 (RA2).'
    ibaddec1 = WHERE(dec1 LT -90. OR dec1 GT 90., nbaddec1)
    IF (nbaddec1 GT 0) THEN $
        MESSAGE, 'spherematch does not accept DEC outside -90 to 90 (DEC1).'
    ibaddec2 = WHERE(dec2 LT -90. OR dec2 GT 90., nbaddec2)
    IF (nbaddec2 GT 0) THEN $
        MESSAGE, 'spherematch does not accept DEC outside -90 to 90 (DEC2).'
    ;
    ; Check for degenerate case.
    ;
    IF (npoints1 EQ npoints2 AND npoints2 EQ 1L) THEN BEGIN
        gcirc, 2, ra1, dec1, ra2, dec2, as12
        odistance12 = as12/3.6D3
        IF (odistance12[0] LT matchlength) THEN BEGIN
            match1 = [0L]
            match2 = [0L]
            IF ARG_PRESENT(distance12) THEN distance12=[odistance12]
        ENDIF
        RETURN
    ENDIF
    ;
    ; Continue with usual processing.
    ;
    soname = FILEPATH('libspheregroup.'+idlutils_so_ext(), $
                        root_dir=GETENV('IDLUTILS_DIR'), SUBDIRECTORY='lib')
    onmatch = KEYWORD_SET(estnmatch) ? LONG(estnmatch) : 0L
    onmatch_save = onmatch
    ;
    ; First pass on matching C code.
    ;
    omatch1     = LONARR(onmatch > 1)
    omatch2     = LONARR(onmatch > 1)
    odistance12 = DBLARR(onmatch > 1)
    retval = CALL_EXTERNAL(soname, 'spherematch', $
                            LONG(npoints1), DOUBLE(ra1), DOUBLE(dec1), $
                            LONG(npoints2), DOUBLE(ra2), DOUBLE(dec2), $
                            DOUBLE(matchlength), DOUBLE(chunksize), $
                            LONG(omatch1), LONG(omatch2), $
                            DOUBLE(odistance12), LONG(onmatch))
    IF onmatch EQ 0 THEN RETURN
    ;
    ; Second pass, if we did not allocate enough space before
    ;
    IF onmatch GT onmatch_save THEN BEGIN
        omatch1     = LONARR(onmatch)
        omatch2     = LONARR(onmatch)
        odistance12 = DBLARR(onmatch)
        retval = CALL_EXTERNAL(soname, 'spherematch', $
                                LONG(npoints1), DOUBLE(ra1), DOUBLE(dec1), $
                                LONG(npoints2), DOUBLE(ra2), DOUBLE(dec2), $
                                DOUBLE(matchlength), DOUBLE(chunksize), $
                                LONG(omatch1), LONG(omatch2), $
                                DOUBLE(odistance12), LONG(onmatch))
    ENDIF
    ;
    ; trim padding in output arrays
    ;
    IF onmatch LT N_ELEMENTS(omatch1) THEN BEGIN
        omatch1     = omatch1[0:onmatch-1]
        omatch2     = omatch2[0:onmatch-1]
        odistance12 = odistance12[0:onmatch-1]
    ENDIF
    ;
    ; Retain only desired matches
    ;
    sorted = SORT(odistance12)
    IF (maxmatch GT 0L) THEN BEGIN
        gotten1 = LONARR(npoints1)
        gotten2 = LONARR(npoints2)
        nmatch=0L
        FOR i = 0L, onmatch-1L DO BEGIN
            IF ((gotten1[omatch1[sorted[i]]] LT maxmatch) AND $
                (gotten2[omatch2[sorted[i]]] LT maxmatch)) THEN BEGIN
                gotten1[omatch1[sorted[i]]] = gotten1[omatch1[sorted[i]]]+1L
                gotten2[omatch2[sorted[i]]] = gotten2[omatch2[sorted[i]]]+1L
                nmatch = nmatch+1L
            ENDIF
        ENDFOR
        gotten1    = LONARR(npoints1)
        gotten2    = LONARR(npoints2)
        match1     = LONARR(nmatch)
        match2     = LONARR(nmatch)
        distance12 = DBLARR(nmatch)
        nmatch=0L
        FOR i = 0L, onmatch-1L DO BEGIN
            IF ((gotten1[omatch1[sorted[i]]] LT maxmatch) AND $
                (gotten2[omatch2[sorted[i]]] LT maxmatch)) THEN BEGIN
                gotten1[omatch1[sorted[i]]] = gotten1[omatch1[sorted[i]]]+1L
                gotten2[omatch2[sorted[i]]] = gotten2[omatch2[sorted[i]]]+1L
                match1[nmatch] = omatch1[sorted[i]]
                match2[nmatch] = omatch2[sorted[i]]
                distance12[nmatch] = odistance12[sorted[i]]
                nmatch = nmatch+1L
            ENDIF
        ENDFOR
    ENDIF ELSE BEGIN
        nmatch=onmatch[sorted]
        onmatch=0
        match1=omatch1[sorted]
        omatch1=0
        match2=omatch2[sorted]
        omatch2=0
        IF ARG_PRESENT(distance12) THEN distance12=odistance12[sorted]
        odistance12=0
    ENDELSE
    ;
    ;
    ;
    RETURN
END
