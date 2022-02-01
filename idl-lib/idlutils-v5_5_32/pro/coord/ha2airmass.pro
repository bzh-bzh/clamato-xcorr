;+
; NAME:
;   ha2airmass
;
; PURPOSE:
;   Compute airmass.
;
; CALLING SEQUENCE:
;   airmass = ha2airmass( dec, [ ha=, $
;    longitude=, latitude=, altitude= ] )
;
; INPUTS:
;   dec            - Declination [degrees]
;   ha             - Hour angle (degrees); default to 0
;
; OPTIONAL KEYWORDS:
;   longitude      - Longitude of observatory;
;                    default to (360-105.820417) deg for APO
;   latitute       - Latitude of observatory; default to 32.780361 deg for APO
;   altitude       - Altitude of observatory; default to 2788 m for APO
;
; OUTPUTS:
;   airmass        - Airmass; 1.0 for zenith
;
; OPTIONAL OUTPUTS:
;   ipa            - Position angle for image rotator (degrees)
;
; COMMENTS:
;   This routine only returns sec(z) for the airmass.
;   Formula from Smart, Spherical Astronomy.
;
; EXAMPLES:
;
; BUGS:
;   EQUINOX does nothing except for the IPA calculation!
;   ALTITUDE is unused!
;
; PROCEDURES CALLED:
;   ct2lst
;   ll2uv()
;   precess
;
; REVISION HISTORY:
;   10-May-2000  Written by D. Schlegel, Princeton, & D. Hogg, IAS
;   02-Jun-2000  Fixed minor bugs, Schlegel
;   05-Nov-2000  Added HA keyword
;-
;------------------------------------------------------------------------------
function tai2air_crossprod, aa, bb

  if (n_elements(A) GT 3) or (n_elements(B) GT 3) then begin
     print, 'Sorry - only 3-vectors, one at a time'
     return, 0
  endif

  A = reform(aa, 3) & B=reform(bb, 3)
  C = (A-A)*1.                       ; zero vector of type float or double
  C[0] = determ([[1, 0, 0], [A], [B]])
  C[1] = determ([[0, 1, 0], [A], [B]])
  C[2] = determ([[0, 0, 1], [A], [B]])

  return, C
end
;------------------------------------------------------------------------------
; angle (Degrees) between any vectors A, B (not necessarily unit vectors)

function tai2air_ang, a, b
  c = (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
  c = ((c/sqrt(total(a^2)*total(b^2))) < 1) > (-1)
  theta = acos(c)*180./!dpi

  return, theta
end
;------------------------------------------------------------------------------
function ha2airmass, dec, equinox1, ha=hadeg, $
 longitude=longitude, latitude=latitude, altitude=altitude

   ; Default to location of Apache Point Observatory
   if (NOT keyword_set(longitude)) then longitude = 360. - 105.820417d0
   if (NOT keyword_set(latitude)) then latitude = 32.780361d0
   if (NOT keyword_set(altitude)) then altitude = 2788.

   DRADEG = 180.d0 / !DPI

   decrad = dec / DRADEG
   if (keyword_set(hadeg)) then harad = hadeg / DRADEG $
    else harad = 0.d0
   latrad = latitude / DRADEG

   ;----------
   ; Compute airmass with spherical trig

   cosz = sin(decrad)*sin(latrad)+ cos(decrad)*cos(harad)*cos(latrad)
   airmass = 1.d0 / cosz

   return, airmass
end
;------------------------------------------------------------------------------
