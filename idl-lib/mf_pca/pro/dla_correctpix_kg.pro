;+
; NAME:
;   dla_correctpix_kg
;
; PURPOSE:
;   For an input wavelength vector, return a vector
;   of exp(tau) to allow corrections for the damping wings of a
;   DLA. The core of the DLA should already be masked.
;
; CALLING SEQUENCE:
;   dla_correctpix_kg, wave_in, n_hi, [z_abs=z_abs, /log_column]
;
; INPUTS:
;   wave_in   - Wavelength vector, either in observed frame or DLA
;               restframe. Z_ABS needs to be input for the former.
;   n_hi      - Neutral hydrogen column density
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;   z_abs       - Redshift of DLA. Input only if input wavelength is in
;                 observed frame. 
;   log_column  - Set if input N_HI is in log10
;
; OUTPUTS:
;   exptauvec   - Vector corresponding to input WAVE_IN vector, with the
;                 value of exp(tau)
;
;
; INTERNAL SUPPORT ROUTINES:
;   tau_lorentz
;
;
; MODIFICATION HISTORY:
;
;  KG Lee   23/3/2011   - Ensured that single-precision number is
;                         returned 
;-

FUNCTION TAU_LORENTZ, deltawave, n_hi
; Return optical depth at some rest wavelength separation DELTAWAVE
; (in angstroms) for a column density of N_HI

e = 4.80320d-10
m_e = 9.109389d-28
c=2.99792d10
gamma_lambda = 7616.d
lambda_cm = 1.21567d-5 
f_lu = 0.4164d

; Make sure than n_hi is a scalar
if (size(n_hi))[0] GT 0 then n_hi = (n_hi)[0] 

n_hi = replicate(n_hi, n_elements(deltawave))

emc_factor = e^2/m_e/c^3

tau = emc_factor/(4.*!pi) * f_lu * gamma_lambda * lambda_cm $
  * n_hi * (1215.67)^2 / deltawave^2 

RETURN, float(TAU)
END 

FUNCTION DLA_CORRECTPIX_KG, wave_in, n_hi, z_abs=z_abs, $
                            log_column=log_column

if keyword_set(z_abs) then wave = wave_in/(1.+z_abs) $
  else wave = wave_in

if keyword_set(log_column) then nhi = 10.^n_hi $
  else nhi = n_hi

dwave = wave - 1215.67

exptau = exp(tau_lorentz(dwave, nhi))

RETURN, exptau
END 
