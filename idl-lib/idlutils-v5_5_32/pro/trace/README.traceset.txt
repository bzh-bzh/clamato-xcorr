
DEFINITION

traceset (or "trace set"):

A set of coefficient vectors defining a set of functions
over a common independent-variable domain specified by "xmin"
and "xmax" values. The functions in the set are defined in
terms of a linear combination of basis functions (such as
Legendre of Chebyshev polynonials) up to a specified maximum
order, weighted by the values in the coefficient vectors, and
evaluated using a suitable affine rescaling of the
dependent-variable domain (such as [xmin, xmax] -> [-1, 1]
for Legendre polynomials). For evaluation purposes, the
domain is assumed by default to be a zero-based integer
baseline from xmin to xmax such as would correspond to a
digital detector pixel grid.

Etymology: from the original motivating use case of
encoding the "traces" of multiple spectra across the
detector in a multi-fiber spectrograph.

(Definition from A. Bolton, U. of Utah, June 2015)

