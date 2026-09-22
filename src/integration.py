"""Native-bin number and lung-deposited surface area integration."""
import numpy as np


def delta_log10(diameter_nm):
    """Study convention: forward log10 spacing; repeat the last spacing."""
    d = np.asarray(diameter_nm, dtype=float)
    if d.ndim != 1 or len(d) < 2 or np.any(d <= 0) or np.any(np.diff(d) <= 0):
        raise ValueError('At least two positive, strictly increasing native diameters are required.')
    widths = np.diff(np.log10(d))
    return np.append(widths, widths[-1])


def integrate_ldsa(values, diameter_nm, deposition_fraction, lower_nm=17., upper_nm=500., representation='delta_n'):
    """Return um2 cm-3. Missing channels are omitted; wholly missing rows stay NaN.

    Select native channel centres within the inclusive target interval. No
    extrapolation, interpolation, missing-bin renormalization or boundary-bin
    fractioning is applied. Pass only the channels of one native configuration.
    """
    x = np.asarray(values, dtype=float)
    d = np.asarray(diameter_nm, dtype=float)
    f = np.asarray(deposition_fraction, dtype=float)
    if x.ndim != 2 or d.ndim != 1 or f.shape != d.shape or x.shape[1] != len(d):
        raise ValueError('Input must be observations x channels with aligned diameter and DF vectors.')
    if np.any(~np.isfinite(d)) or np.any(d <= 0) or np.any(np.diff(d) <= 0):
        raise ValueError('Diameters must be finite, positive and strictly increasing.')
    if np.any(~np.isfinite(f)) or np.any((f < 0) | (f > 1)):
        raise ValueError('Deposition fractions must be finite and within [0, 1].')
    if np.any(x[np.isfinite(x)] < 0) or np.any(np.isinf(x)):
        raise ValueError('Apply source quality control before integration.')
    if representation == 'dndlog10dp':
        x = x * delta_log10(d)
    elif representation != 'delta_n':
        raise ValueError('representation must be delta_n or dndlog10dp')
    sel = (d >= lower_nm) & (d <= upper_nm)
    if not sel.any():
        return np.full(len(x), np.nan)
    y = x[:, sel] * np.pi * (d[sel] / 1000.)**2 * f[sel]
    result = np.nansum(y, axis=1)
    result[np.isnan(y).all(axis=1)] = np.nan
    return result
