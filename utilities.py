"""
This file is for utility routines.
"""
from astropy.table import Table
import os
import numpy as np

data_dir = 'data'

# List of entries to include/update.
new_columns = [
    'ishape_hsm_regauss_derived_sigma_e',
    'ishape_hsm_regauss_derived_rms_e',
    'ishape_hsm_regauss_derived_shape_weight',
    'ishape_hsm_regauss_derived_shear_bias_m',
    'ishape_hsm_regauss_derived_shear_bias_c1',
    'ishape_hsm_regauss_derived_shear_bias_c2'
    ]

def create_calibs(catalog):
    """
    This utility creates calibration factors and other information for HSC shape catalogs, starting
    from a catalog from the LSST Science Pipelines (catalog), and returns a new astropy.Table with
    the additional information.
    """
    # Make output data structure.
    output_catalog = Table()
    for column in new_columns:
        output_catalog[column] = np.zeros(len(catalog))

    output_catalog['ishape_hsm_regauss_derived_sigma_e'] = get_sigma_e_model(catalog)
    output_catalog['ishape_hsm_regauss_derived_rms_e'] = get_erms_model(catalog)
    output_catalog['ishape_hsm_regauss_derived_shape_weight'] = get_weight_model(catalog)

    output_catalog['ishape_hsm_regauss_derived_shear_bias_m'] = get_m_model(catalog)
    model_a = get_a_model(catalog)
    psf_e1, psf_e2 = get_psf_ellip(catalog)
    output_catalog['ishape_hsm_regauss_derived_shear_bias_c1'] = model_a*psf_e1
    output_catalog['ishape_hsm_regauss_derived_shear_bias_c2'] = model_a*psf_e2
    return output_catalog

def get_snr(catalog):
    """
    This utility computes the S/N for each object in the catalog.
    It does not impose any cuts and returns NaNs for invalid S/N values.
    """

    snr = catalog['iflux_cmodel'] / catalog['iflux_cmodel_err']
    return snr

def get_res(catalog):
    """
    Returns the resolution
    """
    return catalog['ishape_hsm_regauss_resolution']

def get_psf_ellip(catalog):
    """
    This utility gets the PSF ellipticity (distortion) from a data or
    sims catalog.
    It does not impose flag cuts, but rather assumes that has already been done
    """
    if 'shape_sdss_psf' in catalog.dtype.names:
        psf_mxx, psf_myy, psf_mxy = catalog['shape_sdss_psf'].T
    else:
        psf_mxx = catalog['ishape_sdss_psf_ixx']
        psf_myy = catalog['ishape_sdss_psf_iyy']
        psf_mxy = catalog['ishape_sdss_psf_ixy']
    return (psf_mxx - psf_myy) / (psf_mxx + psf_myy), 2. * \
        psf_mxy / (psf_mxx + psf_myy)

def get_sigma_e_model(catalog):
    """
    This utility returns a model for the shape measurement uncertainty as a function of SNR and
    resolution.  It uses the catalog directly to get the SNR and resolution values.
    """
    # Get the necessary quantities.
    snr = get_snr(catalog)
    log_snr = np.log10(snr)
    res = get_res(catalog)
    # Fix interpolation for NaN/inf
    m = np.isnan(log_snr) | np.isinf(log_snr)
    log_snr[m] = 1. # set to minimum value that passes nominal cut
    snr[m] = 10.

    # Build the baseline model for sigma_e.
    sigma_e = 0.264 * ((snr/20.)**-0.891) * ((res/0.5)**-1.015)

    # Get the corrections from interpolation amongst saved values.
    data_file = "sigmae_correction_v4-rgc-nr.dat"
    data_file = os.path.join(data_dir, data_file)
    dat = np.loadtxt(data_file).transpose()
    saved_snr = dat[0,:]
    log_saved_snr = np.log10(saved_snr)
    saved_res = dat[1,:]
    saved_corr = dat[2,:]

    # Interpolate the corrections (which multiply the power-law results).
    result = grid_interpolate(log_saved_snr, saved_res, saved_corr, log_snr, res)
    return result*sigma_e

def get_erms_model(catalog):
    """
    This utility returns a model for the RMS ellipticity as a function of SNR and
    resolution.  It uses the catalog directly to get the SNR and resolution values.
    """
    # Get the necessary quantities.
    snr = get_snr(catalog)
    log_snr = np.log10(snr)
    res = get_res(catalog)
    # Fix interpolation for NaN/inf
    m = np.isnan(log_snr) | np.isinf(log_snr)
    log_snr[m] = 1. # set to minimum value that passes nominal cut
    snr[m] = 10.

    # Get saved model values.
    data_file = "erms_v4-rgc-nr.dat"
    data_file = os.path.join(data_dir, data_file)
    dat = np.loadtxt(data_file).transpose()
    saved_snr = dat[0,:]
    log_saved_snr = np.log10(saved_snr)
    saved_res = dat[1,:]
    saved_model = dat[2,:]

    # Interpolate the e_rms values and return them.
    result = grid_interpolate(log_saved_snr, saved_res, saved_model, log_snr, res)
    return result

def get_weight_model(catalog):
    """
    This utility returns a model for the shape measurement weight as a
    function of SNR and resolution.  It relies on two other routines
    to get models for the intrinsic shape RMS and measurement error.
    """
    sigmae_meas = get_sigma_e_model(catalog)
    erms = get_erms_model(catalog)
    return 1./(sigmae_meas**2 + erms**2)

def get_m_model(catalog, weight_bias=True):
    """
    Routine to get a model for calibration bias m given some input snr
    and resolution values or arrays.  By default, includes the weight bias contribution, which is
    appropriate if you plan to use the weights provided with the calibrations.
    """
    # Get the necessary quantities.
    snr = get_snr(catalog)
    log_snr = np.log10(snr)
    res = get_res(catalog)
    # Fix interpolation for NaN/inf
    m = np.isnan(log_snr) | np.isinf(log_snr)
    log_snr[m] = 1. # set to minimum value that passes nominal cut
    snr[m] = 10.

    m = -0.1408*((snr/20.)**-1.23)*((res/0.5)**1.76) - 0.0214
    data_file = "model_m_a_correction_v4-rgc-nr.dat"
    data_file = os.path.join(data_dir, data_file)

    dat = np.loadtxt(data_file).transpose()
    saved_snr = dat[0,:]
    log_saved_snr = np.log10(saved_snr)
    saved_res = dat[1,:]
    saved_m = dat[2,:]

    # Interpolate the model residuals and return them.  These are additional
    # biases beyond the power-law model and so should be added to the power-law.
    result = grid_interpolate(log_saved_snr, saved_res, saved_m, log_snr, res)

    if weight_bias:
        result += get_msel_model(catalog)

    return result + m

def get_msel_model(catalog):
    """
    Routine to get a model for calibration bias m due to weight bias, given some input snr
    and resolution values or arrays.
    """
    # Get the necessary quantities.
    snr = get_snr(catalog)
    log_snr = np.log10(snr)
    res = get_res(catalog)
    # Fix interpolation for NaN/inf
    m = np.isnan(log_snr) | np.isinf(log_snr)
    log_snr[m] = 1. # set to minimum value that passes nominal cut
    snr[m] = 10.

    m = -1.31 + (27.26 + (snr/20.)**-1.22) / (res + 20.8)
    data_file = "model_msel_asel_correction_v4-rgc-nr.dat"
    data_file = os.path.join(data_dir, data_file)

    dat = np.loadtxt(data_file).transpose()
    saved_snr = dat[0,:]
    log_saved_snr = np.log10(saved_snr)
    saved_res = dat[1,:]
    saved_m = dat[2,:]

    # Interpolate the model residuals and return them.  These are additional
    # biases beyond the power-law model and so should be added to the power-law.
    result = grid_interpolate(log_saved_snr, saved_res, saved_m, log_snr, res)

    return result + m

def get_a_model(catalog, weight_bias=True):
    """
    Routine to get a model for additive bias coefficient a given some input snr
    and resolution values or arrays.
    """
    # Get the necessary quantities.
    snr = get_snr(catalog)
    log_snr = np.log10(snr)
    res = get_res(catalog)
    # Fix interpolation for NaN/inf
    m = np.isnan(log_snr) | np.isinf(log_snr)
    log_snr[m] = 1. # set to minimum value that passes nominal cut
    snr[m] = 10.

    a = 0.175 * ((snr/20.)**-1.07) * (res - 0.508)
    data_file = "model_m_a_correction_v4-rgc-nr.dat"
    data_file = os.path.join(data_dir, data_file)
    dat = np.loadtxt(data_file).transpose()
    saved_snr = dat[0,:]
    log_saved_snr = np.log10(saved_snr)
    saved_res = dat[1,:]
    saved_a = dat[4,:]

    # Interpolate the model residuals and return them.  These are additional
    # biases beyond the power-law model and so should be added to the power-law.
    result = grid_interpolate(log_saved_snr, saved_res, saved_a, log_snr, res)
    if weight_bias:
        result += get_asel_model(catalog)
    return result + a

def get_asel_model(catalog):
    """
    Routine to get a model for additive bias coefficient `a` due to weight bias given some input snr
    and resolution values or arrays.
    """
    # Get the necessary quantities.
    snr = get_snr(catalog)
    log_snr = np.log10(snr)
    res = get_res(catalog)
    # Fix interpolation for NaN/inf
    m = np.isnan(log_snr) | np.isinf(log_snr)
    log_snr[m] = 1. # set to minimum value that passes nominal cut
    snr[m] = 10.

    a = -0.089 * (res-0.71) * ((snr/20.)**-2.2)
    data_file = "model_msel_asel_correction_v4-rgc-nr.dat"

    data_file = os.path.join(data_dir, data_file)
    dat = np.loadtxt(data_file).transpose()
    saved_snr = dat[0,:]
    log_saved_snr = np.log10(saved_snr)
    saved_res = dat[1,:]
    saved_a = dat[4,:]

    # Interpolate the model residuals and return them.  These are additional
    # biases beyond the power-law model and so should be added to the power-law.
    result = grid_interpolate(log_saved_snr, saved_res, saved_a, log_snr, res)
    return result + a

def grid_interpolate(x, y, z, eval_x, eval_y):
    """
    This is a utility for interpolating a 2D function z(x, y) linearly to values
    (x, y) = (eval_x, eval_y), but also enabling extrapolation beyond the (x, y)
    bounds using the nearest neighbor method.
    """
    from scipy.interpolate import griddata
    result = griddata((x, y), z, (eval_x, eval_y), method='linear')
    nn_result = griddata((x, y), z, (eval_x, eval_y), method='nearest')
    mask = np.isnan(result)
    result[mask] = nn_result[mask]
    return result
