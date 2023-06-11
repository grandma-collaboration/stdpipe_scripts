# Generic imports

import os
import glob
import json
import time as ts
import Calculate_upper as upper
import GetObsinfov2 as getobs
import humanize
import socket
from stdpipe import astrometry, photometry, catalogs, cutouts, templates, subtraction, plots, psf, pipeline, utils
from matplotlib import colors
from astropy.utils.exceptions import AstropyUserWarning
from astropy.wcs import FITSFixedWarning
import warnings
from astropy.table import Table
import astroscrappy
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits as fits
from astropy.wcs import WCS
import numpy as np
from matplotlib import pyplot as plt, cm
plt.rc('image', origin='lower', cmap='Blues_r')


# Disable some annoying warnings from astropy
warnings.simplefilter(action='ignore', category=FITSFixedWarning)
warnings.simplefilter(action='ignore', category=AstropyUserWarning)


# import PIL

# ==========================================================================
# ==========================================================================
# Def
# ==========================================================================
# ==========================================================================


def check_withinregion(target_ra, target_dec, target_radius, src_ra, src_rec):

    dist = astrometry.spherical_distance(
        src_ra, src_rec, target_ra, target_dec)
    if dist < target_radius:
        return 1.0
    return 0.0


def create_science_results(repo_photo):

    os.system("rm -r " + repo_photo)
    os.system("mkdir -p " + os.path.join(repo_photo, "Cross-Match"))
    os.system("mkdir -p " + os.path.join(repo_photo, "Soustraction"))
    return repo_photo


def detection_source(filename, repo_photo, detection, obs, verbose, mask_simage_threshold, show_preprocimg, show_masksimg, cat_ref, stars_threshold, filter_ref, cat_col_mag_color, cat_col_mag1_color, cat_col_mag2_color, show_calib, show_candidate, template_type, template_ref, show_photorefimage, cat_colsous_mag, show_diffimg, sextractor_default, target_ra, target_dec, target_radius, config_data):

    candidates = ""

    image = fits.getdata(filename).astype(np.double)
    header = fits.getheader(filename)

    instrument = obs.instru

    time = utils.get_obs_time(header, verbose=False)
    utc = header.get("OBSDATE")
    mjd = obs.mjd
    # exp = header.get('EXP')
    exp = obs.exp
    fname = obs.filters

    try:
        gain = float(header.get('GAIN'))
    except TypeError:
        print("no gain value found in the header")
        gain = obs.gain
    if verbose is "True":
        print('Processing %s: filter %s gain %.2f at %s' %
              (filename, fname, gain, time))

    # As the image is already dark-subtracted and flatfielded, we need just to mask saturated stars and cosmic rays
    # The mask is just a binary frame with the same size as the image where True means that this pixel should not be used for the analysis
    mask = image > mask_simage_threshold
    # mask = mask * True
    # print(mask)

    cmask, cimage = astroscrappy.detect_cosmics(image, mask, verbose=True)
    # print(cmask)
    if verbose is True:
        print('Done masking cosmics: %d pixels masked' % np.sum(cmask))
    mask |= cmask

    if show_preprocimg == "True":
        norm = ImageNormalize(vmin=3550, vmax=3700, stretch=SqrtStretch())
        plots.imshow(image)
        plt.title('Pre-processed image')
        plt.imshow(image, norm=norm)
        plt.savefig(os.path.join(repo_photo, "Pre-processed image.png"))
        plt.clf()

    if show_masksimg == "True":
        plots.imshow(mask)
        plt.title('Mask')
        plt.savefig(os.path.join(repo_photo, "Mask.png"))
        plt.clf()

    # print("mask",mask)
    mask = mask*True
    # We will detect objects using SExtractor and get their measurements in apertures with 3 pixels radius
    obj = photometry.get_objects_sextractor(
        image, mask=mask, aper=3.0, gain=gain, edge=10)
    if verbose is True:
        print(len(obj), 'objects detected')

    # Rough estimation of average FWHM of detected objects, taking into account only unflagged (e.g. not saturated) ones
    fwhm = np.median(obj['fwhm'][obj['flags'] == 0])
    if verbose is True:
        print('Average FWHM is %.1f pixels' % fwhm)

    # We will pass this FWHM to measurement function so that aperture and background radii will be relative to it.
    # We will also reject all objects with measured S/N < 5
    obj = photometry.measure_objects(
        obj, image, mask=mask, fwhm=fwhm, gain=gain, aper=1.0, bkgann=[5, 7], sn=3, verbose=True)
    if verbose is True:
        print(len(obj), 'objects properly measured')

    # Load initial WCS
    wcs = WCS(header)

    # Get the center position, size and pixel scale for the image
    center_ra, center_dec, center_sr = astrometry.get_frame_center(
        wcs=wcs, width=image.shape[1], height=image.shape[0])
    pixscale = astrometry.get_pixscale(wcs=wcs)
    if verbose is True:
        print('Frame center is %.2f %.2f radius %.2f deg, %.2f arcsec/pixel' %
              (center_ra, center_dec, center_sr, pixscale*3600))

    # Let's get cat objects brighter than r=19 mag
    print("cat_ref", cat_ref, "filter_ref", filter_ref)
    cat = catalogs.get_cat_vizier(
        center_ra, center_dec, center_sr, cat_ref, filters={filter_ref: '<'+"19"})
    if verbose is True:
        print(len(cat), 'catalogue stars')

    # Let's use SCAMP for astrometric refinement.
    wcs = pipeline.refine_astrometry(
        obj, cat, 5*pixscale, wcs=wcs, method='scamp', cat_col_mag=filter_ref, verbose=True)

    if wcs is not None:
        # Update WCS info in the header
        astrometry.clear_wcs(header, remove_comments=True,
                             remove_underscored=True, remove_history=True)
        header.update(wcs.to_header(relax=True))

    # Photometric calibration using 2 arcsec matching radius, r magnitude, g-r color and second order spatial variations
    print("cat_col_mag_color", cat_col_mag_color, "cat_col_mag1_color",
          cat_col_mag1_color, "cat_col_mag2_color", cat_col_mag2_color)
    # ts.sleep(250)
    m = pipeline.calibrate_photometry(obj, cat, sr=3/3600, cat_col_mag=cat_col_mag_color, cat_col_mag1=cat_col_mag1_color,
                                      cat_col_mag2=cat_col_mag2_color, max_intrinsic_rms=0.02, order=2, verbose=True)
    if verbose is True:
        print(m['zero_fn'])
    # ts.sleep(120)
    # The code above automatically augments the object list with calibrated magnitudes, but we may also do it manually
    obj['mag_calib'] = obj['mag'] + m['zero_fn'](obj['x'], obj['y'])
    obj['mag_calib_err'] = np.hypot(
        obj['magerr'], m['zero_fn'](obj['x'], obj['y'], get_err=True))

    # There is a handy plotting function useful for quick checking of fitting quality and uncorrected trends. It marks the stars used for final model fit with red dots, and flagged (e.g. saturated) stars rejected from the fitting from the start - with yellow diagonal crosses.

    # Photometric residuals as a function of catalogue magnitude
    if show_calib == True:
        plt.subplot(211)
        plots.plot_photometric_match(m)
        plt.ylim(-0.5, 0.5)

        # Photometric residuals as a function of catalogue color
        plt.subplot(212)
        plots.plot_photometric_match(m, mode='color')
        plt.ylim(-0.5, 0.5)
        plt.xlim(0.0, 1.5)
        plt.savefig(os.path.join(repo_photo, "Photometric residuals.png"))
        plt.clf()
        # plt.show()

        # Zero point (difference between catalogue and instrumental magnitudes for every star) map
        plt.subplot(121)
        plots.plot_photometric_match(
            m, mode='zero', bins=6, show_dots=True, aspect='equal')
        plt.title('Zero point')
        # plt.show()

        # Fitted zero point model with second-order spatial polynomial term
        plt.subplot(122)
        plots.plot_photometric_match(
            m, mode='model', bins=6, show_dots=True, aspect='equal')
        plt.title('Zero point model')
        plt.savefig(os.path.join(repo_photo, "Zeropoint.png"))
        plt.clf()
        plt.close()

    candidateswithin = []
    if detection == "cat":
        # Filtering of transient candidates
        candidates = pipeline.filter_transient_candidates(
            obj, cat=cat, sr=5/3600, vizier=['vsx'], verbose=True)
        plt.clf()
        for i, cand in enumerate(candidates):
            if check_withinregion(target_ra, target_dec, target_radius, cand['ra'], cand['dec']) == 1.0:
                print('Candidate %d with mag = %.2f +/- %.2f at x/y = %.1f %.1d and RA/Dec = %.6f %.6f' %
                      (i, cand['mag_calib'], cand['mag_calib_err'], cand['x'], cand['y'], cand['ra'], cand['dec']))
                cutout = cutouts.get_cutout(
                    image, cand, 20, mask=mask, header=header)
                cutout['template'] = templates.get_hips_image(
                    template_ref, header=cutout['header'])[0]
                candidateswithin.append([cand['mag_calib'], cand['mag_calib_err'],
                                        cand['ra'], cand['dec'], cand['mag_calib'], cand['mag_calib_err']])
                plots.plot_cutout(cutout, qq=[0.5, 99.9], stretch='linear')
                plt.savefig(os.path.join(repo_photo, "Candidate-"+str(
                    np.round(cand['ra'], 4))+"_"+str(np.round(cand['dec'], 4))+"WITHINTARGET.png"))
                plt.clf()
            else:
                cutout = cutouts.get_cutout(
                    image, cand, 20, mask=mask, header=header)
                cutout['template'] = templates.get_hips_image(
                    'flux-Rp/I/350/gaiaedr3', header=cutout['header'])[0]
                candidateswithin.append([cand['mag_calib'], cand['mag_calib_err'],
                                        cand['ra'], cand['dec'], cand['mag_calib'], cand['mag_calib_err']])
                plots.plot_cutout(cutout, qq=[0.5, 99.9], stretch='linear')
                plt.savefig(os.path.join(repo_photo, "Candidate-"+str(
                    np.round(cand['ra'], 4))+"_"+str(np.round(cand['dec'], 4))+"NOTINTARGET.png"))
            # plt.show()

    if detection == "soustraction":
        tmpl = ""

        # Get r band image from PanSTARRS with science image original resolution and orientation
        # tmpl = templates.get_hips_image('PanSTARRS/DR1/r', wcs=wcs, width=image.shape[1], height=image.shape[0], get_header=False)

        if template_type == "cat":
            # Get "image from template_ref with 2x resolution (as we have large pixel scale) and downsample it
            print('flux-Rp/I/350/gaiaedr3', template_ref, "image_shape1",
                  image.shape[1]*2, "image_shape2", image.shape[0]*2)
            img_ref_tmp = repo_photo+"/tmp.fits"
            print(img_ref_tmp)
            tmpl = templates.get_hips_image(template_ref, wcs=wcs[::0.5, ::0.5], width=image.shape[1]
                                            * 2, height=image.shape[0]*2, get_header=False, local_file_path=img_ref_tmp, verbose=True)
            # print(tmpl)
            # np.save(tmpl)
            tmpl = utils.rebin_image(tmpl, 2)

        if template_type == "obs":
            # Let's load the image and parse its header!
            tmpl = fits.getdata(template_ref).astype(np.double)

        if config_data["show_refimagesous"] == "True":
            plots.imshow(tmpl)
            plt.savefig(os.path.join(repo_photo, "Template_image.png"))
            plt.clf()

        # To demonstrate it, let's perform the photometry on template and match it with catalogue
        tobj = photometry.get_objects_sextractor(
            tmpl, mask=~np.isfinite(tmpl), wcs=wcs, aper=3, gain=1.0)
        print(cat, cat_col_mag_color)
        m_tmpl = pipeline.calibrate_photometry(
            tobj, cat, sr=1/3600, cat_col_mag=cat_col_mag_color, max_intrinsic_rms=0.1, verbose=True)  # cat_col_mag_err='e_RPmag'

        # pipeline.calibrate_photometry(tobj, cat, sr=1/3600, cat_col_mag=cat_colsous_mag, max_intrinsic_rms=0.1, verbose=False) #cat_col_mag_err='e_rmag'

        if config_data["show_photorefimage"] == "True":
            plt.subplot(211)
            plots.plot_photometric_match(m_tmpl, mode='mag')
            plt.text(14, 1.2, 'Unflagged saturated stars!',
                     color='red', fontsize=20)
            plt.text(15.2, 0.4, 'Scatter due to undersampled HiPS template',
                     color='blue', fontsize=12)
            plt.ylim(-2, 2)
            plt.savefig(os.path.join(
                repo_photo, "Photometry_referenceimage.png"))
            plt.clf()

        tmask = templates.mask_template(
            tmpl, cat, cat_col_mag=cat_col_mag_color, cat_saturation_mag=15, wcs=wcs, dilate=3, verbose=True)

        if show_masksimg == "True":
            plots.imshow(tmask)
            plt.title('Template mask')
            plt.savefig(os.path.join(repo_photo, "Mask.png"))
            plt.clf()

        # Run the subtraction getting back all possible image planes, assuming the template to be noise-less, and estimating image noise model from its statistics.

        import photutils

        bg = photutils.Background2D(
            image, 128, mask=mask, exclude_percentile=30).background
        tbg = photutils.Background2D(
            tmpl, 128, mask=tmask, exclude_percentile=30).background
        res = subtraction.run_hotpants(image-bg, tmpl-tbg, mask=mask, template_mask=tmask, get_convolved=True, get_scaled=True,
                                       get_noise=True, image_fwhm=fwhm, template_fwhm=1.5, image_gain=gain, template_gain=1e6, err=True, verbose=True)
        diff, conv, sdiff, ediff = res

        sdiff[tmask] = 0.0

        if show_diffimg == "True":
            plots.imshow(sdiff, vmin=-3, vmax=10)
            plt.title('Noise-scaled difference image')
            plt.savefig(os.path.join(repo_photo, "Difference_image.png"))
            plt.clf()

        dmask = diff == 1e-30  # Bad pixels as marked by HOTPANTS

        # Get PSF model and store it to temporary file
        psf_model = psf.run_psfex(
            image, mask=mask, order=0, gain=gain, psffile='/tmp/psf.psf', verbose=True)

        # Run SExtractor on difference image with custom noise model, returning object footprints and some additional fields
        sobj, segm = photometry.get_objects_sextractor(diff, mask=mask | tmask | dmask, err=ediff, edge=10, wcs=wcs, aper=2.0, extra_params=['CLASS_STAR', 'NUMBER'], extra={
                                                       'SEEING_FWHM': fwhm, 'STARNNW_NAME': sextractor_default}, checkimages=['SEGMENTATION'], psf='/tmp/psf.psf', verbose=True)

        # Perform forced aperture photometry, again with custom noise model and forced zero background level
        sobj = photometry.measure_objects(sobj, diff, mask=mask | tmask | cmask | dmask, fwhm=fwhm, aper=1.0, bkgann=[
                                          5, 7], sn=3, verbose=True, bg=0, err=ediff)

        # The difference is in original image normalization, so we know photometric zero point
        sobj['mag_calib'] = sobj['mag'] + m['zero_fn'](sobj['x'], sobj['y'])
        sobj['mag_calib_err'] = np.hypot(
            sobj['magerr'], m['zero_fn'](sobj['x'], sobj['y'], get_err=True))

        # We may immediately reject flagged objects as they correspond to imaging artefacts (masked regions)
        sobj = sobj[sobj['flags'] == 0]

        print(len(sobj), 'transient candidates found in difference image line 307')

        candidates = pipeline.filter_transient_candidates(
            sobj, sr=2/3600, flagged=True, vizier=['ps1', 'vsx'], skybot=True, time=time, verbose=True)

        withinXRT = np.zeros(len(candidates))
        for i, cand in enumerate(candidates):
            print('Candidate %d with mag = %.2f +/- %.2f at x/y = %.1f %.1d and RA/Dec = %.4f %.4f' %
                  (i, cand['mag_calib'], cand['mag_calib_err'], cand['x'], cand['y'], cand['ra'], cand['dec']))

        for i, cand in enumerate(candidates):

            if show_candidate == "True":

                print('SPREAD_MODEL = %.3f +/- %.3f, CLASS_STAR = %.2f' %
                      (cand['spread_model'], cand['spreaderr_model'], cand['CLASS_STAR']))

                if check_withinregion(target_ra, target_dec, target_radius, cand['ra'], cand['dec']) == 1.0:
                    cutout = cutouts.get_cutout(image, cand, 20, mask=mask | tmask | dmask, diff=diff, template=tmpl, convolved=conv, err=ediff, footprint=(
                        segm == cand['NUMBER']), header=header, filename=filename, time=time)
                    plots.plot_cutout(cutout, ['image', 'template', 'convolved', 'diff', 'footprint', 'mask'], qq=[
                                      0.5, 99.5], stretch='linear')
                    print([cand['mag_calib'], cand['mag_calib_err'], cand['ra'],
                          target_obj['dec'], cand['mag_calib'], cand['mag_calib_err']])
                    candidateswithin.append([cand['mag_calib'], cand['mag_calib_err'], cand['ra'],
                                            target_obj['dec'], cand['mag_calib'], cand['mag_calib_err']])
                    plt.savefig(os.path.join(repo_photo, "Candidate-"+str(
                        np.round(cand['ra'], 5))+"_"+str(np.round(cand['dec'], 5))+"WITHINTARGET.png"))
                    plt.clf()
                else:
                    cutout = cutouts.get_cutout(
                        image, cand, 20, mask=mask, header=header)
                    cutout['template'] = templates.get_hips_image(
                        template_ref, header=cutout['header'])[0]
                    candidateswithin.append([cand['mag_calib'], cand['mag_calib_err'],
                                            cand['ra'], cand['dec'], cand['mag_calib'], cand['mag_calib_err']])
                    plots.plot_cutout(cutout, qq=[0.5, 99.9], stretch='linear')
                    plt.savefig(os.path.join(repo_photo, "Candidate-"+str(
                        np.round(cand['ra'], 5))+"_"+str(np.round(cand['dec'], 5))+"NOTINTARGET.png"))
                    plt.clf()

    # print("candidatewithin",candidateswithin)
    # candidateswithin=np.array(candidateswithin)
    # np.save(os.path.join(repo_photo, "Candidates_"+str(detection)), candidateswithin)
    return obs, candidateswithin


def compare_candidates(repo_photo, fits_file, xrt_ra, xrt_dec, xrt_radius):
    detection = "cat"
    candidates_sous = np.load(os.path.join(
        repo_photo, "Cross-Match", "Candidates_"+str(detection)+".npy"))

    detection = "cat"
    candidates_cat = np.load(os.path.join(
        repo_photo, "Cross-Match", "Candidates_"+str(detection)+".npy"))

    candidates_cross = []
    for i, candcat in enumerate(candidates_cat):
        for j, candsous in enumerate(candidates_sous):
            val = check_withinregion(
                candsous['ra'], candsous['dec'], 1/3600., candcat['ra'], candcat['dec'])
            if val == 1:
                val_xrt = check_withinregion(
                    xrt_ra, xrt_dec, xrt_radius, candcat['ra'], candcat['dec'])
                if val_xrt == 1:
                    print('Candidate with CATALOGUE mag = %.2f +/- %.2f , RA/Dec = %.4f %.4f' %
                          (candcat['mag_calib'], candcat['mag_calib_err'], candcat['ra'], candcat['dec']))
                    print("and")
                    print('Candidate with SOUS mag = %.2f +/- %.2f , RA/Dec = %.4f %.4f' %
                          (candsous['mag_calib'], candsous['mag_calib_err'], candsous['ra'], candsous['dec']))
                    candidates_cross.append([candcat['mag_calib'], candcat['mag_calib_err'], candcat['ra'],
                                            candcat['dec'], candsous['mag_calib'], candsous['mag_calib_err']])
    return candidates_cross


def create_logs(repo_photo, fits_file, obs, candidates_cross, upperlimit, config):
    name = (fits_file.split("/")[-1]).split("_")[0]+"_"+obs.instru + \
        "_"+obs.filters+"_"+str(np.round(obs.mjd, 3))+"_photometry.csv"
    print("repo", repo_photo+"/"+name)
    author = config["author"]
    GRB_time = Time(config["event_time"]).mjd
    with open(repo_photo+"/"+name, 'a+') as output:
        header = "#utc,mjd,instrument,filter,mag,magerr,limiting_mag,exp,magsys,catalogcalib,filtercalib,"+"\n"
        output.write(header)
        if len(candidates_cross) > 0:
            for i in np.arange(len(candidates_cross)):
                diff_time = (obs.mjd-GRB_time)  # in seconds
                sentence = str(obs.obserdate)+","+str(obs.mjd)+","+str(obs.instru)+","+str(obs.filters)+","+str(float(np.round(candidates_cross[i][0], 2)))+","+str(float(round(float(candidates_cross[i][1]), 2)))+","+str(np.round(upperlimit, 1))+","+str(obs.exp)+","+str(obs.magsys)+","+str(obs.catalogcalib)+","+str(obs.filtercalib)+", #"+"RA: "+str(
                    np.round(candidates_cross[i][2], 5))+" dec: "+str(np.round(candidates_cross[i][3], 5))+" ,STDPIPE mag "+str(np.round(candidates_cross[i][4], 2))+","+str(np.round(candidates_cross[i][5], 2))+" , "+str(author)+" , "+str(np.round(diff_time, 3))+" delay, " + str(humanize.precisedelta(diff_time*24.0*3600.0, minimum_unit="hours"))+"\n"
                print(str(float(round(float(candidates_cross[i][1]), 2))))
                output.write(sentence)
        else:
            diff_time = (obs.mjd-GRB_time)  # in seconds
            sentence = str(obs.obserdate)+","+str(obs.mjd)+","+str(obs.instru)+","+str(obs.filters)+","+str(0.0)+","+str(0.0)+","+str(np.round(upperlimit, 1))+","+str(obs.exp)+","+str(obs.magsys)+","+str(obs.catalogcalib)+","+str(obs.filtercalib) + \
                ", #"+"RA: "+str(0.0)+" dec: "+str(0.0)+" ,STDPIPE mag "+str(0.0)+","+str(0.0)+" , "+str(author)+" , "+str(
                    np.round(diff_time, 3))+" delay, " + str(humanize.precisedelta(diff_time*24.0*3600.0, minimum_unit="hours"))+"\n"
            print(sentence)
            output.write(sentence)


def create_html(repo_photo):
    with open(os.path.join(repo_photo, "index_results.html"), "w") as fhtml:
        title = os.path.split(repo_photo)[1]
        html = ""
        html += "<!doctype html>\n"
        html += "<html lang=\"fr\">\n"
        html += "<head>\n"
        html += "  <meta charset=\"utf-8\">\n"
        html += f"  <title>{title}</title>\n"
        html += "</head>\n"
        html += "<body>\n"
        html += f"<h1>{title}</h1>\n"
        # --
        detecs = ["Soustraction", "Cross-Match"]
        for detec in detecs:
            html += f"<h2>{detec}</h2>\n"
            pngfiles = glob.glob(os.path.join(
                repo_photo, detec, "*WITHINTARGET.png"))
            for pngfile in pngfiles:
                # im = PIL.Image.open('whatever.png')
                # width, height = im.size
                fname = detec + "/" + os.path.basename(pngfile)
                html += f"<img src=\"{fname}\" alt=\"{os.path.basename(pngfile)}\"><br>\n"
        html += "<body>\n"
        fhtml.write(html)


def analyse_image(fits_file, target_ra, target_dec, target_radius, GRB_time, repo_res, images_folder, verbose, results_logs, force, config_data):
    # ==========================================================================
    # Common
    # ==========================================================================

    ################### PLOTS ###################
    # Preprocessed scientific images
    show_preprocimg = config_data["show_processing"]

    # Mask applied in the science image
    show_masksimg = config_data["show_masksimg"]

    # Show plots for calibration (with zero points)
    show_calib = config_data["show_calib"]

    # Show plots candidates with catalogue cross matched
    show_candidate = config_data["show_candidate"]

    # Show plots for reference image in the soustraction process
    show_refimagesous = config_data["show_refimagesous"]

    # Show refs
    show_photorefimage = config_data["show_photorefimage"]

    # Show image differencing
    show_diffimg = config_data["show_diffimg"]

    # Show upperlimit calculation
    show_upper = config_data["show_upper"]

    fits_file_comp = os.path.join(repo_res, images_folder, fits_file)
    root_file, extension = os.path.splitext(fits_file)
    root_file = root_file.split("/")[-1]
    repo_photo = repo_res+results_logs+root_file

    obs = getobs.treat_GRBdata(fits_file_comp)
    print("Parameter obs", obs)

    ################### SCIENCE IMAGE TREATEMENT ###################

    # To remove mask saturated stars and cosmic rays in the science image
    mask_simage_threshold = 5000

    # Choose the reference catalog to extract the photometry of stars in the field
    cat_ref = obs.catalogcalib

    # Select stars only brighter than a threshold and a filter
    stars_threshold = 19.0
    filter_ref = obs.filtercalib

    # Add color term
    cat_col_mag_color = obs.filtercalib
    cat_col_mag1_color = obs.filter_ref
    cat_col_mag2_color = obs.filtercalib

    ################### IMAGE SOUSTRACTION ###################

    template_type = "cat"

    if template_type == "cat":
        template_ref = obs.catalogref  # "PanSTARRS/DR1/r"
    if template_type == "obs":
        template_cat_ref = ""

    cat_colsous_mag = obs.filtercalib
    # Create the image analyzing repository ################################"
    create_science_results(repo_photo)

    ################# " CALCULATE UPPER ###################

    upperlimit = -100.

    if show_upper == "True":
        print("########## UPPER ################")
        upperlimit, PLOT_UPPER = upper.calculate_upper(
            fits_file_comp, obs, config_data)
        plt.clf()
        plt.plot(PLOT_UPPER[0], PLOT_UPPER[1], '.', label=PLOT_UPPER[2])
        plt.plot(PLOT_UPPER[3], PLOT_UPPER[4], 'r.',
                 PLOT_UPPER[5], label=PLOT_UPPER[6])
        plt.plot(PLOT_UPPER[7], PLOT_UPPER[8], '-', label=PLOT_UPPER[9])

        plt.yscale('log')
        plt.axhline(PLOT_UPPER[10], ls=':', color='red', label=PLOT_UPPER[11])

        plt.xlabel('Measured magnitude')
        plt.ylabel('Measured signal to noise ratio')

        plt.title(fits_file.split("/")[-1])

        plt.legend()

        plt.xlim(14, 23)
        # plt.show()
        # print(repo_photo)
        plt.savefig(os.path.join(repo_photo, obs.instru+"_" +
                    obs.filters+"_"+str(np.round(obs.mjd, 3))+"_upperlimit.png"))
        plt.show()
        plt.clf()

    candidates_cross = []
    if force == "True":
        print("########## FORCE ################")
        candidates_cross = force_photometry(target_ra, target_dec, target_radius, fits_file_comp, obs,
                                            mask_simage_threshold, cat_ref, cat_col_mag_color, cat_col_mag1_color, cat_col_mag2_color, filter_ref)

    # CROSS-MATCH ################################"
    else:
        if config_data["cross_match"] == "True":
            print("########## CROSS ################ =" *
                  20 + " Cross-Match " + "="*20)
            detection = "cat"
            repo_cat = os.path.join(repo_photo, "Cross-Match")
            res = detection_source(fits_file_comp, repo_cat, detection, obs, verbose, mask_simage_threshold, show_preprocimg, show_masksimg, cat_ref, stars_threshold, filter_ref, cat_col_mag_color, cat_col_mag1_color,
                                   cat_col_mag2_color, show_calib, show_candidate, template_type, template_ref, show_photorefimage, cat_colsous_mag, show_diffimg, config_data["sextractor"], target_ra, target_dec, target_radius, config_data)
            obs, candidates_cross = res
            # candidates_cross = compare_candidates(repo_photo,fits_file,target_ra,target_dec,target_radius)

    # SOUSTRACTION ################################"
    if config_data["soustraction"] == "True":
        print("="*20 + " Soustraction " + "="*20)
        detection = "soustraction"
        repo_sous = os.path.join(repo_photo, "Soustraction")
        res = detection_source(fits_file_comp, repo_sous, detection, obs, verbose, mask_simage_threshold, show_preprocimg, show_masksimg, cat_ref, stars_threshold, filter_ref, cat_col_mag_color, cat_col_mag1_color,
                               cat_col_mag2_color, show_calib, show_candidate, template_type, template_ref, show_photorefimage, cat_colsous_mag, show_diffimg, config_data["sextractor"], target_ra, target_dec, target_radius, config_data)

        utc, mjd, instrument, filter_name, candidates_sous, exp, magsys, catalogcalib, filtercalib = res

    """
    ################ COMPARE ################################"
    if verbose is True:
    print( "="*20 + " Compare " + "="*20)
    #exp="0.0"
    candidates_cross = compare_candidates(repo_photo,fits_file,target_ra,target_dec,target_radius)
    """
    # Create LOGS & HTML ################################"
    # print(candidates_cross)
    create_logs(repo_photo, fits_file, obs,
                candidates_cross, upperlimit, config_data)
    create_html(repo_photo)
    if verbose is True:
        print("="*20 + " Finished " + "="*20)


def force_photometry(target_ra, target_dec, target_radius, filename, obs, mask_simage_threshold, cat_ref, cat_col_mag_color, cat_col_mag1_color, cat_col_mag2_color, filter_ref):

    # Target
    candidates_cross = []

    image_old = fits.getdata(filename).astype(np.double)
    header_old = fits.getheader(filename)

    if header_old.get('CTYPE2') == 'DEC---TAN':
        header_old['CTYPE2'] = 'DEC--TAN'
    for _ in ['CDELTM1', 'CDELTM2', 'XPIXELSZ', 'YPIXELSZ']:
        header_old.remove(_, ignore_missing=True)
    if header_old.get('CTYPE1') == 'RA---TAN':
        for _ in ['PV1_1', 'PV1_2']:
            header_old.remove(_, ignore_missing=True)

    wcs = WCS(header_old)
    target = SkyCoord(ra=target_ra, dec=target_dec, unit="deg")
    target_obj = Table({'ra': [target.ra.deg], 'dec': [target.dec.deg]})

    target_obj['x'], target_obj['y'] = wcs.all_world2pix(
        target_obj['ra'], target_obj['dec'], 0)

    # print(image_old.shape[0],target_obj['x'],target_obj['y'])
    image, header = cutouts.crop_image_centered(
        image_old, target_obj['x'], target_obj['y'], image_old.shape[0]/3.0, header=header_old)

    # print(image.shape[0])
    # stop
    # print(header)
    wcs = WCS(header)
    time = utils.get_obs_time(header, verbose=False)
    fname = obs.filters
    gain = float(obs.gain)
    mask = image > mask_simage_threshold
    cmask, cimage = astroscrappy.detect_cosmics(image, mask, verbose=True)
    print('Done masking cosmics: %d pixels masked' % np.sum(cmask))
    mask |= cmask
    mask = mask*False
    obj = photometry.get_objects_sextractor(
        image, mask=mask, aper=5.0, gain=gain, edge=10)
    fwhm = np.median(obj['fwhm'][obj['flags'] == 0])
    obj = photometry.measure_objects(
        obj, image, mask=mask, fwhm=fwhm, gain=gain, aper=1.0, bkgann=[5, 7], sn=3, verbose=True)
    print(len(obj), 'objects properly measured')

    wcs = WCS(header)
    center_ra, center_dec, center_sr = astrometry.get_frame_center(
        wcs=wcs, width=image.shape[1], height=image.shape[0])
    pixscale = astrometry.get_pixscale(wcs=wcs)
    print('Frame center is %.2f %.2f radius %.2f deg, %.2f arcsec/pixel' %
          (center_ra, center_dec, center_sr, pixscale*3600))

    cat = catalogs.get_cat_vizier(
        center_ra, center_dec, center_sr, cat_ref, filters={filter_ref: '<'+"21"})
    print(len(cat), 'catalogue stars')
    wcs = pipeline.refine_astrometry(
        obj, cat, 5*pixscale, wcs=wcs, method='scamp', cat_col_mag='rmag', verbose=True)
    if wcs is not None:
        # Update WCS info in the header
        astrometry.clear_wcs(header, remove_comments=True,
                             remove_underscored=True, remove_history=True)
        header.update(wcs.to_header(relax=True))
    # stop
    m = pipeline.calibrate_photometry(obj, cat, sr=1/3600, cat_col_mag=cat_col_mag_color, cat_col_mag1=cat_col_mag1_color,
                                      cat_col_mag2=cat_col_mag2_color, max_intrinsic_rms=0.02, order=2, verbose=True)

    wcs = WCS(header)

    target = SkyCoord(ra=target_ra, dec=target_dec, unit="deg")
    target_obj = Table({'ra': [target.ra.deg], 'dec': [target.dec.deg]})

    target_obj['x'], target_obj['y'] = wcs.all_world2pix(
        target_obj['ra'], target_obj['dec'], 0)

    # Forced photometry
    target_obj = photometry.measure_objects(target_obj, image, mask=mask, fwhm=fwhm, aper=1.0, bkgann=[
                                            5, 7], sn=None, verbose=verbose, gain=gain)
    print(target_obj)
    target_obs = np.array(target_obj)

    # Assign calibrated magnitude
    target_obj['mag_calib'] = target_obj['mag'] + \
        m['zero_fn'](target_obj['x'], target_obj['y'], target_obj['mag'])
    target_obj['mag_calib_err'] = np.hypot(target_obj['magerr'], m['zero_fn'](
        target_obj['x'], target_obj['y'], target_obj['mag'], get_err=True))

    # 5-sigma limit
    target_obj['mag_limit'] = -2.5*np.log10(5*target_obj['fluxerr']) + m['zero_fn'](
        target_obj['x'], target_obj['y'], target_obj['mag'])

    candidates_cross.append([target_obj['mag_calib'], target_obj['mag_calib_err'], target_obj['ra'],
                            target_obj['dec'], target_obj['mag_calib'], target_obj['mag_calib_err']])
    # print(candidates_cross)
    # plt.title('GRB2210')
    # plots.imshow(images)
    # plt.show()
    # plt.clf()
    # stop
    return candidates_cross
# ==========================================================================
# ==========================================================================
# Main
# ==========================================================================
# ==========================================================================


def load_config(filename):
    with open(filename) as config_file:
        config_data = json.load(config_file)
    return config_data


def print_configuration(hostname, author, sextractor_default, res_folder, images_folder, fits_file, res_logs):
    print("=" * 20 + " Configuration " + "=" * 20)
    print(f"hostname = {hostname}")
    print(f"author = {author}")
    print(f"sextractor_default = {sextractor_default}")
    print(f"repo_res = {res_folder}")
    print(f"images_folder = {images_folder}")
    print(f"fits_file = {fits_file}")
    print(f"results_logs = {res_logs}")


def process_image(fits_file, target_ra, target_dec, target_radius, event_time, res_folder, images_folder, verbose, res_logs, force, config):

    analyse_image(fits_file, target_ra, target_dec, target_radius, event_time,
                  res_folder, images_folder, verbose, res_logs, force, config)
    try:
        # analyse_image(fits_file, target_ra, target_dec, target_radius, GRB_time, repo_res, images_folder, verbose, results_logs, force)
        print("NO TRY PASS")
    except:
        pass


def main():
    config_data = load_config('config.json')

    target_ra = config_data["target_ra"]
    target_dec = config_data["target_dec"]
    target_radius = config_data["target_radius"]
    event_time = Time(config_data["event_time"]).mjd

    res_folder = config_data["res_folder"]
    images_folder = config_data["images_folder"]
    res_logs = config_data["res_logs"]
    select_img = config_data["select_img"]
    verbose = config_data["verbose"]
    if verbose == "True":
        verbose = True
    else:
        verbose = False

    force = config_data["force"]

    hostname = socket.gethostname()

    author = config_data["author"]
    upperlimit = 0.0
    sextractor_default = config_data["sextractor"]

    fits_list = glob.glob(images_folder + "*.fit*", recursive=True)

    for fits_file in fits_list:
        if select_img in fits_file:
            if verbose:
                print_configuration(hostname, author, sextractor_default,
                                    res_folder, images_folder, fits_file, res_logs)
            process_image(fits_file, target_ra, target_dec, target_radius, event_time,
                          res_folder, images_folder, verbose, res_logs, force, config_data)


if __name__ == '__main__':
    main()
