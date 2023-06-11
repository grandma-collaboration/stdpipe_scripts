# Generic imports

from matplotlib import pyplot as plt,cm
plt.rc('image', origin='lower', cmap='Blues_r')

import numpy as np
import glob, os
from tqdm import tqdm

from astropy.wcs import WCS
from astropy.io import fits as fits

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time

from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from scipy import stats

import astroscrappy

# Disable some annoying warnings from astropy
import warnings
from astropy.wcs import FITSFixedWarning
warnings.simplefilter(action='ignore', category=FITSFixedWarning)
from astropy.utils.exceptions import AstropyUserWarning
warnings.simplefilter(action='ignore', category=AstropyUserWarning)
from matplotlib import colors
from astropy.coordinates import SkyCoord, search_around_sky

from stdpipe import astrometry, photometry, catalogs, cutouts, templates, subtraction, plots, psf, pipeline, utils

from astropy.visualization.mpl_normalize import ImageNormalize

import reproject
from astropy.table import Table


def spherical_match(ra1, dec1, ra2, dec2, sr=1/3600):
    """Positional match on the sphere for two lists of coordinates.

    Aimed to be a direct replacement for :func:`esutil.htm.HTM.match` method with :code:`maxmatch=0`.

    :param ra1: First set of points RA
    :param dec1: First set of points Dec
    :param ra2: Second set of points RA
    :param dec2: Second set of points Dec
    :param sr: Maximal acceptable pair distance to be considered a match, in degrees
    :returns: Two parallel sets of indices corresponding to matches from first and second lists, along with the pairwise distances in degrees

    """

    idx1,idx2,dist,_ = search_around_sky(SkyCoord(ra1, dec1, unit='deg'), SkyCoord(ra2, dec2, unit='deg'), sr*u.deg)

    dist = dist.deg # convert to degrees

    return idx1, idx2, dist


def calculate_upper(filename,obs,config_data):


	target_ra = config_data["target_ra"]
	target_dec = config_data["target_dec"]
	target_radius = config_data["target_radius"]
	
		
	image = fits.getdata(filename).astype(np.double)
	header = fits.getheader(filename)
	
	# Fix PRISM headers
	if header.get('CTYPE2') == 'DEC---TAN':
		header['CTYPE2'] = 'DEC--TAN'
	for _ in ['CDELTM1', 'CDELTM2', 'XPIXELSZ', 'YPIXELSZ']:
		header.remove(_, ignore_missing=True)
	if header.get('CTYPE1') == 'RA---TAN':
		for _ in ['PV1_1', 'PV1_2']:
		    header.remove(_, ignore_missing=True)
	wcs = WCS(header)
	
	
	#target = SkyCoord(ra = target_ra, dec = target_dec, unit="deg")
	#target_obj = Table({'ra':[target.ra.deg], 'dec':[target.dec.deg]})
    
	#target_obj['x'],target_obj['y'] = wcs.all_world2pix(target_obj['ra'], target_obj['dec'], 0)
		
	#print(image_old.shape[0],target_obj['x'],target_obj['y'])
	#image, header = cutouts.crop_image_centered(image_old, target_obj['x'], target_obj['y'], image_old.shape[0]/4.0, header=header_old)

	verbose=0.0
	

	#header = fits.getheader(filename)
	#image = fits.getdata(filename).astype(np.double)
	
	if verbose==0.0:
		print('Processing %s: %d x %d' % (os.path.split(filename)[1], image.shape[1], image.shape[0]))

	# Image gain, e/ADU
	gain = obs.gain
	print("gain", gain)
	#gain=100.

	# Create mask of bad pixels
	mask = image > 50000 # Rough saturation level

	# Cosmics
	cmask, cimage = astroscrappy.detect_cosmics(image, mask, verbose=False)
	mask |= cmask
	print('Done masking cosmics')

	# Initial WCS
	#f = fits.open(filename)
	wcs = WCS(header)

	# Extract objects
	obj = photometry.get_objects_sextractor(image, mask=mask, r0=1, aper=5.0, wcs=wcs, gain=gain)

	if verbose==0.0:
		print(len(obj), 'objects found')

	ra0,dec0,sr0 = astrometry.get_frame_center(wcs=wcs, width=image.shape[1], height=image.shape[0])
	pixscale = astrometry.get_pixscale(wcs=wcs)

	print("obs.catalogcalib",obs.catalogcalib,obs.filtercalib,obs.filtercalib)
	cat = catalogs.get_cat_vizier(ra0, dec0, sr0, obs.catalogcalib, filters={f'obs.filtercalib':'<22'})
	
	
	if verbose==0.0:
		print(len(cat), 'catalogue stars')


	cat_col_mag = obs.filtercalib
	cat_color_mag1 = obs.filter_ref
	cat_color_mag2 = obs.filtercalib

	# WCS refinement
	wcs = pipeline.refine_astrometry(obj, cat, 5*pixscale, wcs=wcs, order=0, cat_col_mag=cat_col_mag, verbose=True)
	if wcs is None or not wcs.is_celestial:
		print('WCS refinement failed')

	# Update WCS info in the header
	astrometry.clear_wcs(header, remove_comments=True, remove_underscored=True, remove_history=True)
	header.update(wcs.to_header(relax=True))

	# Photometric calibration
	m = pipeline.calibrate_photometry(obj, cat, pixscale=pixscale, cat_col_mag=cat_col_mag, cat_col_mag1=cat_color_mag1, cat_col_mag2=cat_color_mag2, order=0, verbose=True)

	zero_fn = m['zero_fn'] # Function to get the zero point as a function of position on the image

	obj['mag_calib'] = obj['mag'] + zero_fn(obj['x'], obj['y'])


	plots.imshow(image)
	plt.title('Original image')


	# We may roughly estimage the effective gain of the image from background mean and rms as gain = mean/rms**2
	bg,rms = photometry.get_background(image, mask=mask, get_rms=True)

	if verbose==0.0:
		print('Effective gain is %.2f' % np.median(bg/rms**2))

	# We do not have enough stars to study PSF spatial variance, so we use order=0 here
	psf_model,psf_snapshots = psf.run_psfex(image, mask=mask, gain=gain,checkimages=['SNAPSHOTS'], order=0, verbose=True)

	#stop
	plots.imshow(psf_snapshots)
	plt.title('PSF snapshots')
	#plt.show()

	sims = [] # will hold simulated stars


	# We will repeatedly inject the stars, detect objects, and 
	# match them in order to see whether simulated stars are detectable

	for _ in tqdm(range(100)):
		image1 = image.copy()

		# Simulate 20 random stars 
		sim = pipeline.place_random_stars(image1, psf_model, nstars=20, minflux=0.05, maxflux=1e6, wcs=wcs, gain=gain, saturation=50000)
		#stop
		
		sim['mag_calib'] = sim['mag'] + zero_fn(sim['x'], sim['y'])
		sim['detected'] = False
		sim['mag_measured'] = np.nan
		sim['magerr_measured'] = np.nan
		sim['flags_measured'] = np.nan

		mask1 = image1 >= 50000
		
		obj1 = photometry.get_objects_sextractor(image1, mask=mask|mask1, r0=1, aper=5.0, wcs=wcs, gain=gain, minarea=3, sn=0.05)
		obj1['mag_calib'] = obj1['mag'] + zero_fn(obj1['x'], obj1['y'])

		# Positional match within FWHM/2 radius
		oidx,sidx,dist = spherical_match(obj1['ra'], obj1['dec'], sim['ra'], sim['dec'], pixscale*np.median(obj1['fwhm'])/2)
		# Mark matched stars
		sim['detected'][sidx] = True
		# Also store measured magnitude, its error and flags
		sim['mag_measured'][sidx] = obj1['mag_calib'][oidx]
		sim['magerr_measured'][sidx] = obj1['magerr'][oidx]
		sim['flags_measured'][sidx] = obj1['flags'][oidx]
		#print(sim)

		sims.append(sim)

	#     break
		
	from astropy.table import vstack
	sims = vstack(sims)

	plt.clf()
	# Only show unflagged detections
	idx = sims['detected'] & (sims['flags_measured'] == 0)

	plt.errorbar(sims['mag_calib'][idx], (sims['mag_measured'] - sims['mag_calib'])[idx], sims['magerr_measured'][idx], fmt='.', capsize=0, alpha=0.5)

	plt.axhline(0, ls=':', color='black')
	plt.ylim(-1, 1)

	plt.xlabel('Injected magnitude')
	plt.ylabel('Measured - Injected')

	plt.clf()

	h0,b0,_ = plt.hist(sims['mag_calib'], range=[12,22], bins=50, alpha=0.3, label='All simulated stars');
	h1,b1,_ = plt.hist(sims['mag_calib'][sims['detected']], range=[12,22], bins=50, alpha=0.3, label='Detected');
	h2,b2,_ = plt.hist(sims['mag_calib'][idx], range=[12,22], bins=50, alpha=0.3, label='Detected and unflagged');

	plt.legend()

	plt.xlabel('Injected magnitude')

	plt.clf()

	#np.save("test_upper",[sims['mag_calib'], 1/sims['magerr_measured']])
	mag_calib=sims['mag_calib']
	SNR=1/sims['magerr_measured']
	indice=np.where(np.isnan(SNR)==False)
	#print(len(indice))
	mag_calib=mag_calib[indice]
	#print(indice)
	SNR=SNR[indice]
	#print(SNR)
	#plt.plot(mag_calib, SNR, '.', label='Simulated stars')
	#plt.xlim(5, 20)
	#plt.show()

	indice_filter=np.empty(0,dtype=int)
	step=0.2
	mags_bins=np.arange(13,24,step)
	for mag in mags_bins:
		indice_cutmag=np.where((mag_calib >= mag) & (mag_calib < (mag+step)))[0]
		if len(indice_cutmag) > 1:
			SNR_cutmag=SNR[indice_cutmag]
			percentile50=np.percentile(SNR_cutmag, 50)
			percentile25=np.percentile(SNR_cutmag, 25)
			cut_percentile=np.where((SNR_cutmag<percentile50) & (SNR_cutmag>percentile25))[0]
			indice_filter=np.concatenate((indice_filter,indice_cutmag[cut_percentile]),axis=0)
		

	#print("lenindicefilter",indice_filter)		
	mag_calib_cut=mag_calib[indice_filter]
	SNR_cut=SNR[indice_filter]
	#print(mag_calib)
	print("ratio",len(SNR_cut),len(SNR))
	z, res, _, _, _ = np.polyfit(mag_calib_cut, np.log10(SNR_cut), 2,  full = True)
	#print("residus",res)
	p = np.poly1d(z)
	propre_x=np.arange(13,25,0.1)
	indice_value=np.where((pow(10,p(propre_x))>5.0))

	
	PLOT_UPPER=[mag_calib,SNR,'Simulated stars',mag_calib_cut,SNR_cut,obs.filters,'Cut',propre_x,pow(10,p(propre_x)),"Upper lim. "+str(np.round(propre_x[indice_value][-1],1)),5,'Detection limit S/N=5','Measured magnitude','Measured signal to noise ratio']
	

	return np.round(propre_x[indice_value][-1],1),PLOT_UPPER
