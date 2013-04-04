#!/usr/bin/env python
#
# Created 12-Mar-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
#
# usage: 

# Import libraries
import math
import os
import random
from itertools import cycle

import numpy as np
import pyfits as pf
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as p
	
import yanny

BSR = '/clusterfs/riemann/raid006/bosswork/boss/spectro/redux'

redVersion = 'v5_6_0'
blueVersion = 'test/dmargala/redux/v5_6_0'

colors = ['k','b','g','r','c','y','m']
linestyles = ['-','--','-.',':']

def plotImage(image, filename="",logScale=True,colorScaleLimit=4,\
	xlabel="", ylabel="", title="",dpi=200,aspect=4,cmap='seismic'):
	if not filename:
		print "plotImage: ouput filename not specifed!"
		return 0
	# Set up figure
	p.ioff()
	p.clf()
	fig = p.figure()
	ax = fig.add_subplot(111)
	# Draw the figure
	if logScale:
		logScaleLimit = np.log10(colorScaleLimit)
		cax = ax.imshow(np.log10(image), cmap=cmap, aspect=aspect, vmin=-logScaleLimit, vmax=logScaleLimit)
		nTicks = np.round(colorScaleLimit)
		tickValues = np.concatenate((sorted(1./(np.arange(nTicks)+1)),np.arange(nTicks-1)+2))
		tickLabels = ["%.2f"%tickValue for tickValue in tickValues]
		cbar = fig.colorbar(cax, ticks=np.log10(tickValues))
		cbar.ax.set_yticklabels(tickLabels)
	else:
		cax = ax.imshow(image, cmap=cmap, aspect=aspect, vmin=0, vmax=colorScaleLimit)
		cbar = fig.colorbar(cax)
	# Add labels
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_title(title)
	# Save the figure
	fig.savefig(filename, dpi=dpi)
	return fig

def getExposureIDs(plate,mjd,fiber='',camera='',verbose=False):
	plan2dFilename = os.path.join(BSR,redVersion,plate,'spPlan2d-%s-%s.par'%(plate,mjd))
	plan2d = yanny.yanny(plan2dFilename)
	spexp = plan2d['SPEXP']
	exposureIDs = []
	for i in range(len(spexp['flavor'])):
		if spexp['flavor'][i] == 'science':
			if camera == '' and fiber == '':
				for c in range(4):
					exposureIDs.append('-'.join(spexp['name'][i][c].split('-')[1:3]).split('.')[0])
			else:
				exposureIDs.append('-'.join(spexp['name'][i][0].split('-')[1:3]).split('.')[0])
	return exposureIDs

def plotspFluxcalibRatios(plate,mjd,fiberlist=[],verbose=False):
	# Look up exposure ids
	if verbose:
		print 'Looking up exposureIDs for %s-%s...' % (plate,mjd)
	exposureIDs = getExposureIDs(plate,mjd,camera='b')
	# Iterate over exposure ids
	if verbose:
		print 'Found %d frames...' % len(exposureIDs)
	for exposureID in exposureIDs:
		bossFluxcalibname = os.path.join(BSR,redVersion,plate,'spFluxcalib-%s.fits.gz'%exposureID)
		blueFluxcalibname = os.path.join(BSR,blueVersion,plate,'spFluxcalib-%s.fits.gz'%exposureID)
		# Read boss standard calibrated spFluxcalib
		bosshdu = pf.open(bossFluxcalibname)
		npix = bosshdu[0].header['naxis1']
		wave = np.arange(npix)
		bossflux = bosshdu[0].data
		bossheader = bosshdu[0].header
		bosshdu.close()
		# Read blue standard calibrated spFluxcalib
		bluehdu = pf.open(blueFluxcalibname)
		bluenpix = bluehdu[0].header['naxis1']
		blueflux = bluehdu[0].data
		bluehdu.close()
		if bluenpix != npix:
			print '!! warning bossnpix != bluenpix: %d != %d' % (npix,bluenpix)
			print '!! skipping %s-%s...' % (plate,mjd)
			return -1
		# Plot ratio
		savedir = os.path.join(BSR,blueVersion,plate)
		# imageFilename = os.path.join(savedir,'spFluxcalibRatioImage-%s.png'%exposureID)
		# plotImage(bossflux/blueflux,filename=imageFilename,\
		# 	xlabel='Pixel',ylabel='Fiberid',\
		# 	title='spFluxcalib-%s Ratio (%s-%s)\nDefault Pipeline / Pipeline Using Blue Standards'%(exposureID,plate,mjd),\
		# 	aspect=6)

		# # Plot spectra
		# p.ioff()
		# p.clf()
		# fig = p.figure(1)
		# maxFiberid = 500
		# randomFibers = random.sample(set(range(maxFiberid)),maxFiberid/10)
		# for i in randomFibers if not fiberlist else (np.array(fiberlist) % 500):
		# 	color = colors[i%len(colors)]
		# 	alpha = .5
		# 	lw = 1
		# 	# Draw the spectra with the appropriate line properties
		# 	lines = p.plot(wave,bossflux[i]/blueflux[i])
		# 	p.setp(lines,color=color,alpha=alpha,lw=lw)
		# # Save the figure
		# p.xlabel('Pixels')
		# p.ylabel('Default Pipeline / Pipeline Using Blue Standards')
		# p.title('spFluxcalib-%s Ratio (%s-%s)'%(exposureID,plate,mjd))
		# p.axis([wave[0],wave[-1],0,2])
		# p.grid(True)
		# savename = os.path.join(savedir,'spFluxcalibRatio-%s.png'%exposureID)
		# fig.savefig(savename, dpi=200)
		
def plotspFluxcorrRatios(plate,mjd,fiberlist=[],verbose=False):
	# Look up exposure ids
	exposureIDs = getExposureIDs(plate,mjd)
	# Iterate over exposure ids
	if verbose:
		print 'Found %d frames...' % len(exposureIDs)
	for exposureID in exposureIDs:
		bossFluxcorrname = os.path.join(BSR,redVersion,plate,'spFluxcorr-%s.fits.gz'%exposureID)
		blueFluxcorrname = os.path.join(BSR,blueVersion,plate,'spFluxcorr-%s.fits.gz'%exposureID)
		# Read boss standard calibrated spFluxcalib
		bosshdu = pf.open(bossFluxcorrname)
		npix = bosshdu[0].header['naxis1']
		wave = np.arange(npix)
		bossflux = bosshdu[0].data
		bosshdu.close()
		# Read blue standard calibrated spFluxcalib
		bluehdu = pf.open(blueFluxcorrname)
		bluenpix = bluehdu[0].header['naxis1']
		blueflux = bluehdu[0].data
		bluehdu.close()
		if bluenpix != npix:
			print '!! warning bossnpix != bluenpix: %d != %d' % (npix,bluenpix)
			print '!! skipping %s-%s...' % (plate,mjd)
			return -1
		# Plot Ratio
		savedir = os.path.join(BSR,blueVersion,plate)
		imageFilename = os.path.join(savedir,'spFluxcorrRatioImage-%s.png'%exposureID)
		plotImage(bossflux/blueflux,filename=imageFilename,\
			xlabel='Pixel',ylabel='Fiberid',\
			title='spFluxcorr-%s Ratio (%s-%s)\nDefault Pipeline / Pipeline Using Blue Standards'%(exposureID,plate,mjd),\
			aspect=6)
		# # Plot spectra
		# p.ioff()
		# p.clf()
		# fig = p.figure(1)
		# maxFiberid = 500
		# randomFibers = random.sample(set(range(maxFiberid)),maxFiberid/10)
		# for i in randomFibers if not fiberlist else (np.array(fiberlist) % 500):
		# 	color = colors[i%len(colors)]
		# 	alpha = 1
		# 	lw = 1
		# 	# Draw the spectra with the appropriate line properties
		# 	lines = p.plot(wave,bossflux[i]/blueflux[i])
		# 	p.setp(lines,color=color,alpha=alpha,lw=lw)
		# # Save the figure
		# p.xlabel('Pixels')
		# p.ylabel('Default Pipeline / Pipeline Using Blue Standards')
		# p.title('spFluxcorr-%s Ratio (%s-%s)'%(exposureID,plate,mjd))
		# p.axis([wave[0],wave[-1],0,2])
		# p.grid(True)
		# savename = os.path.join(savedir,'spFluxcorrRatio-%s.png'%exposureID)
		# fig.savefig(savename, dpi=200)

def plotspFluxdistortRatios(plate,mjd,fiberlist=[],verbose=False):
	# Look up exposure ids
	if verbose:
		print 'Loading spFluxdistort file for %s-%s...' % (plate,mjd)
	bossFluxdistortFilename = os.path.join(BSR,redVersion,plate,'spFluxdistort-%s-%s.fits'%(plate,mjd))
	blueFluxdistortFilename = os.path.join(BSR,blueVersion,plate,'spFluxdistort-%s-%s.fits'%(plate,mjd))
	# Read boss standard calibrated spFluxdistort
	bosshdu = pf.open(bossFluxdistortFilename)
	bossdistort = bosshdu[0].data
	bosshdu.close()
	# Read blue standard calibrated spFluxdistort
	bluehdu = pf.open(blueFluxdistortFilename)
	bluedistort = bluehdu[0].data
	bluehdu.close()

	npix = len(bossdistort[0])
	pixels = range(npix)
	if len(bluedistort[0]) != npix:
		print '!! warning bluenpix != bossnpix: %d != %d' % (len(bluedistort[0]),npix)
		print '!! skipping %s-%s...' % (plate,mjd)
		return -1
	# Plot Ratio
	savedir = os.path.join(BSR,blueVersion,plate)
	imageFilename = os.path.join(savedir,'spFluxdistortRatioImage-%s-%s.png'%(plate,mjd))
	plotImage(bossdistort/bluedistort,filename=imageFilename,\
		xlabel='Resampled Wavelength Bin (~3500-10500 Ang)',ylabel='Fiberid',\
		title='spFluxdistort-%s-%s Ratio\nDefault Pipeline / Pipeline Using Blue Standards'%(plate,mjd))
	# # Plot spectra
	# p.ioff()
	# p.clf()
	# fig = p.figure(1)
	# maxFiberid = 1000
	# randomFibers = random.sample(set(range(maxFiberid)),maxFiberid/10)
	# for i in randomFibers if not fiberlist else fiberlist:
	# 	color = colors[i%len(colors)]
	# 	alpha = 1
	# 	lw = 1
	# 	# Draw the spectra with the appropriate line properties
	# 	lines = p.plot(pixels,bossdistort[i]/bluedistort[i])
	# 	p.setp(lines,color=color,alpha=alpha,lw=lw)
	# # Label the figure
	# p.xlabel('Pixels')
	# p.ylabel('Default Pipeline / Pipeline Using Blue Standards')
	# p.title('spFluxdistort-%s-%s Ratio'%(plate,mjd))
	# p.axis([0,npix,0,2])
	# p.grid(True)
	# # Save the figure
	# savename = os.path.join(savedir,'spFluxdistortRatio-%s-%s.png'%(plate,mjd))
	# fig.savefig(savename, dpi=200)
def plotspFrameRatios(plate,mjd,fiberlist=[],verbose=False):
	# Look up exposure ids
	if verbose:
		print 'Looking up exposureIDs for %s-%s...' % (plate,mjd)
	exposureIDs = getExposureIDs(plate,mjd)
	# Iterate over exposure ids
	if verbose:
		print 'Found %d frames...' % len(exposureIDs)
	for exposureID in exposureIDs:
		bossCFramename = os.path.join(BSR,redVersion,plate,'spFrame-%s.fits.gz'%exposureID)
		blueCFramename = os.path.join(BSR,blueVersion,plate,'spFrame-%s.fits.gz'%exposureID)
		# Read standard boss spCFrame
		bosshdu = pf.open(bossCFramename)
		npix = bosshdu[0].header['naxis1']
		#wave = 10**(bosshdu[3].data[0]
		bossflux = bosshdu[0].data
		bossobjtype = bosshdu[5].data['objtype']
		alt = bosshdu[0].header['alt']
		seeing = bosshdu[0].header['SEEING50']
		airmass = bosshdu[0].header['airmass']
		#bossivar = bosshdu[1].data
		bosshdu.close()
		# Read blue standard calibrated spCFrame
		bluehdu = pf.open(blueCFramename)
		blueflux = bluehdu[0].data
		blueobjtype = bluehdu[5].data['objtype']
		#blueivar = bluehdu[1].data
		bluehdu.close()
		# Plot Ratio 
		savedir = os.path.join(BSR,blueVersion,plate)
		imageFilename = os.path.join(savedir,'spFrameRatioImage-%s.png'%exposureID)
		plotImage(bossflux/blueflux,filename=imageFilename,\
			xlabel='Pixel',ylabel='Fiberid',\
			title='spFrame-%s Ratio (%s-%s)\nAirmass = %.2f, Seeing50 = %.2f\nDefault Pipeline / Pipeline Using Blue Standards'%(exposureID,plate,mjd,airmass,seeing),\
			aspect=6)

def plotspCFrameRatios(plate,mjd,fiberlist=[],verbose=False):
	# Look up exposure ids
	if verbose:
		print 'Looking up exposureIDs for %s-%s...' % (plate,mjd)
	exposureIDs = getExposureIDs(plate,mjd)
	# Iterate over exposure ids
	if verbose:
		print 'Found %d frames...' % len(exposureIDs)
	for exposureID in exposureIDs:
		bossCFramename = os.path.join(BSR,redVersion,plate,'spCFrame-%s.fits'%exposureID)
		blueCFramename = os.path.join(BSR,blueVersion,plate,'spCFrame-%s.fits'%exposureID)
		# Read standard boss spCFrame
		bosshdu = pf.open(bossCFramename)
		npix = bosshdu[0].header['naxis1']
		wave = 10**bosshdu[3].data[0]
		bossflux = bosshdu[0].data
		bossobjtype = bosshdu[5].data['objtype']
		alt = bosshdu[0].header['alt']
		seeing = bosshdu[0].header['SEEING50']
		airmass = bosshdu[0].header['airmass']
		#bossivar = bosshdu[1].data
		bosshdu.close()
		# Read blue standard calibrated spCFrame
		bluehdu = pf.open(blueCFramename)
		blueflux = bluehdu[0].data
		blueobjtype = bluehdu[5].data['objtype']
		#blueivar = bluehdu[1].data
		bluehdu.close()
		# Plot Ratio 
		savedir = os.path.join(BSR,blueVersion,plate)
		# imageFilename = os.path.join(savedir,'spCFrameRatioImage-%s.png'%exposureID)
		# plotImage(bossflux/blueflux,filename=imageFilename,\
		# 	xlabel='Pixel',ylabel='Fiberid',\
		# 	title='spCFrame-%s Ratio (%s-%s)\nAirmass = %.2f, Seeing50 = %.2f\nDefault Pipeline / Pipeline Using Blue Standards'%(exposureID,plate,mjd,airmass,seeing),\
		# 	aspect=6)

		# Plot Individual Spectra
		# p.ioff()
		# p.clf()
		# fig = p.figure()
		# maxFiberid = 500
		# randomFibers = random.sample(set(range(maxFiberid)),maxFiberid/10)
		# colorcycler = cycle(colors)
		# for i in randomFibers if not fiberlist else (np.array(fiberlist) % maxFiberid):
		# 	color = next(colorcycler)
		# 	alpha = 1
		# 	lw = 1
		# 	# Draw the spectra with the appropriate line properties
		# 	lines = p.plot(wave,bossflux[i]/blueflux[i])
		# 	p.setp(lines,color=color,alpha=alpha,lw=lw)
		# # Save the figure
		# p.xlabel('Wavelength (Angstroms)')
		# p.ylabel('Default Pipeline / Pipeline Using Blue Standards')
		# p.title('spFluxcorr-%s Ratio (%s-%s)\nAlt = %.2f'%(exposureID,plate,mjd,alt))
		# p.axis([wave[0],wave[-1],0,2])
		# p.grid(True)
		# savename = os.path.join(savedir,'spCFrameRatio-%s.png'%exposureID)
		# fig.savefig(savename,dpi=200)


def plotspPlateRatios(plate,mjd,fiberlist=[],verbose=False):
	# Look up exposure ids
	if verbose:
		print 'Loading plate file for %s-%s...' % (plate,mjd)
	bossPlateFilename = os.path.join(BSR,redVersion,plate,'spPlate-%s-%s.fits'%(plate,mjd))
	bluePlateFilename = os.path.join(BSR,blueVersion,plate,'spPlate-%s-%s.fits'%(plate,mjd))
	# Read boss standard calibrated spPlate
	bosshdu = pf.open(bossPlateFilename)
	c0 = bosshdu[0].header['coeff0']
	c1 = bosshdu[0].header['coeff1']
	npix = bosshdu[0].header['naxis1']
	wave = 10.**(c0 + c1 * np.arange(npix))
	bossflux = bosshdu[0].data
	bossobjtype = bosshdu[5].data['OBJTYPE']
	bossivar = bosshdu[1].data
	bosshdu.close()
	# Read blue standard calibrated spPlate
	bluehdu = pf.open(bluePlateFilename)
	blueflux = bluehdu[0].data
	bluenpix = bluehdu[0].header['naxis1']
	blueobjtype = bluehdu[5].data['OBJTYPE']
	blueivar = bluehdu[1].data
	bluehdu.close()
	if bluenpix != npix:
		print '!! warning bossnpix != bluenpix: %d != %d' % (npix,bluenpix)
		print '!! skipping %s-%s...' % (plate,mjd)
		return -1
	# Plot Ratio
	savedir = os.path.join(BSR,blueVersion,plate)
	imageFilename = os.path.join(savedir,'spPlateRatioImage-%s-%s.png'%(plate,mjd))
	plotImage(bossflux/blueflux,filename=imageFilename,\
		xlabel='Resampled Wavelength Bin (~3500-10500 Ang)',ylabel='Fiberid',\
		title='spPlate-%s-%s Ratio\nDefault Pipeline / Pipeline Using Blue Standards'%(plate,mjd))

	# # Plot spectra
	# p.ioff()
	# p.clf()
	# fig = p.figure(1)
	# maxFiberid = 1000
	# randomFibers = random.sample(set(range(maxFiberid)),maxFiberid/10)
	# for i in randomFibers if not fiberlist else fiberlist:
	# 	color = colors[i%len(colors)]
	# 	alpha = .3
	# 	lw = 1
	# 	# if bossobjtype[i] == 'QSO':
	# 	# 	color = 'b'
	# 	# 	alpha = 0.1
	# 	# 	continue
	# 	# elif bossobjtype[i] == 'GALAXY':
	# 	# 	color = 'r'
	# 	# 	alpha = 0.0
	# 	# 	continue
	# 	# elif bossobjtype[i] == 'SPECTROPHOTO_STD':
	# 	# 	color = 'g'
	# 	# 	alpha = 1
	# 	# 	continue
	# 	# elif blueobjtype[i] == 'SPECTROPHOTO_STD':
	# 	# 	color = 'm'
	# 	# 	alpha = .3
	# 	# else:
	# 	# 	continue
	# 	# Draw the spectra with the appropriate line properties
	# 	lines = p.plot(wave,bossflux[i]/blueflux[i] * (blueivar[i] > 0))
	# 	p.setp(lines,color=color,alpha=alpha,lw=lw)
	# # Label the figure
	# p.xlabel('Wavelength (Angstroms)')
	# p.ylabel('Default Pipeline / Pipeline Using Blue Standards')
	# p.title('spPlate-%s-%s Ratio'%(plate,mjd))
	# p.axis([3500,6500,0,2])
	# p.grid(True)
	# # Save the figure
	# savename = os.path.join(savedir,'spPlateRatio-%s-%s.png'%(plate,mjd))
	# fig.savefig(savename, dpi=200)

if __name__ == '__main__':
	plateMJDlist = []
	for line in open('bluePlateMJDList.txt'):
		(plate, mjd) = line.strip().split()
		plateMJDlist.append((plate,mjd))
	#plateMJDlist = [('6298','56208')]
	maxFiberid = 500
	fiberlist = random.sample(set(range(maxFiberid)),20)
	print fiberlist

	plateMJDlist = [('6298','56208')]
	for plate,mjd in plateMJDlist:
		print 'Creating plots for %s-%s...' % (plate,mjd)
		plotspFluxcalibRatios(plate,mjd,fiberlist=fiberlist)
		#plotspFluxcorrRatios(plate,mjd,fiberlist=fiberlist)
		#plotspFluxdistortRatios(plate,mjd,fiberlist=fiberlist)
		plotspCFrameRatios(plate,mjd,fiberlist=fiberlist)
		#plotspFrameRatios(plate,mjd,fiberlist=fiberlist)
		#plotspPlateRatios(plate,mjd,fiberlist=fiberlist)


