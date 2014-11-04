;+
; NAME:
;   plate_tpcorr
;
; PURPOSE:
;   Calculate throughput corrections
;
; CALLING SEQUENCE:
;   plate_tpcorr, plateid, ha, fwhm
;
; INPUTS:
;   plateid - plate number
;   ha - exposure hour angle
;   fwhm - exposure psf fwhm
;
; OPTIONAL INPUTS: 
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   22-Oct-2014  Written by Daniel Margala (dmargala@uci.edu), UC Irvine.
;-
pro plate_tpcorr, plateid, mjd, mapmjd, ha, fwhm, outfilename

; Assume guiding for 5400 Angstroms light
guideon=5400.
; Assume this is pointing number one and no offset
pointing=1L
offset=0L

; path to plPlugMap file
plugmapname = speclog+'/'+mapmjd+'/plPlugMapM-'+plateid+'-'+mapmjd+'-01.par'
plugmap = yanny_readone(plugmapname)

; Set path to directory for the specified plate number
platedir= plate_dir(plateid)

; Construct path to input plateHoles file, contains fiber positions and other relevant information
fullfile= platedir+'/plateHolesSorted-'+strtrim(string(f='(i6.6)',plateid),2)+'.par'
splog, 'Opening plateHoles file: '+fullfile
check_file_exists, fullfile, plateid=plateid

; Parse contents of input file, 
; the input file has a global header and a data entry with multiple fields for each fiber
full= yanny_readone(fullfile, hdr=phdr, /anon)
definition= lines2struct(phdr)
default= definition

; Read design hour angle, temperature
; Temperature is set per plate in the platePlans.par file
design_ha=float(strsplit(definition.ha, /extr))
temp=float(definition.temp)

design_platescale_alt=float(definition.design_platescale_alt)
mm_to_arcsec = 3600./design_platescale_alt
fiber_diameter = 2.0 

; Calculate ra and dec for center of the plate
; The values are accessible through 'racen' and 'deccen' after this
; This is necessary for handling multiple pointing plates (and 
; we keep it general for multi-offset plates).
plate_center, definition, default, pointing, offset, $
              racen=racen, deccen=deccen

; Select "good" targets/fibers, this is a list of array indices
igood=where(full.target_ra ne 0. or full.target_dec ne 0. and $
            full.pointing eq pointing, ngood)
; Collect ra/dec/lambda/x/y info for each target
ra= full[igood].target_ra
dec= full[igood].target_dec
lambda= full[igood].lambda_eff ;; e.g., 5400 for LRGs, 4000 for QSOs, 16600 for APOGEE targets
xforig= full[igood].xfocal
yforig= full[igood].yfocal

fiberids = plugmap[igood].fiberid

; Calculate xfocal and yfocal for this pointing (should be similar 
; to xforig/yforig up to round-off)
plate_ad2xy, definition, default, pointing, offset, ra, dec, $
             lambda, xf=xfocal, yf=yfocal, lst=racen+design_ha[pointing-1L], $
             airtemp=temp

; Grab all fibers with with lambda_eff at 5400. These targets are used to fit 
; for guiding parameters: rotation, scale, xshift, yshift.
ifit= where(full[igood].lambda_eff eq guideon, nfit)

;; Calculate xtmp, ytmp (all igood targets) at the exposure hour angle
plate_ad2xy, definition, default, pointing, offset, ra, dec, $
            lambda, xf=xtmp, yf=ytmp, lst=racen+ha, $
            airtemp=temp

;; Fit rotation, scale, shift parameters in guide targets
ha_fit, xfocal[ifit], yfocal[ifit], xtmp[ifit], ytmp[ifit], $
       xnew=xtmp2, ynew=ytmp2, rot=rot, scale=scale, $
       xshift=xshift, yshift=yshift

splog, format='(%"rot, scale, xshift, yshift: %f, %f, %f, %f",$)', rot, scale, xshift, yshift

; Only calculate offsets for targets with 4000 lambda_eff
focuson= 4000.
i4000= where(full[igood].lambda_eff eq focuson, n4000)
xfocal4000= xfocal[i4000]
yfocal4000= yfocal[i4000]

fiberids4000 = fiberids[i4000]

;; Apply rotation, scale, shift adjustments (i4000 targets)
ha_apply, xfocal4000, yfocal4000, xnew=xfocal4000_exp, ynew=yfocal4000_exp, rot=rot, scale=scale, $
         xshift=xshift, yshift=yshift

; Calculate xfocal and yfocal at the reference wavelength, where the targets
; would be if lambda_eff was 5400
lambda_ref= replicate(5400., n4000) 
plate_ad2xy, definition, default, pointing, offset, ra[i4000], dec[i4000], $
             lambda_ref, xf=xfocalref, yf=yfocalref, lst=racen+design_ha[pointing-1L], $
             airtemp=temp

;; Apply rotation, scale, shift adjustments (i4000 targets)
ha_apply, xfocalref, yfocalref, xnew=xfocalref_exp, ynew=yfocalref_exp, rot=rot, scale=scale, $
         xshift=xshift, yshift=yshift

; set up wavelengths to calculate offsets at
nlambda = 71L
lambda_min = 3500.
lambda_max = 10500.
lambda_array = lambda_min+(lambda_max-lambda_min)*(findgen(nlambda)/float(nlambda-1L))

; Create empty arrays to store guiding corrections and position offsets at each wavelength
tpcorr= fltarr(n4000, nlambda)

for i=0L, nlambda-1L do begin
  new_lambda = replicate(lambda_array[i], n4000) 
  ;; Calculate xtmp, ytmp at this wavelength
  plate_ad2xy, definition, default, pointing, offset, ra[i4000], dec[i4000], $
               new_lambda, xf=xtmp, yf=ytmp, lst=racen+ha, $
               airtemp=temp
  ;; Apply rotation, scale, shift adjustments
  ha_apply, xtmp, ytmp, xnew=xnew, ynew=ynew, rot=rot, scale=scale, $
            xshift=xshift, yshift=yshift
  ;; Calculate throughput for this wavelength at the reference hole position
  ref_offsets = mm_to_arcsec*sqrt((xnew-xfocalref_exp)^2+(ynew-yfocalref_exp)^2)
  tpref = fiberfraction(fwhm, ref_offsets, fiber_diameter)
  ;; Calculate throughput for this wavelength at the actual hole position
  offsets = mm_to_arcsec*sqrt((xnew-xfocal4000_exp)^2+(ynew-yfocal4000_exp)^2)
  tp = fiberfraction(fwhm, offsets, fiber_diameter)
  ;; Save throughput correction for this wavelength
  tpcorr[*,i] = tpref/tp
endfor


; open output file
splog, 'Saving throughput corrections to file: '+outfilename

openw, 1, outfilename

; iterate over targets
for i=0L, n4000-1L do begin
  ; print x,y location
	printf, 1, format='(%"%d %f %f ",$)', fiberids4000[i], xfocal4000[i], yfocal4000[i]
  ; print tabulated throughput corrections
	for j=0L, nlambda-1L do begin
 		   printf, 1, format='(%"%f ",$)', tpcorr[i,j]
	endfor
	printf, 1, format='(%"")'
endfor

close, 1

end
;------------------------------------------------------------------------------
