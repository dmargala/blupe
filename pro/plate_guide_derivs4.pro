;+
; NAME:
;   plate_guide_derivs
; PURPOSE:
;   Calculate guide derivatives as a function of HA
; CALLING SEQUENCE:
;   plate_guide_derivs, plateid [, pointing, guideon= ]
; INPUTS:
;   plateid - plate number
; OPTIONAL INPUTS:
;   pointing - pointing number
;   guideon - wavelength to guide on in Angstroms (default 5400.)
; COMMENTS:
;   Writes output to:
;    - plateGuideAdjust-XXXXXX-pP-lGUIDEON.par
;    - plateGuideOffsets-XXXXXX-pP-lGUIDEON.par
;   If there are no holes designed for the "guideon" wavelength,
;     then it bombs
; REVISION HISTORY:
;   10-Jun-2008  MRB, NYU
;-
pro plate_guide_derivs3, in_plateid, pointing, guideon=guideon

common com_plate_guide_derivs3, plateid, full, definition, default, phdr

; By default, optimize guiding for 5400 Angstroms light
; Assume this is pointing number one
if(NOT keyword_set(guideon)) then guideon=5400.
if(NOT keyword_set(pointing)) then pointing=1L
offset=0L

; If this plateid is the same as the last one, we wont read
; in the plate information again; but if it is a different one,
; we reset these variables (and will therefore read in new files 
; below)
if(keyword_set(plateid)) then begin
   if(in_plateid ne plateid) then begin
      full=0L
      definition=0L
      default=0L
   endif
endif 
plateid= in_plateid

; Set path to directory for the specified plate number
platedir= plate_dir(plateid)

; Construct paths for two output files:
;    ...plateGuideAdjust... : guider corrections for optimal guiding
;    ...plateGuideOffsets... : target positions on plate as function of HA
;post=string(f='(i6.6)', plateid)+ $
;     '-p'+strtrim(string(pointing),2)+ $
;     '-l'+strtrim(string(guideon, f='(i5.5)'),2)
;outdir='~/guiding'
;adjustfile= outdir+'/plateGuideAdjust-'+post+'.par'
;offsetfile=outdir+'/plateGuideOffsets-'+post+'.par'

; Construct path to input plateHoles file, contains fiber positions and other relevant information
fullfile= platedir+'/plateHolesSorted-'+ $
          strtrim(string(f='(i6.6)',plateid),2)+'.par'
check_file_exists, fullfile, plateid=plateid

; Parse contents of input file, 
; the input file has a global header and a data entry with multiple fields for each fiber
if(n_tags(full) eq 0) then begin
   full= yanny_readone(fullfile, hdr=phdr, /anon)
   definition= lines2struct(phdr)
   default= definition
endif

; Read design hour angle, temperature
; Temperature is set per plate in the platePlans.par file
ha=float(strsplit(definition.ha, /extr))
temp=float(definition.temp)

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
nextlambda= lambda + 1400.0
xforig= full[igood].xfocal
yforig= full[igood].yfocal

; Calculate xfocal and yfocal for this pointing (should be similar 
; to xforig/yforig up to round-off)
plate_ad2xy, definition, default, pointing, offset, ra, dec, $
             lambda, xf=xfocal, yf=yfocal, lst=racen+ha[pointing-1L], $
             airtemp=temp

; I'm pretty sure this would grab all fibers with with lambda_eff at 5400, so basically
;  everything except the QSOs
; I propose something like:
;  guidefibers = where(full[goodfibers].holetype eq 'GUIDE', nguide)
; MRB: I would check how stable this is first; it may be that my code
;  is actually worse behaved for low N than the guider is, and that is 
;  worth checking before using these for corrections; this is a
;  critical issue for APOGEE and can't be changed without
;  thorough testing.
ifit= where(full[igood].lambda_eff eq guideon, nfit)
if(nfit eq 0) then begin
   file_delete, adjustfile, /allow
   file_delete, offsetfile, /allow
   splog, 'No holes with LAMBDA_EFF='+strtrim(string(guideon),2)
   return
endif

; Set up hour angle window for guiding optimization and offset calculations
if(ha[pointing-1L] lt -120. OR $
   ha[pointing-1L] gt  120.) then begin
   message, 'HA desired is more than 120 deg! I refuse.'
endif

nha=17L
; 2012-01-05 reduce hour angle range for plots temporarily
minha= (ha[pointing-1L]-45.)
maxha= (ha[pointing-1L]+45.)
;minha= (ha[pointing-1L]-30.)>(-80.)
;maxha= (ha[pointing-1L]+30.)<(80.)
hatest= minha+(maxha-minha)*(findgen(nha)/float(nha-1L))


; set up wavelengths to calculate offsets at
nLambda = 25L
minLambda = 3500.
maxLambda = 10500.
lambdaArray = minLambda+(maxLambda-minLambda)*(findgen(nLambda)/float(nLambda-1L))

xmin = -350.
xmax = 350.
xnum = 10L

ngridpoints = xnum*xnum

xgrid = fltarr(ngridpoints)
ygrid = fltarr(ngridpoints)

for i=0L, xnum-1L do begin
  for j=0L, xnum-1L do begin
    xgrid[i + j*xnum] = xmin + i*(xmax-xmin)/float(xnum-1L)
    ygrid[i + j*xnum] = xmin + j*(xmax-xmin)/float(xnum-1L)
  endfor
endfor

; Create empty arrays to store guiding corrections and position offsets at each hour angle value
xfall= fltarr(ngridpoints, nha, nLambda)
yfall= fltarr(ngridpoints, nha, nLambda)
;xfall2= fltarr(ngood, nha)
;yfall2= fltarr(ngood, nha)

rot=fltarr(nha)
scale=fltarr(nha)
xshift=fltarr(nha)
yshift=fltarr(nha)

; Iterate over HA array
for i=0L, nha-1L do begin
   ;; Calculate xtmp, ytmp (all igood targets) at this hour angle
   plate_ad2xy, definition, default, pointing, offset, ra, dec, $
                lambda, xf=xtmp, yf=ytmp, lst=racen+hatest[i], $
                airtemp=temp
   ;; Fit rotation, scale, shift parameters in guide targets
   ha_fit, xfocal[ifit], yfocal[ifit], xtmp[ifit], ytmp[ifit], $
           xnew=xtmp2, ynew=ytmp2, rot=rottmp, scale=scaletmp, $
           xshift=xshifttmp, yshift=yshifttmp
   ;; Save rotation, scale, shift parameters at this hour angle
   rot[i]=rottmp
   scale[i]=scaletmp
   xshift[i]=xshifttmp
   yshift[i]=yshifttmp
   ;; Apply rotation, scale, shift adjustments (all igood targets)
   ha_apply, xtmp, ytmp, xnew=xnew, ynew=ynew, rot=rot[i], scale=scale[i], $
             xshift=xshift[i], yshift=yshift[i]
   for k=0L, nLambda-1L do begin
         
      ;; Save x,y position (all igood targets) at this hour angle
      ;xfall[*,i,k]= xnew
      ;yfall[*,i,k]= ynew
      nextLambda = lambda
      nextLambda[igood] = lambdaArray[k] 
   ;; Calculate xtmp, ytmp (all igood targets) at this hour angle
      plate_ad2xy, definition, default, pointing, offset, ra, dec, $
                   nextLambda, xf=xtmp, yf=ytmp, lst=racen+hatest[i], $
                   airtemp=temp
   ;; Fit rotation, scale, shift parameters in guide targets
   ;ha_fit, xfocal[ifit], yfocal[ifit], xtmp[ifit], ytmp[ifit], $
   ;        xnew=xtmp2, ynew=ytmp2, rot=rottmp, scale=scaletmp, $
   ;        xshift=xshifttmp, yshift=yshifttmp
   ;; Save rotation, scale, shift parameters at this hour angle
   ;rot[i]=rottmp
   ;scale[i]=scaletmp
   ;xshift[i]=xshifttmp
   ;yshift[i]=yshifttmp
   ;; Apply rotation, scale, shift adjustments (all igood targets)
      ha_apply, xgrid, ygrid, xnew=xnew, ynew=ynew, rot=rot[i], scale=scale[i], $
                xshift=xshift[i], yshift=yshift[i]
   ;; Save x,y position (all igood targets) at this hour angle
      xfall[*,i,k]= xnew
      yfall[*,i,k]= ynew
   endfor
endfor

; open output file
outdir='~/guiding'
outfilename= outdir+'/plate-'+strtrim(string(plateid),2)+'.dat'
openw, 1, outfilename

; iterate over targets
for i=0L, ngridpoints-1L do begin
	printf, 1, $
		format='(%"%f %f ",$)', xgrid[i], ygrid[i]
   	; iterate over hour angles
   	for j=0L, nha-1L do begin
		for k=0L, nLambda-1L do begin
   		   ;altdirection = atan(xfall2[i,j]-xfall[i,j],yfall2[i,j]-yfall[i,j])
   		   printf, 1, $
   			   format='(%"%f %f %f %f ",$)', $
   			   hatest[j], lambdaArray[k], xfall[i,j,k]-xgrid[i], yfall[i,j,k]-ygrid[i]
		endfor
	endfor
	printf, 1, format='(%"")'
endfor

close, 1

end
;------------------------------------------------------------------------------
