FUNCTION ilc_wrap,axesArg,sunLatArg,sunDWLongArg,obsLatArg,obsWLong=obsWLong,$  
  nPixels=nPixels,spots=spots,bgAlbedo=bgAlbedo,scaling=scaling,debug=debugArg,$
  dmag=dmag,err=err,save=save

; this function is just a wrapper for the ilc.pro function.  it provides a few 
; extra options that may be useful to observers, including the option to convert the 
; light curve to a change in magnitude, add random gaussian noise, and save the output 
; images and plots to combined or individual files.  
; 
; INPUT:
; 
; axesArg =       ellipsoid's semimajor axes (x,y,z), arbitrary units
; sunLatArg =     sub-solar latitude (degrees)
; sunDWLongArg =  sub-solar delta longitude (west) wrt the observer (degrees);
;                 if sunLatArg=obsLatArg, then this value is the solar phase angle
; obsLatArg =     sub-observer latitude (degrees)
; obsWLong =      sub-ovserver longitude (west), default = 0
; nPixels =       size of the output images of the ellipsoid, default = 512
;                 
; OPTIONAL INPUT:
; 
; bgAlbedo =      albedo of the background
; spots =         an nspot x 4 array containing information on the surface spots;
;                 each column (spot) should contain 4 rows, each with a different piece of info: 
;                 1) spot center's latitude on the ellipsoid (degrees)
;                 2) spot center's longitude (west) wrt observer (degrees)
;                 3) spot radius, in units of the axes defined in axesArg
;                 4) spot albedo
; scaling =       if this keyword is set, the .img files will be normalized to each other, 
;                 just for display purposes
; dmag =          if this keyword is set, the fluxes will be converted to a differential
;                 magnitude relative to the mean
; err =           if the user wants to add gaussian noise, set this keyword to the desired
;                 error level
; save =          if the user wants to save the images and plots to a file, set this keyword;
;                 for a combined image-light curve plot, "save" should only have 1 element
;                 containing the desired filename of the plot; for 2 output files (one with
;                 the images and another with the light curve), "save" shoudl be a 2-element
;                 string array with the names of the output image and light curve files, resp.
;                 
; OUTPUT:
; 
; a database with the following information:
; 
;   { dib:          Total Disk-Integrated Brightnesses or (if dmag is set) differential magnitudes
;   , img:          Square array (image) of Disk-Resolved Brighnesses
;   , obsLatDeg:    Sub-observer latitude, degrees
;   , obsWLongDeg:  Sub-observer WEST longitude, degrees
;   , sunLatDeg:    Sub-solar latitude, degrees
;   , sunWLongDeg:   Sub-solar WEST longitude, degrees
;   }
;   
; EXAMPLES:
; 
;;; create a rotating ellipsoid with:
;;;   - dimensions 8x7x6,  
;;;   - background albedo of 0.6, 
;;;   - 37 realizations for one rotation (rotational resolution...)
;;;   - observer latitude of 10
;;;   - solar latitude of 0
;;;   - sun positioned 2-degrees west of the observer
;;;   - image size of 64x64 pixels
;;;   - 2 spots, with separate parameters: 
;;;     - -10deg lat., -30deg long. (W of observer), 1.5 units spot radius, albedo of 0.7
;;;     - 10deg lat., 40deg long. (W of observer), 0.8 units spot radius, albedo of 0.8
;;; 
; 
; .compile ilc    ;;; So createSpot will be available
; 
; spots = [ createSpot( -10, -30, 1.5, br=0.7),$
;           createspot( 10, 40, 0.8, br=0.8)]
; 
; data=ilc_wrap([8,7,6],0,2,10,obsWLong=indgen(37)*10,nPixels=64,bgAlbedo=0.7,$
;             spots=spots,/scaling)
;             
;;; create a rotating ellipsoid without spots, but with the following properties:
;;;   - dimensions 8x7x6
;;;   - 46 realizations for one rotation
;;;   - observer latitude of 0
;;;   - solar latitude of 0
;;;   - sun positioned 2-degrees west of the observer
;;;   - image size of 75x75 pixels
;;;   - convert the output to differential magnitudes
;;;   - add gaussian error of 0.01 magnitudes
;;;   - save the output as two separate files, one named "images" with the model images, and
;;;     one named "lightcurve" with the differential magnitude light curve
;
; data=ilc_wrap([8,7,6],0,2,0,obsWLong=indgen(45)*(360./45.),nPixels=75,/dmag,err=0.01,$
;             save=['images','lightcurve'])
;             
; WRITTEN: 4/24/13
; 

tv,make_array(!D.X_SIZE*2,!D.Y_SIZE*2,val=255b)

rtn=ilc(axesArg,sunLatArg,sunDWLongArg,obsLatArg,obsWLong=obsWLong,nPixels=nPixels,$
        bgAlbedo=bgAlbedo,spots=spots,scaling=scaling)

nLCPoints = n_elements(obsWLong)

if nLCPoints gt 0L then begin

  grayImg = make_array( 64, 64, value=127b )
  obsWLongDegs = rtn.obsWLongDeg
  for i=0,35 do begin
    wLong = i * 10d0
    tmp = min( abs(obsWLongDegs - wLong), iw)
    xOffset = 64 + (i MOD 9)*66
    yOffset = 64 + (i / 9)*66
    if tmp lt 5d0 then begin
      tvscl, rebin(rtn[iw].img,64,64), xOffset, yOffset
    endif else begin
      tv, dummy, xOffset, yOffset
    endelse
  endfor

endif

blue = ishft(255L,16)
!x.tickinterval = 90


if (n_elements(dmag) gt 0) then begin

  meanDib = mean(rtn.dib)

  dmags= -2.5*alog10(rtn.dib/meanDib)

  if (n_elements(err) gt 0) then begin
    rgauss=randomn(seed,n_elements(dmags))
    dmags += rgauss*err
    dmerr=make_array(n_elements(dmags),val=err)
  endif

  rtn.dib=dmags

  yTitle = 'Change in Magnitude'

endif else begin

  yTitle = 'Disk-Integrated Brightness (arbitrary units)'

  ptr_free,ptr_new(dmerr,/no_copy) ;;; Ensure dmerr is undefined

endelse

thisIsGdl,isGdl
if keyword_set(isGdl) then begin
  siSfx = '.jpg'
  siJpeg = 1
  ptr_free,ptr_new(siTiff,/no_copy)
endif else begin
  siSfx = '.tiff'
  siTiff = 1
  ptr_free,ptr_new(siJpeg,/no_copy)
endelse

if (n_elements(save) le 1) then begin

  tmp = ilc_wrap_plot( rtn.obsWLongDeg,rtn.dib,yTitle,dmerr,COLOR=blue,/NOERASE,/NOCLIP )

endif

if (n_elements(save) ge 1) then begin

  if keyword_set(isGdl) then begin
    write_jpeg,strtrim(save[0],2)+siSfx,tvrd(order=0,true=1)
  endif else begin
    saveimage,strtrim(save[0],2)+siSfx,tiff=siTiff,jpeg=siJpeg
  endelse

endif

if (n_elements(save) eq 2) then begin

  tmp = ilc_wrap_plot( rtn.obsWLongDeg,rtn.dib,yTitle,dmerr,COLOR=blue )

  set_plot,'ps'
  device,filename=strtrim(save[1],2)+'.ps',/landscape,/color,bits_per_pixel=8
  tmp = ilc_wrap_plot( rtn.obsWLongDeg,rtn.dib,yTitle,dmerr,COLOR=blue )
  device,/close
  set_plot,'x'

endif

return,rtn

END
