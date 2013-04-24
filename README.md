See ilc.pro for more detail

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ILC.PRO - IDL LightCurve generator
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; FUNCTION ILC
;;;
;;; Calculate Disk-Resolved and Disk-Integrated Brightness (dib) for
;;; ellipsoid given lighting and viewing geometries, assuming Sun and
;;; observer are at infinite distances.
;;;
;;; N.B. All input longitudes are degrees WEST!
;;;
;;; Input arguments:
;;;
;;;   axesArg        Ellipsoid semi-major axes, XYZ
;;;   sunLatArg      Sub-solar latitude
;;;   sunDWLongArg   Sub-solar delta West longitude wrt observer
;;;   obsLatArg      Sub observer latitude
;;;
;;; Input keywords:
;;;
;;;   obsWLong       Sub-observer West longitude(s), default=0
;;;   nPixel         Size of image, default=512
;;;   albedoMap      Albedo map (not yet implemented)
;;;   PhotFunc       String; name of photometric func
;;;   Arg2PhotFunc   Optional argument to PhotFunc
;;;   debug         
;;;
;;; If multiple observer longitudes are supplied (argument obsWLong), then
;;; function ILC returns information for each observer longitude as
;;; an element in an array.
;;;
;;; Returned datum (for each input observer longitude) is an IDL structure:
;;;
;;;   { dib:          Total Disk-Integrated Brightnesses
;;;   , img:          Square array (image) of Disk-Resolved Brighnesses
;;;   , obsLatDeg:    Sub-observer latitude, degrees
;;;   , obsWLongDeg:  Sub-observer WEST longitude, degrees
;;;   , sunLatDeg:    Sub-solar latitude, degrees
;;;   , sunWLongDeg:   Sub-solar WEST longitude, degrees
;;;   }
;;;
;;; Example:
;;;
;;; - Get lightcurve of [5,3,1] triaxial ellipsoid at 37 sub-observer
;;;   West Longitudes at 10-degree increments (0,10,20,...,360), with 
;;;   - sub-solar latitude = -10deg (10S)
;;;   - sub-solar West longitude = 40deg West of sub-observer point
;;;   - sub-observer latitude = 20deg (20N)
;;;   - 64x64 pixel images
;;;
;;;     rtn = ilc([5,3,1],-10,40,20 ,obsWLong=indgen(37)*10 ,nPixels=64)
;;;
;;; - Display the first of those images of Disk-Resolved Brightnesses:
;;;
;;;     tvscl, rtn[0].img
;;;
;;; - Plot lightcurve:  Observer West Longitude vs. Disk-Integrated
;;;                     Brightnesses
;;;
;;;     plot, rtn.obsWLong, rtn.dib
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
