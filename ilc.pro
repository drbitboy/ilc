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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Utility functions:

pro rpdCommonSet
common rpdCommon, rpd, dpr
  if n_elements(rpd) eq 1L then return
  rpd = !dpi / 180d0
  dpr = 180d0 / !dpi
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Utility vector[3] functions

;;; - Return dot product of two vectors
function vdot,v1,v2
  return, total(v1*v2,1)
end

;;; - Return length squared of a vector
function vnormsq,v
  return, vdot(v,v)
end

;;; - Return length of a vector
function vnorm,v
  return, sqrt(vnormsq(v))
end

;;; - Return unit vector parallel to and in same direction as a vector
function vhat,v
  v2 = vnormsq(v)
  if v2 gt 0d0 then return, v/sqrt(v2)
  return, v
end

;;; N.B. v[3] X v2[3] cross product is already part of IDL library (crossp)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Utility array[N,3] functions
;;;
;;; The functions with names that have "nx3_" as a prefix work on
;;; arguments that each contain mulitiple 3-D vectors stored in arrays
;;; of shape [N,3], where N is the number of vectors in the array.
;;;
;;; For example, if you had four vectors, [1,2,3], [4,5,6], [7,8,9],
;;; [10,11,12], then they would be arranged in an array of shape
;;; [4,3] like this:
;;;
;;;       1       4       7      10
;;;       2       5       8      11
;;;       3       6       9      12
;;;
;;; The first index (0 through 3) selects the vector.
;;; The second index (0 through 2) selects the X, Y or Z component.
;;; Each vector is an array of shape [1,3].
;;; 
;;; These nx3_* functions treat the vectors in the array arguments
;;; independently.
;;;
;;; If there is more than one argument, e.g. v1 and v2 for the dot (or
;;; cross) products, then all arguments must have the same shape [N,3], and
;;; the vector dot (or cross) product of the single vectors v1[i,0:2] and
;;; v2[i,0:2] will be element i of the returned [N] array
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; - Return dot products of vector pairs
;;;   Inputs:  N 3-D vector pairs as two [N,3] arrays
;;;   Return value:   N dot products as one [N] array
function nx3_vdot,v1,v2
  return, total(v1*v2,2)
end

;;; - Return lengths squared of vectors
;;;   Inputs:  N 3-D vectors as one [N,3] array
;;;   Return value:   N lengths squared as one [N] array
function nx3_vnormsq,v
  return, nx3_vdot(v,v)
end

;;; - Return lengths of vectors
;;;   Inputs:  N 3-D vectors as one [N,3] array
;;;   Return value:   N lengths as one [N] array
function nx3_vnorm,v
  return, sqrt(nx3_vnormsq(v))
end

;;; - Return unit vectors of vectors, from [N,3] array argument,
;;;   Inputs:  N 3-D vectors as one [N,3] array
;;;   Return value:   N 3-D unit vectors as one [N,3] array
function nx3_vhat,v
  vRtn = v
  v2 = nx3_vnormsq(v)
  iw = where( v2 gt 0d0, ct)
  if ct gt 0 then vRtn[iw,*] = v[iw,*] / ([1d0,1,1] ## sqrt(v2[iw]))
  return, vRtn
end

;;; - Return cross products of vector pairs
;;;   Inputs:  N 3-D vector pairs as two [N,3] arrays
;;;   Return value:   N 3-D vectors as one [N,3] array
function nx3_crossp,v1,v2
  return, [ [ (v1[*,1] * v2[*,2]) - (v1[*,2]*v2[*,1]) ] $
          , [ (v1[*,2] * v2[*,0]) - (v1[*,0]*v2[*,2]) ] $
          , [ (v1[*,0] * v2[*,1]) - (v1[*,1]*v2[*,0]) ] $
          ]
end

function createSpot, Lat, WLong, Radius, brightness=bArg
  return, { LatDeg: double(Lat[0]) $
          , WLongDeg: double(WLong[0]) $
          , Radius: double(Radius[0]) $
          , Brightness: n_elements(bArg) eq 1L ? double(bArg[0]) : 0d0 $
          , Center: [0d0,0,0] $
          }
end


;;; Convert spot center's body-centric Latitude and West Longitude to
;;; the ellipsoid surface point using the ellipsoid axes

pro spotFiddle, axesIn, spotIn
common rpdCommon
  rpdCommonSet

  rWLong = rpd * spotIn.WLongDeg
  rLat = rpd * spotIn.LatDeg
  cosLat = cos(rLat)
  spotXYZ = [cos(rWLong)*cosLat, -sin(rWLong)*cosLat, sin(rLat)]
  spotIn.Center = spotXyz / sqrt( vnormsq(spotXYZ/axesIn) )
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; FUNCTION ILC - see description above, at start of file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function ilc $
, axesArg $               ;;; Ellipsoid semi-major axes, XYZ
, sunLatArg $             ;;; Sub-solar latitude
, sunDWLongArg $          ;;; Sub-solar delta longitude (West) wrt obs
, obsLatArg $             ;;; Sub observer latitude
, obsWLong=obsWLongArg $  ;;; Sub observer longitude(s) (West), dflt=0
, nPixels=nPixelsArg $    ;;; Size of image, dflt=512
, albedoMap=albedoMap $   ;;; Albedo map (not yet implemented)
, PhotFunc=PhotFunc $     ;;; String; name of photometric func
, Arg2PhotFunc=Args2PhotFunc $  ;;; Optional argument to PhotFunc
, spots=spots $           ;;; Surface spots
, bgAlbedo=bgAlbedo $     ;;; Background albedo (0 to 1)
, scaling=scaling $       ;;; if this keyword is set, the .img files will be normalized to each other
, debug=debugArg

common rpdCommon

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; If the observer longitude argument, obsWLongArg, has multiple values,
  ;;; then call self (ilc) recursively with each value, and return a vector
  ;;; of results

  nLongitudes = n_elements(obsWLongArg)
  if nLongitudes gt 1L then begin

    for iLongitude=0L,nLongitudes-1L do begin
      iRtn = ilc( axesArg $
                , sunLatArg,sunDWLongArg,obsLatArg $
                , obsWLong=obsWLongArg[iLongitude] $
                , nPixels=nPixelsArg $
                , albedoMap=albedoMap $
                , PhotFunc=PhotFunc $
                , Arg2PhotFunc=Args2PhotFunc $
                , spots=spots $
                , bgAlbedo=bgAlbedo $
                , scaling=scaling $
                , debug=debugArg $
                )
      if iLongitude eq 0L then begin
        rtn = replicate( iRtn, nLongitudes )
      endif else begin
        rtn[iLongitude] = iRtn
      endelse
    endfor

    return,rtn

  endif

  ;;; N.B. See below for detailed explanation of PhotFunc and Arg2PhotFunc

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Set rpd and dpr (Radians Per Degree and reciprocal)
  rpdCommonSet
  dodebug = keyword_set(debugArg)


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Get ellipsoid semi-major axes' lengths, maximum axis length
  axes = double(axesArg[0:2])
  maxes = max(axes)

  ;;; Scaling from unit sphere system to ellipsoid system, and the reverse
  sph2ell = axes
  ell2sph = 1d0 / sph2ell


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; This section of calculatations are in the ellipsoid coordinate
  ;;; system; most vector and quantity variables start with an 'e'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Sub-observer and -solar West longitude and latitude
  obsWLong = rpd * double( n_elements(obsWLongArg) eq 1L ? obsWLongArg[0] : 0)
  obsLat = rpd * double( obsLatArg[0] )
  sunLat = rpd * double( sunLatArg[0] )
  ;;; Sub-solar longitude argument is an offset from sub-observer longitude
  sunWLong = rpd * double( sunDWLongArg[0] ) + obsWLong

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Convert West longitude and latitudes into unit 3-D vectors
  eObs = [ [cos(obsWLong),-sin(obsWLong)]*cos(obsLat), sin(obsLat) ]
  eSun = [ [cos(sunWLong),-sin(sunWLong)]*cos(sunLat), sin(sunLat) ]


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Pixels are represented as rays in space; this routine finds the
  ;;; intersection of each pixel's ray with the ellipsoid and calculates
  ;;; the brightness at that center of each pixel.  A ray comprises a
  ;;; position in space and a direction in space; both are vectors; the
  ;;; position's origin is the reference frame origin i.e the center of the
  ;;; ellipsoid; the direction is a unit vector whose "origin" is the
  ;;; position of the ray.
  ;;; - The observer is modeled as being at an infinite distance from
  ;;;   the ellipsoid so all pixel rays' directions are parallel i.e.
  ;;;   all are -eObs

  ;;; Build a grid of 3-D positions on a plane to represent image pixels:
  ;;; - The plane will
  ;;;   - be normal to eObs, and
  ;;;   - pass through the origin
  ;;; - There will be two unit grid vectors in that plane used to build
  ;;;   the grid:  one representing Up; one representing Right

  ;;; - Initially attempt to define "UP" normal to a plane defined by
  ;;;   three points:  observer; ellipsoid; sun.
  eUp = crossp( eObs, eSun )

  ;;; - Handle special case of zero phase (collinear obs.-ellipsoid-sun)
  ;;;   - I.e. eUp, as calculated above is zero-length
  ;;;   - make Up coincident with +Z in image plane
  if vnormsq(eUp) eq 0d0 then begin
    right = crossp( [0d0,0,1], eObs )
    eUp = crossp( eObs, right )
  endif
  ;;;   - OR, if observer and sun are on the Z axis, then set +X to the
  ;;;     right
  if vnormsq(eUp) eq 0d0 then begin
    eUp = crossp( eObs, [1d0,0,0] )
  endif

  ;;; - eRight is cross product of eUp and eObs (observer) vectors
  ;;; - scale eUp and eRight into unit vectors
  eUp = vhat(eUp)
  eRight = vhat( crossp(eUp,eObs) )


  ;;; - Convert nPixels argument into even number of pixels in square array
  nPixels = ( long( n_elements(nPixelsArg) eq 1L ? nPixelsArg[0] : 512 ) / 2L ) * 2L
  nPixelSqs = nPixels * nPixels

  ;;; - Create [N*N] vector to expand single 3-vectors into [N*N,3] arrays
  nPixelSqOnes = make_array( nPixelSqs, value=1d0)

  ;;; - Values and vectors for NxN pixel image will be a linear (N*N) vector
  ;;; - will be broken into a 2-D NxN array later

  ;;; - calculate offsets (N/2-.5) up, down, left, right from origin
  ;;;   e.g. if nPixels=512, max offset will be 255.5
  np2 = (nPixels-1) / 2d0

  ;;; - get scaling factor vectors for eUp and eRight

  eupScales = ( (lindgen(nPixelSqs) / nPixels) - np2 ) * maxes / np2
  erightScales = ( (lindgen(nPixelSqs) mod nPixels) - np2 ) * maxes / np2

  ;;; - Create [N*N,3] array of pixel positions, as each passes through a
  ;;;   plane, which is perpendicular to the observer direction and passes
  ;;;   through the origin:

  ePixels = (eUp ## eupScales) + (eRight ## erightScales)

  ;;; N.B. The ## operator is a matrix multiplication operator in IDL;
  ;;;      it multiplies the columns in the first array by the rows in
  ;;;      the second.
  ;;;
  ;;;      1) E.g. a general example:
  ;;;
  ;;;        [1 2]      [ 7  8  9 10]      [ 29  32  35  38]
  ;;;        [3 4]  ##  [11 12 13 14]  =>  [ 65  72  79  86]
  ;;;        [5 6]                         [101 112 123 134]
  ;;;
  ;;;      - the 29 in the result is (1*7 + 2*11)
  ;;;        - 1st arg 1st row & 2nd arg 1st column
  ;;;
  ;;;      - the 79 in the result is (3*8 + 4*13)
  ;;;        - 1st arg 2nd row & 2nd arg third column
  ;;;
  ;;;      2) E.g. a specific example applicable to the code in this file:
  ;;;
  ;;;        [1]                       [ 4   5   6   7.5]  ;; Each column
  ;;;        [2]  ##  [4 5 6 7.5]  =>  [ 8  10  12  15  ]  ;; is a multiple
  ;;;        [3]                       [12  15  18  22.5]  ;; of [1,2,3]
  ;;;
  ;;;         That shows that the ## operator can by used to scale a single
  ;;;         3-D column vector into multiple 3-D column vectors each of
  ;;;         which is a multiple of the single vector.
  ;;;
  ;;;      3) The use of the ## operation in this file takes advantage of
  ;;;         an apparent feature in IDL whereby, if the first array is
  ;;;         is not a 2-D array but has only one row, and the second
  ;;;         array has only one row, then it automatically "does the right
  ;;;         thing" by transposing the first array before doing the "##"
  ;;;         matrix multiplcation e.g. for a 3-D vector [X,Y,Z]:
  ;;;
  ;;;                                     [aX  bX  cX  dX]
  ;;;         [X Y Z]  ##  [a b c d]  =>  [aY  bY  cY  dY]
  ;;;                                     [aZ  bZ  cZ  dZ]


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Convert to sphere coordinate system; variables start with 's'
  ;;; - Simple linear scalings, by reciprocals of ellipsoid semi-major
  ;;;   axes XYZ lengths, of XYZ coordinates in ellipsoid system
  ;;;   maintains lines as lines but in a system where original
  ;;;   ellipsoid surface is the surface of a unit sphere, and for which
  ;;;   the intersection calculations are straightforward

  ;;; - convert pixel locations
  sPixels = ePixels * (ell2sph ## nPixelSqOnes)

  ;;; - convert [3] unit vector direction to observer to unit vector
  ;;;   in sphere system
  sobs = vhat(eObs*ell2sph)


  ;;; - vector to observer (sobs) and plane of pixels (sPixels) may no
  ;;;   longer be perpendicular; subtract offsets from sPixels so that will
  ;;;   be true again; offsets are dot products of unit vectors toward
  ;;;   observer and vectors to each pixel in non-perpendicular plane

  sPixels -= sobs ## nx3_vdot( sobs ## nPixelSqOnes, sPixels )


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Determine which pixels are inside the unit sphere
  ;;; - Calculate lengths squared of pixel vectors from origin
  sNearLengthSqs = nx3_vnormsq( sPixels )

  ;;; - Those pixels with lengths squared less than one (the radius of a
  ;;;   unit sphere) must intersect the surface of the sphere:
  iwInside = where( sNearLengthSqs lt 1d0, ctInside )


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Two ellipsoid system vectors, eTangent0 and eTangent1, tangent to
  ;;; ellipsoid surface at each pixel position, will be used to calculate
  ;;; ellipsoid surface normals.  Set sphere system analogs initially to
  ;;; zero; non-zero tangent vectors at points of intersection with pixels
  ;;; vectors will be set below.

  sTangent0 = make_array( nPixelSqs, 3, val=0d0)
  sTangent1 = sTangent0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Perform calculations for the pixels which were inside the unit sphere

  if ctInside gt 0L then begin

    ;;; - Calculate offsets toward actual points of intersections via
    ;;;   Pythagorus and add them to sPixels to get intersections
    sOffsets = sqrt( 1d0 - sNearLengthSqs[iwInside])
    sIntersects = sPixels[iwInside,*] + (sobs ## sOffsets)

    ;;; - We need the normal at each intersected point, not the point
    ;;;   - N.B. iwInside contains the placement of each intersected point
    ;;;          within the array of pixels
    ;;; - The normal is parallel to the ray from the origin to the point in
    ;;;   the sphere system, but *NOT* in the ellipsoid system, so a simple
    ;;;   conversion of the sphere normal back to the ellipsoid system will
    ;;;   not work.  Instead, vectors tangent to the sphere surface will be
    ;;;   perpendicular to the normal in both the sphere and ellipsoid
    ;;;   systems, and the cross product of two such tangent vectors, when
    ;;;   converted from the sphere to the ellipsoid system, will yield the
    ;;;   ellipsoid surface normal.

    ;;; - Start with cross products of the point, as a stand-in for the
    ;;;   normal, with the unit X and Y vectors; these will be the tangents
    ;;;   used for the case where the Z component is the largest component

    sInsideTangent0 = [ [ make_array(ctInside,val=0d0) ] $
                      , [ sIntersects[*,2] ] $
                      , [-sIntersects[*,1] ] $
                      ]

    sInsideTangent1 = [ [-sIntersects[*,2] ] $
                      , [ make_array(ctInside,val=0d0)] $
                      , [ sIntersects[*,0] ] $
                      ]

    ;;; - for the cases where Z is not the largest component, replace
    ;;;   either the cross product with X in sInsideTangent0 or the cross
    ;;;   product with Y in sInsideTangent1 with the cross product with Z

    sIntAbs = abs(sIntersects)
    iwNotXY = where( sIntAbs[*,2] lt max(sIntAbs[*,0:1], dim=2), ctNotXY )

    if ctNotXY gt 0L then begin

      iwXltY = where( sIntAbs[iwNotXY,0] lt sIntAbs[iwNotXY,1], ctXltY, complement=iwYleX )

      if ctXltY gt 0L then begin
        sInsideTangent1[iwNotXY[iwXltY],*] = [ [ sIntersects[iwNotXY[iwXltY],1] ] $
                                             , [-sIntersects[iwNotXY[iwXltY],0] ] $
                                             , [ make_array(ctXltY,val=0d0) ] $
                                             ]
      endif

      if ctXltY lt ctNotXY then begin
        sInsideTangent0[iwNotXY[iwYleX],*] = [ [ sIntersects[iwNotXY[iwYleX],1] ] $
                                             , [-sIntersects[iwNotXY[iwYleX],0] ] $
                                             , [ make_array(ctNotXY-ctXltY,val=0d0) ] $
                                             ]
      endif

    endif  ;;; if ctNotXY gt 0L

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;; Cross product of tangents will be parallel to normals, but some may
    ;;; be inverted; dot products of [cross products of tangents] and
    ;;; intersection points will be negative for those that need inverting.

    sNormVdots = nx3_vdot( nx3_crossp(sInsideTangent0, sInsideTangent1), sIntersects )

    ;;; - Invert sInsideTangent0 for any negative dot products

    iwNeg = where( sNormVdots lt 0d0, ctNeg )
    if ctNeg gt 0L then begin
      sInsideTangent0[iwNeg,*] = - sInsideTangent0[iwNeg,*]
    endif

    sTangent0[iwInside,*] = sInsideTangent0
    sTangent1[iwInside,*] = sInsideTangent1

  endif ;;; if ctInside gt 0L


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Convert sphere sytem tangent vectors back to ellipsoid system

  toEll = sph2ell ## nPixelSqOnes
  eTangent0 = sTangent0 * toEll
  eTangent1 = sTangent1 * toEll


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Cross product of tangent vectors will yield ellipsoid surface
  ;;  normals at each pixel; pixels which do not intersect the ellipsoid
  ;;  will have zero normals
  ;;  - convert normals to unit vectors
  ;;  - zero-length tangents [0,0,0] will yield zero-length normals

  eNormals = nx3_vhat( nx3_crossp( eTangent0, eTangent1 ) )


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Dot product of unit normal with unit direction to Sun is the cosine
  ;;; of the incidence angle, MuNaught (Mu0).
  ;;; Dot product of unit normal with unit direction to Observer is the
  ;;; cosine of the emission angle, Mu
  ;;; - Zero-length normals will yield zeros for MuNaught and Mu
  ;;; - Set negative values to zero (see "> 0d0" at the end)

  eMuNaught = nx3_vdot( eSun ## nPixelSqOnes, eNormals ) > 0d0
  eMu = nx3_vdot( eObs ## nPixelSqOnes, eNormals ) > 0d0


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Ensure all MuNaught values are zero where Mu values are zero

  iwZero = where( eMuNaught eq 0d0, ctZero )
  if ctZero gt 0L then begin
    eMu[iwZero] = 0d0
  endif


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Photometric function:  argument PhotFunc
  ;;; - PhotFunc is a string that is the name of an IDL function
  ;;;
  ;;; - Function get three arguments:
  ;;;   - eMuNaught - vector of cosine(incidenceAngle) values
  ;;;   - eMu       - vector of cosine(emissionAngle) values
  ;;;   - eAlpha    - scalar phase angle, radians
  ;;; - eMuNaught and eMu should be vectors of the same length
  ;;; - PhotFunc returns vector of brightnesses (I/F or arbitrary scale)
  ;;;   of the same length as eMuNaught and eMu
  ;;;
  ;;; - Function also takes one optional argument:
  ;;;   - Arg2PhotFunc - additional parameter (e.g. Hapke parameters)
  ;;;   - may be a structure or vector or anything; contents are
  ;;;     dependent on PhotFunc and outside the scope of these comments
  ;;;   
  ;;; - Albedo Map is not yet implemented; probably will be handled in
  ;;;   sphere system section above

  if n_elements(PhotFunc) eq 1L then begin
    eAlpha = acos( (vdot(eObs,eSun) < 1d0) > (-1d0) )
    eBrightnesses = call_function( PhotFunc0, eMuNaught, eMu, eAlpha, Arg2PhotFunc)
  endif else begin
    ;;; Default to Lambert surface if no photometric function supplied
    eBrightnesses = eMuNaught
  endelse

  if (n_elements(bgAlbedo) gt 0) then eBrightnesses = eBrightnesses * double(bgAlbedo)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Adjust for "spots" on ellipsoid surface
  ;;; Only adjust if there are input spots and where brightnesses are not
  ;;; zero

  nSpots = n_elements(spots)

  if ctInside gt 0L then begin
    iwNonZero = where( eBrightnesses[iwInside] gt 0d0 and nx3_vnormsq(sIntersects) gt 0d0, ctNonZero)
  endif else begin
    ctNonZero = 0L
  endelse

  if ctNonZero gt 0L and nSpots gt 0L then begin

    ;;; Convert sPixels with non-zero brightnesses from sphere to ellipsoid
    ;;; system

    ePixelsNonZero = sIntersects[iwNonZero,*] * (sph2ell ## nPixelSqOnes[0:ctNonZero-1L])

    ;;; Loop over input spots

    for iSpot=0L,nSpots-1L do begin
      spot = spots[iSpot]

      ;;; Calculate this spot center on ellipsoid surface per ellipsoid
      ;;; axes

      spotFiddle, axes, spot

      ;;; Calculate offsets from this spot's center to pixel positions,
      ;;; find which pixels are in the spot i.e. within this spot's radius
      ;;; of its center

      eSpotOffsets = ePixelsNonZero - (spot.Center ## nPixelSqOnes[0:ctNonZero-1L])
      iwInSpot = where( nx3_vnormsq( eSpotOffsets) le (spot.Radius^2), ctInSpot, complement=iwOutSpot )
      ;;;iwOutSpot=where(nx3_vnormsq(eSpotoffsets) gt (spot.Radius^2),ctOutSpot)

;      print,spot.center
;      pause=get_kbrd()

      ;;; Scale brightnesses of pixels in this spot

      eBrightnesses[iwInside[iwNonZero[iwInSpot]]] *= spot.Brightness

      if ctInSpot lt ctNonZero and n_elements(bgAlbedo) gt 0 then begin
        eBrightnesses[iwInside[iwNonZero[iwOutSpot]]] = eBrightnesses[iwInside[iwNonZero[iwOutSpot]]]*bgAlbedo
      endif

    endfor ;;; for iSpot=0L,nSpots-1L
  endif ;;; if ctNonZero gt 0L and nSpots gt 0L

  if n_elements(bgAlbedo) gt 0 AND nSpots eq 0L then begin
    eBrightnesses=eBrightnesses*bgAlbedo
  endif

  if keyword_set(scaling) then begin
    array = [0d0]
    if n_elements(bgalbedo) gt 0L then array = [ array, bgalbedo ]
    if nSpots gt 0L then array = [ array, spots.Brightness ]
    eBrightnesses[0]=max(array)
  endif

  img = make_array( nPixels, nPixels, val=0d0 )
  img[*] = eBrightnesses


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; N.B. Abbreviation DIB means Disk-Integrated Brightness

  return, { dib: total(eBrightnesses) $
          , img: img $
          , obsLatDeg: obsLat*dpr $
          , obsWLongDeg: obsWLong*dpr $
          , sunLatDeg: sunLat*dpr $
          , sunWLongDeg: sunWLong*dpr $
          }
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Test code:  .run ilc.pro
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!quiet=1

;;; Get lightcurve of [5,3,1] triaxial ellipsoid at 37 sub-observer
;;; West Longitudes at 10-degree increments (0,10,20,...,360), with 
;;; - sub-solar latitude = -10deg (10S)
;;; - sub-solar longitude = 40deg West of sub-observer
;;; - sub-observer latitude = 20deg (20N)
;;; - 64x64 pixel images



;catcherr=0L
;;;;catch,catcherr
;if catcherr eq 0L then begin
;  print,!error_state
;  device,ret=2
;  erase, color=-1L                       ;;; Clear display to white
;  catch,cancel
;endif else begin
;  catch,cancel
  tv,make_array(!D.X_SIZE*2,!D.Y_SIZE*2,val=255b)
;endelse

if keyword_set(doBrian) then begin
  ;;; Brian's test

  spots = [ createSpot( 12, 170, 1.5) $
          , createSpot( 12, 100, 2.0) $
          ]

  rtn = ilc([7,5,3],2,-70,20 ,obsWLong=indgen(37)*10 ,nPixels=64, spots=spots)

  for i=0,35 do $
    tvscl, rtn[i].img $
         , 64+(i MOD 9)*66, 64+(i/9)*66  ;;; display 36 images

  !x.tickinterval = 90
  plot,rtn.obsWLongDeg,rtn.dib $
      , XSTYLE=1 $
      , PSYM=-6 $
      , SYMSIZE=1.5 $
      , XTITLE='Sub-observer longitude, degW' $
      , YTITLE='Disk-Integrated Brightness, arbitrary units' $
      , COLOR=ishft(255L,16) $
      , /NOERASE, /NOCLIP  ;;; Over-plot lightcurve in blue

endif else begin
  ;;; OR Sarah's test

  spots = [ createSpot( 20, 10, 2.8, br=0.7) ]

  rtn=ilc([4,4,4],0,10,0,obsWLong=indgen(37)*10,nPixels=64,bgAlbedo=0.6 $
         ,spots=spots[0:0],/scaling)

  for i=0,35 do tvscl, rtn[i].img , 64+(i MOD 9)*66, 64+(i/9)*66  ;;; display 36 images

  tots=rtn.dib
  dmags=make_array(n_elements(tots),value=0d)
  for i=0,n_elements(tots)-1 do dmags[i]=-2.5*alog10(tots[i]/median(tots,/even)) 

  y1=max(dmags)+(0.1*(max(dmags)-min(dmags)))
  y2=min(dmags)-(0.1*(max(dmags)-min(dmags)))
  !x.tickinterval = 90
  plot,rtn.obsWLongDeg,dmags,yr=[y1,y2],/ystyle $
      , XSTYLE=1 $
      , PSYM=-6 $
      , SYMSIZE=1.5 $
      , XTITLE='Sub-observer longitude, degW' $
      , YTITLE='Disk-Integrated Brightness, arbitrary units' $
      , COLOR=ishft(255L,16) $
      , /NOERASE, /NOCLIP  ;;; Over-plot lightcurve in blue
endelse


;;; To write a PNG of the plot:
;;;
;;;   write_png,'example.png',tvrd(true=1)

end
