;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Utility plotting function
function ilc_wrap_plot, Xs, Ys $
       , YTITLE $
       , dmerr $
       , COLOR=COLOR $
       , NOERASE=NOERASE $
       , NOCLIP=NOCLIP

  yMax = max(Ys,min=yMin)
  y2 = yMax*1.1 - yMin*.1
  y1 = yMin*1.1 - yMax*.1

  if !d.name eq 'PS' then begin
    charSize=1.4
    thick=4
    font=10
  endif else begin
    charSize=1.4
    thick=1.5
  endelse

  !x.tickinterval = 90

  plot,Xs,Ys,yr=[y1,y2] $
      , charsize=charSize $
      , charthick=thick $
      , thick=thick $
      , xthick=thick $
      , ythick=thick $
      , XSTYLE=1 $
      , /YSTYLE $
      , PSYM=-6 $
      , SYMSIZE=1.5 $
      , XTITLE='Sub-observer longitude, degW' $
      , YTITLE=YTITLE $
      , COLOR=COLOR $
      , NOERASE=NOERASE $
      , NOCLIP=NOCLIP

  thisIsGdl,isGdl
  if n_elements(dmerr) eq n_elements(Xs) and not keyword_set(isGdl) then begin
    catcherr = 0L
    catch,catcherr
    if catcherr eq 0L then begin
      plotyerr,Xs,Ys,dmerr,/hat,thick=thick
    endif else begin
      message,/continue,'Skipping PLOTYERR ...'
    endelse
    catch,/cancel
  end

  return, 0b
end
