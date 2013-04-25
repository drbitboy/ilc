pro thisIsGdl,isGdl
common thisIsGdl_common, saveIsGdl
  if n_elements(saveIsGdl) ne 1L then begin
    saveIsGdl = 0b
    catcherr = 0L
    catch,catcherr
    if catcherr eq 0L then if keyword_set(!gdl) then saveIsGdl = 1b
    catch,/cancel
  endif
  isGdl = saveIsGdl
  return
end
