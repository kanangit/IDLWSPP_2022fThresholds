;v2.1
function print2arrays, filename, x, y, firstinteger = firstinteger
nelemx = N_ELEMENTS(x)
nelemy = N_ELEMENTS(y)
nelem = MIN([nelemx,nelemy])
comma = ','
GET_LUN, lun
OPENW, lun, filename
if (KEYWORD_SET(firstinteger)) then begin
  FOR i=0, nelem -1 DO BEGIN
    PRINTF, lun,format = '(I18,A1,I18)',x[i],comma,y[i]
  ENDFOR
endif  else begin
  FOR i=0, nelem -1 DO BEGIN
    PRINTF, lun,format = '(E18.10,A1,E18.10)',x[i],comma,y[i]
  ENDFOR
endelse
CLOSE, lun
FREE_LUN, lun
RETURN, 0
END