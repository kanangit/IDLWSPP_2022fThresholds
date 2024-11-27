FUNCTION print6arrays, filename, x, y, z, za, zb, zc, firstinteger = firstinteger
nelemx = N_ELEMENTS(x)
nelemy = N_ELEMENTS(y)
nelemz = N_ELEMENTS(z)
nelemza = N_ELEMENTS(za)
nelemzb = N_ELEMENTS(zb)
nelemzc = N_ELEMENTS(zc)
nelem = MIN([nelemx,nelemy,nelemz,nelemza,nelemzb,nelemzc])
comma = ','
GET_LUN, lun
OPENW, lun, filename
if (KEYWORD_SET(firstinteger)) then begin
  FOR i=0, nelem -1 DO BEGIN
    PRINTF, lun,format = '(I18,A1,I18,A1,E18.10,A1,E18.10,A1,E18.10,A1,E18.10)',x[i],comma,y[i],comma,z[i],comma,za[i],comma,zb[i],comma,zc[i]
  ENDFOR
endif  else begin
  FOR i=0, nelem -1 DO BEGIN
    PRINTF, lun,format = '(E18.10,A1,E18.10,A1,E18.10,A1,E18.10,A1,E18.10,A1,E18.10)',x[i],comma,y[i],comma,z[i],comma,za[i],comma,zb[i],comma,zc[i]
  ENDFOR
endelse
CLOSE, lun
FREE_LUN, lun
RETURN, 0
END