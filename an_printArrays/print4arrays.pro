FUNCTION print4arrays, filename, x, y, z, za
nelemx = N_ELEMENTS(x)
nelemy = N_ELEMENTS(y)
nelemz = N_ELEMENTS(z)
nelemza = N_ELEMENTS(za)
nelem = MIN([nelemx,nelemy,nelemz,nelemza])
GET_LUN, lun
OPENW, lun, filename
FOR i=0, nelem -1 DO BEGIN
PRINTF, lun, x[i],',',y[i],',',z[i],',',za[i]
ENDFOR
CLOSE, lun
FREE_LUN, lun
RETURN, 0
END