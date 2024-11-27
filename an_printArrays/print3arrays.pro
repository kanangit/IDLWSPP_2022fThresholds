function print3arrays, filename, x, y, z
nelemx = N_ELEMENTS(x)
nelemy = N_ELEMENTS(y)
nelemz = N_ELEMENTS(z)
nelem = MIN([nelemx,nelemy,nelemz])
GET_LUN, lun
OPENW, lun, filename
FOR i=0, nelem -1 DO BEGIN
PRINTF, lun, x[i],',',y[i],',',z[i]
ENDFOR
CLOSE, lun
FREE_LUN, lun
RETURN, 0
END