FUNCTION get_polygon, N, R

nelements = N

  s = {X:DBLARR(nelements),Y:DBLARR(nelements)}

X = R * COS(DINDGEN(N) * 2.0d * !DPI / DOUBLE(N) + !DPI / 2.0d)
Y = R * SIN(DINDGEN(N) * 2.0d * !DPI / DOUBLE(N) + !DPI / 2.0d)
s.X = X
X = 0
s.Y = Y
Y = 0
  

  
  
  ;asdfg
  
  
  
  
  return, s

polygon

END