FUNCTION s_select, s, ind
nelem = N_ELEMENTS(s.X)
returnStructur = {iParticle:LONARR(nelements),iFrame:LONARR(nelements), $
  area:DBLARR(nelements), X:DBLARR(nelements),Y:DBLARR(nelements),error:DBLARR(nelements)}
returnStructur.iParticle = s.iParticle[ind]
returnStructur.iFrame = s.iFrame[ind]
returnStructur.area = s.area[ind]
returnStructur.X = s.X[ind]
returnStructur.Y = s.Y[ind] 
returnStructur.error = s.error[ind]
RETURN, returnStructur
END