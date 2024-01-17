PRO driver_an_2022fThresholds_pmap
;test

  datapath = 'C:\Users\kanton\OneDrive - University of Iowa\bDocs\expAnalysisBackup\c_14226_vid56\20230628forP_2022fThresholds_pmap\01_an_2022fThresholds_pmap\'

  curDate='20230628'
  print, curDate
  print, datapath
  coreName = STRCOMPRESS('particles_map_' + STRING(curDate) + 'ff_', /REMOVE_ALL)
  ;ROI
  leftBorder = 600.0d
  rightBorder= 1000.0d;
  yMin = 0.0d;
  yMax = 1000.0d

  ;start and end frames
  iBegin = 1
  iEnd =  660


  ;start and end frames for pulse postition fitting:
  iBegin_ppulse = 513
  iEnd_ppulse = 614


  CD, datapath
  CD, 'inputs'

  filenam = DIALOG_PICKFILE(/READ, FILTER = '*.sav')
  RESTORE, filenam


  CD, '..'
  FILE_MKDIR, 'outputs'
  CD, 'outputs'


  ;exclude all the bad elements of the data:
  indGood = WHERE(FINITE(s.iparticle) $
    AND FINITE(s.area) AND FINITE(s.X) AND FINITE(s.Y) $
    AND FINITE(s.error))
  s_good = s_select(s, indGood)



  ;select the region of interest:
  indROI = WHERE(s_good.X LE rightBorder AND s_good.X GE leftBorder AND s_good.Y LE ymax AND s_good.Y GE yMin)
  s_ROI = s_select(s_good,indROI)
  y_temp = yMax - (s_ROI.Y - yMin) ;because the vertical screen coordinates
  ;are from top to bottom, we make this change of variables
  ;switch x and y coordinates:
  s_ROI.Y = s_ROI.X
  s_ROI.X = y_temp


  myFrame = 1;

  FOR myFrame = iBegin, iEnd DO BEGIN
    fname = STRCOMPRESS(coreName+string(myFrame,FORMAT='(I04)')+'.tif', /REMOVE_ALL)
    indMyFrame = WHERE(s_ROI.iFrame eq myFrame)
    
    
    p = plot(s_ROI.X[indMyFrame], s_ROI.Y[indMyFrame], LINESTYLE = 'none',  $
      SYMBOL = 'dot', ASPECT_RATIO = 1, SYM_SIZE = 3, SYM_FILLED = 1)

    stop

  ENDFOR

END