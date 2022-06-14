;2022-06-08
;based on driver_an_2022fThresholds_vrn. Plotting only one frame for Figure.

PRO driver_an_2022fThresholds_vrn_forFig
  ;

  datapath = 'e:\OneDrive - University of Iowa\bDocs\expAnalysisBackup\c_14226_vid56\20220614forP_2022fThresholds_fig\03_code_an_2022fThresholds_vrn_forFig\'

  curDate='20220614'
  print, curDate
  print, datapath
  coreName = STRCOMPRESS('voronoiMap_' + STRING(curDate) + 'ff_', /REMOVE_ALL)
  ;ROI
  leftBorder = 600.0d
  rightBorder= 1000.0d;
  yMin = 0.0d;
  yMax = 1000.0d

  ;frame of interest
  myFrame = 326


  ;start and end frames for pulse postition fitting:
  iBegin_ppulse = 521
  iEnd_ppulse = 624


  CD, datapath
  CD, 'inputs'

  filenam = DIALOG_PICKFILE(/READ, FILTER = '*.sav')
  RESTORE, filenam

  s_pulsePos = read_pulse_pos()
  indPpulseFilt = WHERE(s_pulsePos.time GE iBegin_ppulse AND s_pulsePos.time LE iEnd_ppulse)
  coeffs = POLY_FIT(s_pulsePos.time[indPpulseFilt],s_pulsePos.position[indPpulseFilt],1,/DOUBLE)

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





    fnam = STRCOMPRESS(coreName+string(myFrame,FORMAT='(I04)'), /REMOVE_ALL)
    indMyFrame = WHERE(s_ROI.iFrame eq myFrame)
    pulsePos = coeffs[0] + myFrame * coeffs[1]
    dratio = oldplot_defects_postScript(s_ROI.X[indMyFrame], s_ROI.Y[indMyFrame], pulsePos, fnam)


END