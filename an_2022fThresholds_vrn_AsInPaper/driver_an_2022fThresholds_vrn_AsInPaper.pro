PRO driver_an_2022fThresholds_vrn_AsInPaper

  datapath = 'e:\ag\uiowa2020\oneDrive\bDocs\expAnalysisBackup\c_14226_vid59\20241207forP_2022fThresholds_vrn_corr\01_code_an_2022fThresholds_vrn_AsInPaper\'

  curDate='20241207'
  print, curDate
  print, datapath
  coreName = STRCOMPRESS('voronoiMap_' + STRING(curDate) + 'ff_', /REMOVE_ALL)
  ;ROI
  leftBorder = 600.0d
  rightBorder= 1000.0d;
  yMin = 0.0d;
  yMax = 1200.0d

  ;start and end frames
  iBegin = 1; 1st submission: iBegin = 500
  iEnd =  1000

  ;start and end frames for pulse postition fitting:
  iBegin_ppulse = 661
  iEnd_ppulse = 777

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


  myFrame = 1;

  FOR myFrame = iBegin, iEnd DO BEGIN
    fname = STRCOMPRESS(coreName+string(myFrame,FORMAT='(I04)')+'.tif', /REMOVE_ALL)
    indMyFrame = WHERE(s_ROI.iFrame eq myFrame)
    pulsePos = coeffs[0] + myFrame * coeffs[1]
    dratio = oldplot_defects_asip(s_ROI.X[indMyFrame], s_ROI.Y[indMyFrame], pulsePos, fname)

  ENDFOR
END