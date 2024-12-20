PRO driver_an_2022fThresholds_pmap


  datapath = 'e:\ag\uiowa2020\oneDrive\bDocs\expAnalysisBackup\c_14226_vid56\20241208forP_2022fThresholds_pmap_corr\01_code_an_2022fThresholds_pmap\'

  curDate='20241208'
  print, curDate
  print, datapath
  coreName = STRCOMPRESS('particles_map_' + STRING(curDate) + 'ff_', /REMOVE_ALL)
  ;ROI
  leftBorder = 600.0d
  rightBorder= 1000.0d;
  yMin = 0.0d;
  yMax = 1200.0d
  
  screen_ratio =  DOUBLE(rightBorder - leftBorder) / DOUBLE(yMax - yMin)
  plot_hor_size = yMax - yMin

  ;start and end frames
  
  iBegin = 1;
  iEnd =  800

  ;start and end frames for pulse postition fitting:
  ;iBegin_ppulse = 513
  ;iEnd_ppulse = 614


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
    fnameCSV = STRCOMPRESS(coreName+string(myFrame,FORMAT='(I04)')+'.csv', /REMOVE_ALL)
    indMyFrame = WHERE(s_ROI.iFrame eq myFrame)
    
    curTitle = "particle map for frame"+string(myFrame)
    p = plot(s_ROI.X[indMyFrame], s_ROI.Y[indMyFrame], LINESTYLE = 'none',  $
      XRANGE = [yMin, yMax], YRANGE = [leftBorder, rightBorder], $
      SYMBOL = 'dot', ASPECT_RATIO = 1, SYM_SIZE = 6, SYM_FILLED = 1, $
      Margin = [0.05,0.05,0.05,0.10], $
       DIMENSIONS = [plot_hor_size, plot_hor_size * screen_ratio], $
       TITLE = curTitle, /CURRENT)
    p.save, fname, WIDTH = ROUND(2 * plot_hor_size), LENGTH = ROUND(2 * plot_hor_size * screen_ratio) 
    w = p.WINDOW
    w.erase
    retp = print2arrays(fnameCSV, s_ROI.X[indMyFrame], s_ROI.Y[indMyFrame])

  ENDFOR
  ;p.close

END