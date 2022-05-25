PRO driver_an_2022fThresholds_vrn

  datapath = 'e:\OneDrive - University of Iowa\bDocs\expAnalysisBackup\c_14226_vid56\20220524forP_2022fThresholds\01_code_an_2022fThresholds_vrn\'

  curDate='20220525'
  print, curDate
  print, datapath
  
  ;ROI
  leftBorder = 600.0d
  rightBorder= 1000.0d;
  yMin = 0.0d;
  yMax = 1200.0d
  
  ;start and end frames  
  iBegin = 1
  iEnd =  660
  
  
  CD, datapath
  CD, 'inputs'
  
  filenam = DIALOG_PICKFILE(/READ, FILTER = '*.sav')
  RESTORE, filenam
  
  
  stop;
  CD, '..'
  FILE_MKDIR, 'outputs'
  CD, 'outputs'
  stop;
END