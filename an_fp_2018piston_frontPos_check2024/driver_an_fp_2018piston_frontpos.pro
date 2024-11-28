;+
;v. 0.01. 2019.01.15
;;v. 0.02. 2019.01.15 Debug and all
; :Description:
; takes the output of the an_trimandsave and finds front positions vs frame number
;    ususes a lot of the code from driver_transfshockframerefer2ndderiv.pro
;    v. 0.17
;
;
;
;
;
;
; :Author: Anton Kananovich
;-
PRO driver_an_fp_2018piston_frontPos

  datapath = 'C:\Users\siomau\OneDrive - University of Iowa\bDocs\expAnalysisBackup\c_14226_vid56\20241127forP_2022fThresholds_check_PulsePos\'
  
  curDate='20241127'
  
  
  ;define borders of the region of interest
  binWidth_in_b = 1.0d ; bin width expressed in interpartile distances
  leftBorder = 600.0d
  rightBorder= 1000.0d;
  yMin = 0.0d;
  yMax = 1200.0d ;start and end frames
  plot_hor_size = yMax - yMin ; addition 2024 compatibility
  screen_ratio = 1.0d;
  
  iBegin = 415
  iEnd =  584
  n0 = 0.0011926099d ; equillibrium density (perticles per pixel squared)
  frameRate = 870.0d ;frames per second
  cam_resol = 0.0254303 ;mm/px
  coreName = STRCOMPRESS('frontP_ff'+STRING(iBegin)+'-' + STRING(iEnd) + '_' + STRING(curDate))
  
  aWignerZeiss = 16.337134d
  b_interparticle = aWignerZeiss * SQRT((2.0d*!DPI/SQRT(3.0d)))
  
  ;define the criteria for selecting the interval of the data for
  ;linear fit (for determination of the pulse velocity):
  criterion_amplitude = 2.50d
  criterion_edge = b_interparticle * 3.0
  
  
  CD, datapath
  CD, 'inputs'
  
  dy = binWidth_in_b * b_interparticle
  nB = FLOOR((yMax - yMin)/dY); number of bins
  yBins = DINDGEN(nB)*dY+yMin + dY/2
  yBins = DINDGEN(nB)*dY+yMin + dY/2
  
  filenam = DIALOG_PICKFILE(/READ, FILTER = '*.sav')
  RESTORE, filenam
  
  CD, '..'
  FILE_MKDIR, 'outputs'
  CD, 'outputs'
  
  ;we need only the data inside the region of interest
  ind = WHERE(s.X LE rightBorder AND s.X GE leftBorder $
    AND s.Y LE ymax AND s.Y GE yMin)
  arrlen = N_ELEMENTS(s.X[ind]);
  Yroi = yMax - (s.Y[ind] - yMin) ;because the vertical screen coordinates
  ;are from top to bottom, we make this change of variables
  time = s.iframe[ind] ;array to track time
  
   s_ROI = s_select(s,ind) ; 2024 compatibility
   
   y_temp = yMax - (s_ROI.Y - yMin) ;because the vertical screen coordinates
   ;are from top to bottom, we make this change of variables
   ;switch x and y coordinates:
   s_ROI.Y = s_ROI.X
   s_ROI.X = y_temp
  
  pulseT = DBLARR(iEnd-iBegin+1);
  pulsePos = DBLARR(iEnd-iBegin+1);
  pulseAmp = DBLARR(iEnd-iBegin+1);
  
  
  FOR i = iBegin, iEnd DO BEGIN
    print, 'frame ', i
    myFrame = i ; compatibility 2024
    ;building histogram using cloud-in-cell (cic):
    indf = WHERE(time EQ i)
    Y = Yroi[indf]
    YforHist = (Y - yMin)/(yMax-yMin)*DOUBLE(nB)
    weighNumDens = DBLARR(N_ELEMENTS(YforHist))+1 ; weights for the number density
    histNumDens = CIC(weighNumDens, YforHist, nB, /ISOLATED) / dy / (rightBorder - leftBorder)
    pAmp = histNumDens / n0
    
    maxDen = MAX(pAmp, ind_maxDen)
    pos_maxDen = ybins(ind_maxDen)
    pulseT[i-iBegin] = DOUBLE(i)
    pulsePos[i-iBegin] = pos_maxDen
    pulseAmp[i-iBegin] = maxDen
    
    
    fnameTif = STRCOMPRESS(coreName+ "_histo"+ "_f"+string(myFrame,FORMAT='(I04)')+'.tif', /REMOVE_ALL)
    fnameTifp = STRCOMPRESS(coreName+ "_pmap" +"_f"+string(myFrame,FORMAT='(I04)')+'.tif', /REMOVE_ALL)
    fnameCSV = STRCOMPRESS(coreName+"_f"+string(myFrame,FORMAT='(I04)')+'.csv', /REMOVE_ALL)
   

    curTitle = "histogram for frame"+string(myFrame)
;    p = plot(ybins, pAmp, LINESTYLE = 'none',  $
;      XRANGE = [0, 1200], YRANGE = [0, 1], $
;      SYMBOL = 'dot', ASPECT_RATIO = 1, SYM_SIZE = 6, SYM_FILLED = 1, $
;      Margin = [0.05,0.05,0.05,0.10], $
;      DIMENSIONS = [plot_hor_size, plot_hor_size * screen_ratio], $
;      TITLE = curTitle, /CURRENT)


p = plot(ybins, pAmp, LINESTYLE = 'none',  $
  XRANGE = [0, 1200], YRANGE = [0, 3.1], $
  SYMBOL = 'dot', SYM_SIZE = 6, SYM_FILLED = 1, $
  TITLE = curTitle, /CURRENT)

    p.save, fnameTifp, WIDTH = ROUND(2 * plot_hor_size), LENGTH = ROUND(2 * plot_hor_size * screen_ratio)
    w = p.WINDOW
    w.erase
    retp = print2arrays(fnameCSV, ybins, pAmp)
    
   curTitle = "particles map for frame"+string(myFrame) 
    
    pp = plot(s_ROI.X[indf], s_ROI.Y[indf], LINESTYLE = 'none',  $
      XRANGE = [yMin, yMax], YRANGE = [leftBorder, rightBorder], $
      SYMBOL = 'dot', ASPECT_RATIO = 1, SYM_SIZE = 6, SYM_FILLED = 1, $
      Margin = [0.05,0.05,0.05,0.10], $
      DIMENSIONS = [plot_hor_size, plot_hor_size * screen_ratio], $
      TITLE = curTitle, /CURRENT)
    pp.save, fnameTifp, WIDTH = ROUND(2 * plot_hor_size), LENGTH = ROUND(2 * plot_hor_size * screen_ratio)
    wp = pp.WINDOW
    wp.erase
    retp = print2arrays(fnameCSV, s_ROI.X[indf], s_ROI.Y[indf])
    
    
  ENDFOR
  ;saving the time, max pulse amplitude position and max amplitude
  ;vs time:
  seconds = STRING(SYSTIME(/seconds),FORMAT='(I18)')
  fname = STRCOMPRESS(corename + '_BinWidth' + $
    STRING(binWidth_in_b, format = '(D3.1)') + 'b_' + seconds + '.csv', $
    /REMOVE_ALL)
  fnameForFit = STRCOMPRESS(corename + '_BinWidth' + $
    STRING(binWidth_in_b, format = '(D3.1)') + 'b_' + seconds + '_forFit.csv', $
    /REMOVE_ALL)
  head = ['time','pulsePos', 'amplitude']
  thead = [STRCOMPRESS('time is in frameNumbers. PulsePos is in pixels. Amplitude is in normalized with equilibrium number density')]
  WRITE_CSV, fname, pulseT, pulsePos, pulseAmp, HEADER = head, TABLE_HEADER= thead
  
  ;selecting the data range accoring to the crierion stated in the
  ;beginning of the code for the linear fit:
  
  ind_iwant = WHERE(pulsePos LE yMax - criterion_edge AND $
    pulseAmp GE criterion_amplitude)
  if (N_ELEMENTS(ind_iwant) GT 1) then begin
    pulseT_forFit = pulseT[ind_iwant]
    pulsePos_forFit = pulsePos[ind_iwant]
    pulseAmp_forFit = pulseAmp[ind_iwant]
    
    ;just for a check, calculate the pulse speed we obtain with this data
    coeffs = POLY_FIT(pulseT_forFit,pulsePos_forFit,1,/DOUBLE)
    print, 'pulse velocity = ', coeffs[1] * frameRate * cam_resol
    WRITE_CSV, fnameForFit, pulseT_forFit, pulsePos_forFit, pulseAmp_forFit, HEADER = head, TABLE_HEADER= thead
    
  endif
  
  
END
