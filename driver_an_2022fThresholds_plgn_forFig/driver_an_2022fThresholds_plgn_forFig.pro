pro driver_an_2022fThresholds_plgn_forFig

  datapath = 'e:\ag\uiowa2020\oneDrive\bDocs\expAnalysisBackup\c_14226_vid60\20241206forP_2022fThresholds_plgn_fig2_corr\08_an_2022fThresholds_plgn_forFig\'


  ;start and end frames for pulse postition fitting:

  curDate='20241206'
  print, curDate
  print, datapath
  coreName = STRCOMPRESS('pgnConstr_' + STRING(curDate) + 'ff_', /REMOVE_ALL)
  ;ROI
  leftBorder = 600.0d
  rightBorder= 1000.0d;
  yMin = 0.0d;
  yMax = 1200.0d


  ;frame of interest
  myFrame = 775

  ;start and end frames for pulse postition fitting:
  iBegin_ppulse = 784
  iEnd_ppulse = 935



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




  plot_ybegin = leftBorder
  plot_yend = rightBorder
  plot_xbegin = 0.0d
  plot_xend = yMax


  preshock_offset = 50.0
  postshock_offset = 50.0

  textHeight = 525
  shockPos = coeffs[0] + myFrame * coeffs[1]

  shockPos_indiv = 762.345034846977d

  postshock_right_border = shockPos - preshock_offset
  preshock_left_border = shockPos + postshock_offset

  shockfront_x = [shockPos, shockPos]
  shockFront_y = [textHeight, plot_yend]

  shockfront_indiv_x = [shockPos_indiv, shockPos_indiv]
  shockFront_indiv_y = [textHeight, plot_yend]

  postshock_right_x = [postshock_right_border, postshock_right_border]
  postshock_right_y = [textHeight, plot_yend]
  preshock_left_x = [preshock_left_border, preshock_left_border]
  preshock_left_y = [textHeight, plot_yend]


  screenWidth = 1600
  ratio = DOUBLE(plot_yend - plot_ybegin)/DOUBLE(plot_xend - plot_xbegin)


  ;INPUT5
  angleThreshold = 75.0


  ;INPUT7
  fname_defect_stat = datapath+'\outputs\Defect_stat.txt'

  ;INPUT8
  fname_vertex_stat = datapath+'\outputs\Vertex_stat.txt'


  close,3 ;File unit 3 will be assigned to 'fname_defect_stat'
  close,4 ;File unit 4 will be assigned to 'fname_vertex_stat'


  ;The following 'OPENW' procedure is peculiar to IDL.
  ;OPENW opens a new file for output. If the file exists, its old contents are destroyed.
  ;It is used here to open 'fname_defect_stat' for output.
  openw,3,fname_defect_stat ;File unit 3 is assigned to 'Defect_stat.txt'

  openw,4,fname_vertex_stat ;File unit 4 is assigned to 'Vertex_stat.txt'


  frameNumber = myFrame



  indMyFrame = where(s_ROI.iFrame EQ frameNumber)

  Xtemp= s_ROI.X[indMyFrame] ;X coordinates of particle positions in the current frame.
  Ytemp= s_ROI.Y[indMyFrame] ;Y coordinates of particle positions in the current frame.


  numberOfParticles = size(Xtemp,/N_ELEMENTS) ;Total number of particle coordinates in the current frame.


  ;convert px to mm, shift the origin to zero:
  cam_resol = 0.0254303 ;mm/px
  Xtemp_mm =  px_to_mm(Xtemp, plot_xbegin, cam_resol)
  Ytemp_mm =  px_to_mm(Ytemp, plot_ybegin, cam_resol)
  plot_xbegin_mm = px_to_mm(plot_xbegin, plot_xbegin, cam_resol)
  plot_xend_mm = px_to_mm(plot_xend, plot_xbegin, cam_resol)
  plot_ybegin_mm = px_to_mm(plot_ybegin, plot_ybegin, cam_resol)
  plot_yend_mm = px_to_mm(plot_yend, plot_ybegin, cam_resol)
  shockfront_x_mm = px_to_mm(shockfront_x, plot_xbegin, cam_resol)
  shockfront_y_mm = px_to_mm(shockfront_y, plot_ybegin, cam_resol)
  shockfront_indiv_x_mm = px_to_mm(shockfront_indiv_x, plot_xbegin, cam_resol)
  shockfront_indiv_y_mm = px_to_mm(shockfront_indiv_y, plot_ybegin, cam_resol)


  X = fltarr(numberOfParticles,1)
  Y = fltarr(numberOfParticles,1)

  for particle=0,numberOfParticles-1 do begin
    X(particle,0) = Xtemp_mm[particle]
    Y(particle,0) = Ytemp_mm[particle]
  endfor


  ;The following procedures 'DEVICE', 'LOADCT' and 'WINDOW' are used to set up graphic display options.
  originalDevice = !d.name
  set_plot, 'WIN'
  device, retain=2, decomposed =0
  loadct, 39
  !p.color = 0
  !p.background = 255
  window, 2, xsize = screenWidth, ysize =screenWidth*ratio


  plot,X,Y,psym=3, isotropic=1, xrange = [plot_xbegin_mm,plot_xend_mm], yrange = [plot_ybegin_mm,plot_yend_mm], $
    xstyle = 1, ystyle=1, charsize = 2.0, thick = 8.0, charthick = 2, /NODATA, $
    title = 'polygon construction for frame' + string(frameNumber)


  allAngles = 0
  allSides  = 0
  numTriangles = 0
  sideSelector,X,Y,allAngles,allSides,numTriangles

  track = 0
  maxAngleSelect,allAngles,allSides,numTriangles,track

  sidesToPlot = 0
  selectAndPlotSides,track,allSides,allAngles,angleThreshold,sidesToPlot

  Ulater=0
  polygonTracker=0
  defectClassification=0
  unt=0
  polygonIdentification,sidesToPlot,X,Y,numberOfParticles,Ulater,polygonTracker,defectClassification,unt

  start = 2
  for nu=2,5 do begin
    numsides=nu
    polygonIdentification2,polygonTracker,defectClassification,numsides,start,unt
    start = start+2
  endfor

  out = 0
  TriangleCountandPlot,defectClassification,unt,out
  numTriangles = out

  QuadrilateralCountandPlot,defectClassification,unt,out
  numQuadrilaterals = out

  PentagonsCountandPlot,defectClassification,unt,out
  numPentagons = out

  HexagonCountandPlot,defectClassification,unt,out
  numHexagons = out

  Defect_numbers = [frameNumber,numberOfParticles, numTriangles, numQuadrilaterals, numPentagons,numHexagons]
  printf,3,Defect_numbers

  oplot, X,Y, psym = 3, symsize = 1

  ;plotting the shock front position:
  plots, shockfront_x_mm, shockFront_y_mm, THICK = 3
  ;plots, postshock_right_x, postshock_right_y, COLOR = 47, THICK = 2
  ;plots, preshock_left_x, preshock_left_y, COLOR = 254, THICK = 2


 ; plots, shockfront_indiv_x_mm, shockfront_indiv_y_mm, COLOR = 60, THICK = 2

  scaled = TVRD()
  TVLCT, r, g, b, /Get
  s = Size(scaled, /Dimensions)
  image24 = BytArr(3, s[0], s[1])
  image24[0,*,*] = r[scaled]
  image24[1,*,*] = g[scaled]
  image24[2,*,*] = b[scaled]
  image24 = TVRD(True=1)
  image24 = Reverse(image24, 3)

  ;determine system time, to make a filename unique:
  seconds = STRING(SYSTIME(/seconds),FORMAT='(I18)')

  filename = datapath+ '\outputs\images\' + STRCOMPRESS(coreName+string(frameNumber,FORMAT='(I04)'), /REMOVE_ALL)

  filename_eps = filename+ STRCOMPRESS('_' + seconds + '.eps', /REMOVE_ALL)
  filename_pmap_eps = filename+ '_' + seconds +'_pmapOnly.eps'
  filename_tif = filename+STRCOMPRESS('_'+seconds+'.tif', /REMOVE_ALL)
  filename_emf = filename+ '_' + seconds  +'.emf'



  Write_Tiff,  filename_tif, image24, 1

  VertexClassification,defectClassification,numberOfParticles,Ulater,frameNumber,out
  vertexStatnew = out
  printf,4,vertexStatnew,FORMAT='(F7.1,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1)'
  close,3
  close,4

  set_plot, originalDevice

  ;-----------------------------------------------------------------------------------------------------
  ;-----------------------------------------------------------------------------------------------------
  originalDevice = !d.name
  set_plot, 'ps'
  device, filename=filename_eps
  device, xsize= 1.5 * (3 + 3/8), ysize=  ratio * 1.5 * ((3 + 3/8)), /inches
  device, FONT_SIZE = 12, /TIMES
  device, color=1, bits_per_pixel=24
  device, /encapsulated


  plot,X,Y,psym=3, isotropic=1, xrange = [plot_xbegin_mm,plot_xend_mm], yrange = [plot_ybegin_mm,plot_yend_mm], $
    xstyle = 1, ystyle=1, charsize = 0.8, thick = 8.0, charthick = 1, /NODATA, $
    title = 'polygon construction for frame' + string(frameNumber)


  allAngles = 0
  allSides  = 0
  numTriangles = 0
  sideSelector,X,Y,allAngles,allSides,numTriangles

  track = 0
  maxAngleSelect,allAngles,allSides,numTriangles,track

  sidesToPlot = 0
  selectAndPlotSides,track,allSides,allAngles,angleThreshold,sidesToPlot

  Ulater=0
  polygonTracker=0
  defectClassification=0
  unt=0
  polygonIdentification,sidesToPlot,X,Y,numberOfParticles,Ulater,polygonTracker,defectClassification,unt

  start = 2
  for nu=2,5 do begin
    numsides=nu
    polygonIdentification2,polygonTracker,defectClassification,numsides,start,unt
    start = start+2
  endfor

  out = 0
  TriangleCountandPlot,defectClassification,unt,out
  numTriangles = out

  QuadrilateralCountandPlot,defectClassification,unt,out
  numQuadrilaterals = out

  PentagonsCountandPlot,defectClassification,unt,out
  numPentagons = out

  HexagonCountandPlot,defectClassification,unt,out
  numHexagons = out

  Defect_numbers = [frameNumber,numberOfParticles, numTriangles, numQuadrilaterals, numPentagons,numHexagons]


  oplot, X,Y, psym = 3, symsize = 1

  ;plotting the shock front position:
  plots, shockfront_x_mm, shockFront_y_mm, THICK = 3
  ;plots, postshock_right_x, postshock_right_y, COLOR = 47, THICK = 2
  ;plots, preshock_left_x, preshock_left_y, COLOR = 254, THICK = 2

  ;plots, shockfront_indiv_x_mm, shockfront_indiv_y_mm, COLOR = 60, THICK = 2

  device, /close_file
  set_plot, originalDevice

end

