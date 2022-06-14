pro driver_an_2022fThresholds_plgn_forFig

  datapath = 'e:\OneDrive - University of Iowa\bDocs\expAnalysisBackup\c_14226_vid56\20220614forP_2022fThresholds_fig\05_code_an_2022fThresholds_plgn_forFig_f544\'


  ;start and end frames for pulse postition fitting:
  iBegin_ppulse = 1469
  iEnd_ppulse = 2093

  curDate='20220607'
  print, curDate
  print, datapath
  coreName = STRCOMPRESS('pgnConstr_' + STRING(curDate) + 'ff_', /REMOVE_ALL)
  ;ROI
  leftBorder = 600.0d
  rightBorder= 1000.0d;
  yMin = 0.0d;
  yMax = 1000.0d


  ;frame of interest
  myFrame = 544

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



  forceXlen = 1100
  forceYlen = 1200
  screenWidth = 1200

  plot_ybegin = 600;
  plot_yend = 1000;
  plot_xbegin = 0
  plot_xend = 1100
  screenWidth = forceXlen
  ratio = DOUBLE(plot_yend - plot_ybegin)/DOUBLE(plot_xend - plot_xbegin)


  startFrame = iBegin
  endFrame = iEnd

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


  for frameNumber=startFrame,endFrame do begin

    preshock_offset = 50.0
    postshock_offset = 50.0

    textHeight = 525
    shockPos = coeffs[0] + frameNumber * coeffs[1]
    postshock_right_border = shockPos - preshock_offset
    preshock_left_border = shockPos + postshock_offset

    shockfront_x = [shockPos, shockPos]
    shockFront_y = [textHeight, plot_yend]

    postshock_right_x = [postshock_right_border, postshock_right_border]
    postshock_right_y = [textHeight, plot_yend]
    preshock_left_x = [preshock_left_border, preshock_left_border]
    preshock_left_y = [textHeight, plot_yend]

    indMyFrame = where(s_ROI.iFrame EQ frameNumber)

    Xtemp= s_ROI.X[indMyFrame] ;X coordinates of particle positions in the current frame.
    Ytemp= s_ROI.Y[indMyFrame] ;Y coordinates of particle positions in the current frame.


    numberOfParticles = size(Xtemp,/N_ELEMENTS) ;Total number of particle coordinates in the current frame.

    X = fltarr(numberOfParticles,1)
    Y = fltarr(numberOfParticles,1)

    for particle=0,numberOfParticles-1 do begin
      X(particle,0) = Xtemp[particle]
      Y(particle,0) = Ytemp[particle]
    endfor

    xlen = forceXlen
    ylen = forceYlen


    ;The following procedures 'DEVICE', 'LOADCT' and 'WINDOW' are used to set up graphic display options.
    device, retain=2, decomposed =0
    loadct, 39
    !p.color = 0
    !p.background = 255
    window, 2, xsize = screenWidth, ysize =screenWidth*ratio

 
    plot,X,Y,psym=3, isotropic=1, xrange = [plot_xbegin,plot_xend], yrange = [plot_ybegin,plot_yend], $
      xstyle = 1, ystyle=1, charsize = 2.0, thick = 8.0, charthick = 2, /NODATA, $
      title = 'polygon construction for frame' + string(frameNumber)
    

    allAngles = 0
    allSides  = 0
    numTriangles = 0
    sideSelector,X,Y,allAngles,allSides,numTriangles

    ;User defined procedure 2: 'maxAngleSelect'
    ;
    ;This procedure is used below to perform the following tasks.
    ;
    ;I. Identify duplicate entries of bonds from the matrix 'allSides' and track
    ;   them with the use of a new matrix 'track'
    ;II.Update matrix 'allAngles' so that it contains value of the maximum angle opposite to a given bond.
    ;
    ;The arguments for this procedure are as follows.
    ;
    ; 1.allAngles(matrix): Updated matrix from the procedure 'sideSelector'.
    ;
    ; 2.allSides(matrix): Updated matrix from the procedure 'sideSelector'.
    ;
    ; 3.numTriangles(variable): Updated variable from the procedure 'sideSelector'.
    ;
    ; 4.track(matrix): On exit, contains values 0 or -5
    ;   value of 0 in row 3 (of matrix 'track') for example implies that a bond in row 3 of the matrix 'allSides'
    ;   is found for the first time.
    ;   value of -5 in row 4 (of matrix 'track') for example implies that a triangle side in row 4 has a duplicate entry in
    ;   one of the rows 1-3 of the matrix 'allSides'.
    ;
    ; 5.allAngles(matrix): On exit, contains the maximum angle opposite to a given bond.
    ;

    track = 0
    maxAngleSelect,allAngles,allSides,numTriangles,track

    ;User defined procedure 3: 'selectAngPlotSides'
    ;
    ;This procedure is used below to select the bonds that are opposite angles less than or equal to the threshold value, and plot them.
    ;
    ;The arguments for this procedure are as follows.
    ;
    ;1.track(matrix): Updated matrix from the procedure 'maxAngleSelect'.
    ;
    ;2.allSides(matrix): Updated matrix from the procedure 'maxAngleSelect'.
    ;
    ;3.allAngles(matrix): Updated matrix from the procedure 'maxAngleSelect'.
    ;
    ;4.angleThreshold(variable): INPUT5
    ;
    ;5.sidesToPlot(matrix): On exit, contains the X and Y coordinates of the end points of a
    ;  given bond that has not been removed. The X and Y coordinates for this end points are stored here similar to the
    ;  format in matrix 'allSides'.

    sidesToPlot = 0
    selectAndPlotSides,track,allSides,allAngles,angleThreshold,sidesToPlot


    ;User defined procedure 4: 'polygonIdentification'
    ;
    ;It is used below to identify bonds attached to a given particle position and arrange them according to the ascending
    ;order of angles they make with the positive X axis.
    ;
    ;The arguments for this procedure are as follows.
    ;
    ;1.sidesToPlot(matrix): Updated matrix from the procedure 'selectAndPlotSides'.
    ;
    ;2.X(matrix): X coordinates of the particle positions in the current frame.
    ;
    ;3.Y(matrix: Y coordinates of the particle positions in the current frame.
    ;
    ;4.numberOfParticles(variable): Total number of particle positions in the current frame.
    ;
    ;5.Ulater(matrix): On exit, contains the  number of bonds attached to a given particle position.
    ;
    ;6.polygonTracker(matrix): A temporary matrix.
    ;
    ;7.defectClassification(matrix): On exit, contains the end point coordinates of bonds attached to a given particle position according to
    ; the ascending order of angles the bonds make with the positive X axis.
    ;
    ;8.unt(variable): On exit, contains the number of rows in the matrix 'defectClassification' which has bond information.

    Ulater=0
    polygonTracker=0
    defectClassification=0
    unt=0
    polygonIdentification,sidesToPlot,X,Y,numberOfParticles,Ulater,polygonTracker,defectClassification,unt


    ;User defined procedure 5: 'polygonIdentification2'
    ;
    ;It is used below to identify polygons (triangles,quadrilaterals,pentagons and hexagons).
    ;
    ;The arguments for this procedure are as follows.
    ;
    ;1.polygonTracker(matrix): Updated matrix from the procedure 'polygonIdnetification'.
    ;
    ;2.defectClassification(matrix): Updated matrix from the procedure 'polygonIdnetification'.
    ;
    ;3.numsides(variable): Number of sides of the polygon - 1.
    ;
    ;4.start(variable): First column number of matrix 'defectClassification' to start searching for a particular polygon.
    ;
    ;5.unt(variable): Updated matrix from the procedure 'polygonIdnetification'.

    start = 2
    for nu=2,5 do begin
      numsides=nu
      polygonIdentification2,polygonTracker,defectClassification,numsides,start,unt
      start = start+2
    endfor


    ;User defined procedure 6: 'TriangleCountandPlot'
    ;
    ;It is used below to count and plot triangles.
    ;
    ;The arguments for this procedure are as follows.
    ;
    ;1.defectClassification(matrix): Updated matrix from procedure 'polygonIdentification2'.
    ;
    ;2.unt(variable): Updated variable from procedure 'polygonIdentification2'.
    ;
    ;3.out(variable): On exit, contains the total number of triangles in the polygon construction for current frame.

    out = 0
    TriangleCountandPlot,defectClassification,unt,out
    numTriangles = out

    ;User defined procedure 7: 'QuadrilateralCountandPlot'
    ;
    ;It is used below to count and plot quadrilaterals.
    ;
    ;The arguments for this procedure are as follows.
    ;
    ;1.defectClassification(matrix): Updated matrix from procedure 'polygonIdentification2'.
    ;
    ;2.unt(variable): Updated variable from procedure 'polygonIdentification2'.
    ;
    ;3.out(variable): On exit, contains the total number of quadrilaterals in the polygon construction for current frame.

    QuadrilateralCountandPlot,defectClassification,unt,out
    numQuadrilaterals = out

    ;User defined procedure 8: 'PentagonCountandPlot'
    ;
    ;It is used below to count and plot pentagons.
    ;
    ;The arguments for this procedure are as follows.
    ;
    ;1.defectClassification(matrix): Updated matrix from procedure 'polygonIdentification2'.
    ;
    ;2.unt(variable): Updated variable from procedure 'polygonIdentification2'.
    ;
    ;3.out(variable): On exit, contains the total number of pentagons in the polygon construction for current frame.


    PentagonsCountandPlot,defectClassification,unt,out
    numPentagons = out


    ;User defined procedure 9: 'HexagonCountandPlot'
    ;
    ;It is used below to count and plot hexagons.
    ;
    ;The arguments for this procedure are as follows.
    ;
    ;1.defectClassification(matrix): Updated matrix from procedure 'polygonIdentification2'.
    ;
    ;2.unt(variable): Updated variable from procedure 'polygonIdentification2'.
    ;
    ;3.out(variable): On exit, contains the total number of hexagons in the polygon construction for current frame.


    HexagonCountandPlot,defectClassification,unt,out
    numHexagons = out

    Defect_numbers = [frameNumber,numberOfParticles, numTriangles, numQuadrilaterals, numPentagons,numHexagons]

    ;The following 'PRINTF' procedure is peculiar to IDL.
    ;PRINT procedures perform formatted output
    ;It is used below to write row matrix 'Defect_numbers' to the file unit 3.
    printf,3,Defect_numbers

    oplot, X,Y, psym = 3, symsize = 1

    ;plotting the shock front position:
    plots, shockfront_x, shockFront_y, THICK = 3
    plots, postshock_right_x, postshock_right_y, COLOR = 47, THICK = 2
    plots, preshock_left_x, preshock_left_y, COLOR = 254, THICK = 2

    ;The following 'TVRD', 'TVLCT' and 'WRITETIFF' routines are peculiar to IDL.
    ;They are used below to save tiff images of polygon construction.
    scaled = TVRD()

    TVLCT, r, g, b, /Get
    s = Size(scaled, /Dimensions)
    image24 = BytArr(3, s[0], s[1])
    image24[0,*,*] = r[scaled]
    image24[1,*,*] = g[scaled]
    image24[2,*,*] = b[scaled]
    image24 = TVRD(True=1)
    image24 = Reverse(image24, 3)

    fname = datapath+ '\outputs\images\' + STRCOMPRESS(coreName+string(frameNumber,FORMAT='(I04)')+'.tif', /REMOVE_ALL)

    Write_Tiff, fname, image24, 1

    ;User defined procedure 10: 'VertexClassification'
    ;
    ;It is used below to perform a vertex classification.
    ;
    ;The arguments for this procedure are as follows.
    ;
    ;defectClassification(matrix): Updated matrix from procedure 'HexagonCountandPlot'.
    ;
    ;numberOfParticles(variable): Total number of particle positions in the current frame.
    ;
    ;Ulater(matrix): Updated matrix from procedure 'polygonIdentification'.
    ;
    ;frameNumber(variable): Current frame number.
    ;
    ;out(matrix): On exit, contains the number of vertices for various vertex types.

    VertexClassification,defectClassification,numberOfParticles,Ulater,frameNumber,out
    vertexStatnew = out
    printf,4,vertexStatnew,FORMAT='(F7.1,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1,1X,F7.1)'

  endfor
  close,3
  close,4
end

