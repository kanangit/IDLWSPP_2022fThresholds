pro driver_an_2022fThresholds_plgn

  datapath = 'e:\OneDrive - University of Iowa\bDocs\expAnalysisBackup\c_14226_vid57\20220606forP_2022fThresholds_plgn\'
  
  ;start and end frames
  iBegin = 1
  iEnd =  660
  
  ;start and end frames for pulse postition fitting:
  iBegin_ppulse = 522
  iEnd_ppulse = 624

  curDate='20220606'
  print, curDate
  print, datapath
  coreName = STRCOMPRESS('pgnConstr_' + STRING(curDate) + 'ff_', /REMOVE_ALL)
  ;ROI
  leftBorder = 600.0d
  rightBorder= 1000.0d;
  yMin = 0.0d;
  yMax = 1000.0d



  forceXlen = 1100
  forceYlen = 1200
  screenWidth = 1200

  plot_ybegin = 600;
  plot_yend = 1000;
  plot_xbegin = 0
  plot_xend = 1100
  screenWidth = forceXlen
  ratio = DOUBLE(plot_yend - plot_ybegin)/DOUBLE(plot_xend - plot_xbegin)




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

    ;The following 'PLOT' routine is perculiar to IDL.
    ;The PLOT procedure draws graphs of vector arguments.
    ;Below we use 'PLOT' for plotting particle positions.

    plot,X,Y,psym=3, isotropic=1, xrange = [plot_xbegin,plot_xend], yrange = [plot_ybegin,plot_yend], xstyle = 1, ystyle=1,title = 'polygon construction for frame' + string(frameNumber)
    ;User defined procedure 1: 'sideSelector'
    ;
    ;This procedure is used below to perform the following tasks.
    ;
    ;I. Identify bonds between particle positions by performing a Delaynay triangulation.
    ;II.Compute the angle value opposite to a given bond.
    ;
    ;The input arguments for this procedure are as follows.
    ;
    ; 1.X(row matrix): X coordinates of the particle positions in the current frame.
    ;
    ; 2.Y(row matrix): Y coordinates of the particle positions in the current frame.
    ;
    ; 3.allAngles(column matrix): On exit, contains the angle value opposite to a given bond.
    ;   This matrix contains such angle values for all the bonds.
    ;
    ; 4.allSides(matrix with four columns): On exit, contains the X and Y coordinates of the two end points of a bond.
    ;   This matrix contains such X and Y coordinates for all triangle sides.
    ;
    ;   for example consider row 1 of matrices 'allAngles' and 'allSides'
    ;   allSides (row1)
    ;   54.02   558.2   30.5   547.7
    ;   allAngles(row1)
    ;   65.2
    ;
    ;   This means that two end point coordinates of the bond in row 1 are (54.02,558.2) and (30.5,547.7)
    ;   and the angle opposite to this bond is 65.2 degrees.
    ;
    ;   Both the above matrices include duplicate entries corresponding to a particular bond as it is shared by maximum of
    ;   two triangles identified by the Delaunay triangulation routine.
    ;
    ; 5.numTriangles(variable):  On exit, contains the total number of triangles identified in current frame.

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

;*******User defined procedures ***********
;
;User defined procedure 1

pro sideSelector,X,Y,allAngles,allSides,numTriangles

  ;The following 'TRIANGULATE' routine is peculiar to IDL
  ;The TRIANGULATE procedure constructs a Delaunay triangulation of a planar set of points.
  ;Arguments for this routine are as follows,
  ;  1.X: An array that contains the X coordinates of the particles to be triangulated.
  ;  2.Y: An array that contains the Y coordinates of the points to be triangulated .
  ;  3.TRI: A named variable that, on exit, contains the list of triangles in the Delaunay
  ;    triangulation of the points specified by the X and Y arguments.

  triangulate, X,Y,TRI

  triangleSide = FLTARR(4,3)
  numberOfTriangles = n_elements(TRI[0,*])
  allAngles = FLTARR(1,numberOfTriangles*3)
  allSides = FLTARR(4,numberOfTriangles*3)

  for triangle=0L,numberOfTriangles-1 do begin

    triangleSide(0,0) = X[TRI [0, triangle]]
    triangleSide(1,0) = Y[TRI [0, triangle]]
    triangleSide(2,0) = X[TRI [1, triangle]]
    triangleSide(3,0) = Y[TRI [1, triangle]]

    triangleSide(0,1) = X[TRI [1, triangle]]
    triangleSide(1,1) = Y[TRI [1, triangle]]
    triangleSide(2,1) = X[TRI [2, triangle]]
    triangleSide(3,1) = Y[TRI [2, triangle]]

    triangleSide(0,2) = X[TRI [0, triangle]]
    triangleSide(1,2) = Y[TRI [0, triangle]]
    triangleSide(2,2) = X[TRI [2, triangle]]
    triangleSide(3,2) = Y[TRI [2, triangle]]

    ;The following 'SQRT' function is peculiar to IDL.
    ;The SQRT function computes the square root of a given value.
    ;It is used below to compute the length of each bond.

    side0= sqrt((triangleSide(0,0)-triangleSide(2,0))^2 + (triangleSide(1,0)-triangleSide(3,0))^2)
    side1= sqrt((triangleSide(0,1)-triangleSide(2,1))^2 + (triangleSide(1,1)-triangleSide(3,1))^2)
    side2= sqrt((triangleSide(0,2)-triangleSide(2,2))^2 + (triangleSide(1,2)-triangleSide(3,2))^2)

    ;The following 'ACOS' function is peculiar to IDL.
    ;The ACOS function returns an inverse cosine value, expressed in radians.
    ;It is used below to compute the magnitude of the angle opposite to a bond.

    angle0 = (acos((side1^2+side2^2-side0^2)/(2*side1*side2)))*(180/!pi) ;Angle opposite to side0
    angle1 = (acos((side0^2+side2^2-side1^2)/(2*side0*side2)))*(180/!pi) ;Angle opposite to side1
    angle2 = (acos((side0^2+side1^2-side2^2)/(2*side0*side1)))*(180/!pi) ;Angle opposite to side2

    row = triangle*3

    allAngles(0,row)   = angle0
    allAngles(0,row+1) = angle1
    allAngles(0,row+2) = angle2

    allSides(*,row)   = triangleSide(*,0)
    allSides(*,row+1) = triangleSide(*,1)
    allSides(*,row+2) = triangleSide(*,2)

  endfor

  numTriangles = numberOfTriangles

end

;User defined procedure 2

pro maxAngleSelect,allAngles,allSides,numTriangles,track

  increment = numTriangles*3

  track = fltarr(1,increment)

  for side=0L,increment-1 do begin

    if track(0,side) NE -5 then begin

      tempC1 = allSides(*,side)
      tempAngle1 = allAngles(0,side)

      for side2=side+1,increment-1 do begin

        if track(0,side2) NE -5 then begin

          tempC2 = allSides(*,side2)
          tempAngle2 = allAngles(0,side2)

          if tempC1(0,0) EQ tempC2(0,0)&& tempC1(1,0) EQ tempC2(1,0) && tempC1(2,0) EQ tempC2(2,0) && tempC1(3,0) EQ tempC2(3,0) then begin

            track(0,side2) = -5

            if tempAngle1 GE tempAngle2 then begin

              allAngles(0,side2) = tempAngle1

            endif

            if tempAngle2 GT tempAngle1 then begin

              allAngles(0,side) = tempAngle2

            endif

            break

          endif

          if tempC1(0,0) EQ tempC2(2,0)&& tempC1(1,0) EQ tempC2(3,0) && tempC1(2,0) EQ tempC2(0,0) && tempC1(3,0) EQ tempC2(1,0) then begin

            track(0,side2) = -5

            if tempAngle1 GE tempAngle2 then begin

              allAngles(0,side2) = tempAngle1

            endif

            if tempAngle2 GT tempAngle1 then begin

              allAngles(0,side) = tempAngle2

            endif

            break

          endif

        endif

      endfor

    endif

  endfor

end

;User defined procedure 3

pro selectAndPlotSides,track,allSides,allAngles,angleThreshold,sidesToPlot


  row = where(track NE -5)
  sides = allSides(*,row)
  angles = allAngles(0,row)

  row = where(angles(0,*) LE angleThreshold)
  sidesToPlot = sides(*,row)

  numSides = size(sidesToPlot,/DIMENSIONS)


  for k=0L,numSides(1,0)-1 do begin

    ;The following 'OPLOT' procedure is peculiar to IDL.
    ;OPLOT plots vector data over a previously-drawn plot.
    ;It is used below to plot bonds that are opposite angles less than or equal to the 'angleThreshold'.

    oplot, [sidesToPlot(0, k),sidesToPlot(2, k)],[sidesToPlot(1, k),sidesToPlot(3,k)],color = 0

  endfor

end

;User defined procedure 4

pro polygonIdentification,sidesToPlot,X,Y,numberOfParticles,Ulater,polygonTracker,defectClassification,unt

  sideArrangement = FLTARR(17,numberOfParticles)
  numofEdges = size(sidesToPlot,/DIMENSIONS)
  numofEdges = numofEdges(1,0)
  temp = 0
  xV= 0
  yV= 0
  tempArray = FLTARR(1,numberOfParticles)


  for particle=0L,numberOfParticles-1 do begin
    temp = 0
    xV = X(particle,0)
    yV = Y(particle,0)

    sideArrangement(temp,particle) = xV
    temp=temp+1
    sideArrangement(temp,particle) = yV

    found = 0

    for sides=0L, numofEdges-1 do begin

      x0 =sidesToPlot(0,sides)
      y0 =sidesToPlot(1,sides)
      x1 =sidesToPlot(2,sides)
      y1 =sidesToPlot(3,sides)

      if xV EQ x0 && yV EQ y0 then begin

        if temp+2 LT 16 then begin
          temp=temp+1
          sideArrangement(temp,particle)=x1
          temp=temp+1
          sideArrangement(temp,particle)=y1
        endif

      endif

      if xV EQ x1 && yV EQ y1 then begin

        if temp+2 LT 16 then begin
          temp=temp+1
          sideArrangement(temp,particle)=x0
          temp=temp+1
          sideArrangement(temp,particle)=y0
        endif

      endif

    endfor

    sideArrangement(16,particle)=temp
    tempArray(0,particle) = (temp-1)/2

  endfor

  Ulater = tempArray

  tempAngles = FLTARR(8,numberOfParticles)
  tempnum = 0
  tempvar = 0

  for particle=0L,numberOfParticles-1 do begin

    xV = sideArrangement(0,particle)
    yV = sideArrangement(1,particle)
    tempnum = sideArrangement(16,particle)
    tempvar = 0

    for vert=2, tempnum-1,2 do begin

      x1 = sideArrangement(vert,particle)
      y1 = sideArrangement(vert+1,particle)

      ;The following 'ATAN' function is peculiar to IDL.
      ;The ATAN function returns an inverse tangent value, expressed in radians.
      ;It is used below to compute the magnitude of the angle that a certain bond makes with the positive X axis.

      an = atan((y1-yV),(x1-xV))*180/!pi

      if an LT 0 then begin
        an = an+360
      endif

      tempAngles(tempvar,particle) = an
      tempvar = tempvar+1

    endfor
    tempAngles(7,particle)= tempvar

  endfor

  clockwise = INTARR(7,numberOfParticles)

  for particle=0L,numberOfParticles-1 do begin

    for column=0,6 do begin
      clockwise(column,particle) = -5
    endfor

  endfor

  for particle=0L,numberOfParticles-1 do begin

    tempvar = 1
    tempvar2 = tempAngles(7,particle)

    for ang=0L,tempvar2-1 do begin
      clockwise(ang,particle) = tempvar
      tempvar= tempvar+1
    endfor

  endfor

  tempvar  = 0
  tempvar2 = 0

  for particle=0L,numberOfParticles-1 do begin

    maxValue = tempAngles(7,particle)

    for column=0L,maxValue-1 do begin

      tempvar = tempAngles(column,particle)
      tempvarClock = clockwise(column,particle)

      for column2=column+1,maxValue-1 do begin

        tempvar2 = tempAngles(column2,particle)
        tempvar2Clock = clockwise(column2,particle)

        if tempvar2 LT tempvar then begin

          tempAngles(column,particle) = tempvar2
          clockwise(column,particle) = tempvar2Clock
          tempAngles(column2,particle) = tempvar
          clockwise(column2,particle) = tempvarClock
          tempvar = tempvar2
          tempvarClock = tempvar2Clock

        endif

      endfor

    endfor

  endfor


  sideArrangementClockwise = FLTARR(17,numberOfParticles)

  for particle=0L,numberOfParticles-1 do begin

    tempvar = tempAngles(7,particle)
    sideArrangementClockwise(16,particle) = tempAngles(7,particle)
    tempvar3 = 2
    sideArrangementClockwise(0,particle) = sideArrangement(0,particle)
    sideArrangementClockwise(1,particle) = sideArrangement(1,particle)

    for column=0L,tempvar-1 do begin
      tempvar2 = clockwise(column,particle)*2
      sideArrangementClockwise(tempvar3,particle) = sideArrangement(tempvar2,particle)
      tempvar3 = tempvar3+1
      sideArrangementClockwise(tempvar3,particle) = sideArrangement(tempvar2+1,particle)
      tempvar3 = tempvar3+1
    endfor

  endfor

  startVertexX = 0
  startVertexY = 0
  numBondsAttached = 0
  defectClassification = FLTARR(16,numberOfParticles*6)
  tempvar=0

  for particle=0L,numberOfParticles-1 do begin

    numBondsAttached = sideArrangementClockwise(16,particle)
    startVertexX = sideArrangementClockwise(0,particle)
    startVertexY = sideArrangementClockwise(1,particle)

    for column=1L,numBondsAttached do begin

      x1=sideArrangementClockwise(column*2,particle)
      y1=sideArrangementClockwise(column*2+1,particle)

      defectClassification(0,tempvar) = startVertexX
      defectClassification(1,tempvar) = startVertexY
      defectClassification(2,tempvar) = x1
      defectClassification(3,tempvar) = y1


      if column NE numBondsAttached then begin
        defectClassification(15,tempvar) = -10
      endif else begin
        defectClassification(15,tempvar) = numBondsAttached-1
      endelse
      tempvar = tempvar+1
    endfor

  endfor


  defectClassificationNew = FLTARR(16,tempvar)

  for defect=0L, tempvar-1 do begin
    defectClassificationNew(*,defect) = defectClassification(*,defect)
  endfor

  defectClassification = defectClassificationNew


  row  = 0
  unt = tempvar-1

  track = 0

  polygonTracker = fltarr(4,unt+1)
  polygonTracker(0,*) = defectClassification(0,*)
  polygonTracker(1,*) = defectClassification(1,*)
  polygonTracker(2,*) = defectClassification(2,*)
  polygonTracker(3,*) = defectClassification(3,*)

  for defect=0L, unt do begin

    x0 = defectClassification(0,defect)
    y0 = defectClassification(1,defect)
    x1 = defectClassification(2,defect)
    y1 = defectClassification(3,defect)
    row = defect
    row2 = 0
    found = 0

    trackrow = where(x1 EQ polygonTracker(0,*) and y1 EQ polygonTracker(1,*) and x0 EQ polygonTracker(2,*) and y0 EQ polygonTracker(3,*))
    if trackrow NE -1 && trackrow NE defect then begin
      x0temp = x1
      y0temp = y1
      track = defectClassification(15,trackrow)
      found = 1
      row2 = trackrow


      if (row2+1) LE unt then begin

        if track EQ -10 then begin

          row2= row2+1

        endif else begin

          row2 = row2-track

        endelse

        if found EQ 1 && defectClassification(0,(row2)) EQ x0temp && defectClassification(1,(row2)) EQ y0temp then begin

          defectClassification(4,defect) = defectClassification(2,(row2))
          defectClassification(5,defect) = defectClassification(3,(row2))
          defectClassification(14,defect) = 2

        endif else begin

          defectClassification(14,defect) = 1

        endelse

      endif

    endif

  endfor

end

;User defined procedure 5

pro polygonIdentification2,polygonTracker,defectClassification,numsides,start,unt

  tempvar = 0
  polygonTracker(0,*) = defectClassification(0,*)
  polygonTracker(1,*) = defectClassification(1,*)
  polygonTracker(2,*) = defectClassification(2,*)
  polygonTracker(3,*) = defectClassification(3,*)

  for defect=0, unt do begin
    tempvar=defectClassification(14,defect)

    if tempvar EQ numsides then begin

      x0 = defectClassification(start,defect)
      y0 = defectClassification(start+1,defect)
      x1 = defectClassification(start+2,defect)
      y1 = defectClassification(start+3,defect)
      row = defect
      row2 = 0
      found = 0
      trackrow = where(x1 EQ polygonTracker(0,*) and y1 EQ polygonTracker(1,*) and x0 EQ polygonTracker(2,*) and y0 EQ polygonTracker(3,*))

      if trackrow NE -1 && trackrow NE defect  then begin

        x0temp = x1
        y0temp = y1
        track = defectClassification(15,trackrow )
        found = 1
        row2 = trackrow

        if (row2+1) LE unt then begin

          if track EQ -10 then begin

            row2= row2+1
          endif else begin
            row2 = row2-track
          endelse

          if found EQ 1 && defectClassification(0,(row2)) EQ x0temp && defectClassification(1,(row2)) EQ y0temp then begin

            defectClassification(start+4,defect) = defectClassification(2,(row2))
            defectClassification(start+5,defect) = defectClassification(3,(row2))

            if defectClassification(start+4,defect) EQ defectClassification(0,defect) && defectClassification(start+5,defect) EQ defectClassification(1,defect) then begin

              defectClassification(14,defect) = 10*(numsides+1)+numsides+1
            endif else begin
              defectClassification(14,defect) = numsides+1
            endelse

          endif else begin
            defectClassification(14,defect) = numsides
          endelse
        endif
      endif
    endif
  endfor
end

;User defined procedure 6

pro TriangleCountandPlot,defectClassification,unt,out

  polygonCount = fltarr(6,unt+1)
  polygonCountTrack = fltarr(1,unt+1)
  polygonCount(0,*) = defectClassification(0,*)
  polygonCount(1,*) = defectClassification(1,*)
  polygonCount(2,*) = defectClassification(2,*)
  polygonCount(3,*) = defectClassification(3,*)
  polygonCount(4,*) = defectClassification(4,*)
  polygonCount(5,*) = defectClassification(5,*)

  polygonCountTrack(0,*) = defectClassification(14,*)

  row= where(polygonCountTrack(0,*) EQ 33)
  numPolygon = size(row,/N_ELEMENTS)

  if numPolygon GT 1 then begin
    polygonCount = polygonCount(*,row)
    numPolygon = size(polygonCount)

    numPolygon =  numPolygon(2,0)
    coordinates = fltarr(6,unt)

    for num=0, numPolygon-1 do begin

      Xcs = [polygonCount(0,num),polygonCount(2,num),polygonCount(4,num)]
      Ycs = [polygonCount(1,num),polygonCount(3,num),polygonCount(5,num)]

      ;The following 'SORT' function is peculiar to IDL.
      ;SORT arranges the element in the argument in the ascending order.
      ;It is used below to arrange the X coordinates of the vertices of the polygon in the ascending order.
      ;This is useful in identifying duplicate entries of polygons in the matrix 'defectClassification'.

      Xcs = Xcs(sort(Xcs))
      Ycs = Ycs(sort(Ycs))

      coordinates(0,num) = Xcs(0,0)
      coordinates(1,num) = Ycs(0,0)
      coordinates(2,num) = Xcs(1,0)
      coordinates(3,num) = Ycs(1,0)
      coordinates(4,num) = Xcs(2,0)
      coordinates(5,num) = Ycs(2,0)

    endfor

    track = fltarr(1,numPolygon)

    for num=0,numPolygon-1 do begin
      if track(0,num) NE -5 then begin

        tempCoord = coordinates(*,num)

        row = where(coordinates(0,*) EQ tempCoord(0,0) and coordinates(1,*) EQ tempCoord(1,0) and coordinates(2,*) EQ tempCoord(2,0) and coordinates(3,*) EQ tempCoord(3,0) and coordinates(4,*) EQ tempCoord(4,0) and coordinates(5,*) EQ tempCoord(5,0))

        tempnum = size(row)

        for num2=0,tempnum(1,0)-1 do begin

          track (0,row(num2,0)) = -5

        endfor
        track(0,num) = 5

      endif
    endfor


    row = where(track EQ 5)
    Polygonswithoutrepeats = polygonCount(*,row)
    numTriangles= size(Polygonswithoutrepeats,/N_ELEMENTS)
    numTriangles = numTriangles/6



    for num=0,numTriangles-1 do begin

      Xcs = [Polygonswithoutrepeats(0,num),Polygonswithoutrepeats(2,num),Polygonswithoutrepeats(4,num),Polygonswithoutrepeats(0,num)]
      Ycs = [Polygonswithoutrepeats(1,num),Polygonswithoutrepeats(3,num),Polygonswithoutrepeats(5,num),Polygonswithoutrepeats(1,num)]

      ;The following procedure 'POLYFILL' is peculiar to IDL.
      ;The POLYFILL procedure fills the interior of a region of the display enclosed by an arbitrary two-dimensional polygon.
      ;It is used below to plot and color the identified polygon.

      polyfill, Xcs,Ycs, color=140,linestyle = 0,thick = 0.5 ;color = 140 - green
      oplot, Xcs,Ycs

    endfor
  endif else begin
    numTriangles = 0
  endelse
  out = numTriangles
end

;User defined procedure 7

pro QuadrilateralCountandPlot,defectClassification,unt,out


  polygonCount = fltarr(8,unt+1)
  polygonCountTrack = fltarr(1,unt+1)
  polygonCount(0,*) = defectClassification(0,*)
  polygonCount(1,*) = defectClassification(1,*)
  polygonCount(2,*) = defectClassification(2,*)
  polygonCount(3,*) = defectClassification(3,*)
  polygonCount(4,*) = defectClassification(4,*)
  polygonCount(5,*) = defectClassification(5,*)
  polygonCount(6,*) = defectClassification(6,*)
  polygonCount(7,*) = defectClassification(7,*)

  polygonCountTrack(0,*) = defectClassification(14,*)

  row= where(polygonCountTrack(0,*) EQ 44)
  numPolygon = size(row,/N_ELEMENTS)

  if numPolygon GT 1 then begin
    polygonCount = polygonCount(*,row)
    numPolygon = size(polygonCount)

    numPolygon =  numPolygon(2,0)
    coordinates = fltarr(8,unt)

    for num=0, numPolygon-1 do begin

      Xcs = [polygonCount(0,num),polygonCount(2,num),polygonCount(4,num),polygonCount(6,num)]
      Ycs = [polygonCount(1,num),polygonCount(3,num),polygonCount(5,num),polygonCount(7,num)]

      Xcs = Xcs(sort(Xcs))
      Ycs = Ycs(sort(Ycs))

      coordinates(0,num) = Xcs(0,0)
      coordinates(1,num) = Ycs(0,0)
      coordinates(2,num) = Xcs(1,0)
      coordinates(3,num) = Ycs(1,0)
      coordinates(4,num) = Xcs(2,0)
      coordinates(5,num) = Ycs(2,0)
      coordinates(6,num) = Xcs(3,0)
      coordinates(7,num) = Ycs(3,0)

    endfor

    track = fltarr(1,numPolygon)

    for num=0,numPolygon-1 do begin
      if track(0,num) NE -5 then begin

        tempCoord = coordinates(*,num)

        row = where(coordinates(0,*) EQ tempCoord(0,0) and coordinates(1,*) EQ tempCoord(1,0) and coordinates(2,*) EQ tempCoord(2,0) and coordinates(3,*) EQ tempCoord(3,0) and coordinates(4,*) EQ tempCoord(4,0) and coordinates(5,*) EQ tempCoord(5,0)and coordinates(6,*) EQ tempCoord(6,0) and coordinates(7,*) EQ tempCoord(7,0))

        tempnum = size(row)

        for num2=0,tempnum(1,0)-1 do begin

          track (0,row(num2,0)) = -5

        endfor
        track(0,num) = 5

      endif
    endfor


    row = where(track EQ 5)
    Polygonswithoutrepeats = polygonCount(*,row)
    numQuadrilaterals= size(Polygonswithoutrepeats,/N_ELEMENTS)
    numQuadrilaterals = numQuadrilaterals/8



    for num=0,numQuadrilaterals-1 do begin

      Xcs = [Polygonswithoutrepeats(0,num),Polygonswithoutrepeats(2,num),Polygonswithoutrepeats(4,num),Polygonswithoutrepeats(6,num),Polygonswithoutrepeats(0,num)]
      Ycs = [Polygonswithoutrepeats(1,num),Polygonswithoutrepeats(3,num),Polygonswithoutrepeats(5,num),Polygonswithoutrepeats(7,num),Polygonswithoutrepeats(1,num)]
      polyfill, Xcs,Ycs, color=250,linestyle = 0,thick = 0.5 ;color = 250 - red
      oplot, Xcs,Ycs

    endfor
  endif else begin
    numQuadrilaterals = 0
  endelse
  out = numQuadrilaterals
end

;User defined procedure 8

pro PentagonsCountandPlot,defectClassification,unt,out


  polygonCount = fltarr(10,unt+1)
  polygonCountTrack = fltarr(1,unt+1)
  polygonCount(0,*) = defectClassification(0,*)
  polygonCount(1,*) = defectClassification(1,*)
  polygonCount(2,*) = defectClassification(2,*)
  polygonCount(3,*) = defectClassification(3,*)
  polygonCount(4,*) = defectClassification(4,*)
  polygonCount(5,*) = defectClassification(5,*)
  polygonCount(6,*) = defectClassification(6,*)
  polygonCount(7,*) = defectClassification(7,*)
  polygonCount(8,*) = defectClassification(8,*)
  polygonCount(9,*) = defectClassification(9,*)

  polygonCountTrack(0,*) = defectClassification(14,*)

  row= where(polygonCountTrack(0,*) EQ 55)
  numPolygon = size(row,/N_ELEMENTS)

  if numPolygon GT 1 then begin
    polygonCount = polygonCount(*,row)
    numPolygon = size(polygonCount)

    numPolygon =  numPolygon(2,0)
    coordinates = fltarr(10,unt)

    for num=0, numPolygon-1 do begin

      Xcs = [polygonCount(0,num),polygonCount(2,num),polygonCount(4,num),polygonCount(6,num),polygonCount(8,num)]
      Ycs = [polygonCount(1,num),polygonCount(3,num),polygonCount(5,num),polygonCount(7,num),polygonCount(9,num)]

      Xcs = Xcs(sort(Xcs))
      Ycs = Ycs(sort(Ycs))

      coordinates(0,num) = Xcs(0,0)
      coordinates(1,num) = Ycs(0,0)
      coordinates(2,num) = Xcs(1,0)
      coordinates(3,num) = Ycs(1,0)
      coordinates(4,num) = Xcs(2,0)
      coordinates(5,num) = Ycs(2,0)
      coordinates(6,num) = Xcs(3,0)
      coordinates(7,num) = Ycs(3,0)
      coordinates(8,num) = Xcs(4,0)
      coordinates(9,num) = Ycs(4,0)

    endfor

    track = fltarr(1,numPolygon)

    for num=0,numPolygon-1 do begin
      if track(0,num) NE -5 then begin

        tempCoord = coordinates(*,num)

        row = where(coordinates(0,*) EQ tempCoord(0,0) and coordinates(1,*) EQ tempCoord(1,0) and coordinates(2,*) EQ tempCoord(2,0) and coordinates(3,*) EQ tempCoord(3,0) and coordinates(4,*) EQ tempCoord(4,0) and coordinates(5,*) EQ tempCoord(5,0)and coordinates(6,*) EQ tempCoord(6,0) and coordinates(7,*) EQ tempCoord(7,0)and coordinates(8,*) EQ tempCoord(8,0)and coordinates(9,*) EQ tempCoord(9,0))

        tempnum = size(row)

        for num2=0,tempnum(1,0)-1 do begin

          track (0,row(num2,0)) = -5

        endfor
        track(0,num) = 5

      endif
    endfor


    row = where(track EQ 5)
    Polygonswithoutrepeats = polygonCount(*,row)
    numPentagons= size(Polygonswithoutrepeats,/N_ELEMENTS)
    numPentagons = numPentagons/10



    for num=0, numPentagons-1 do begin

      Xcs = [Polygonswithoutrepeats(0,num),Polygonswithoutrepeats(2,num),Polygonswithoutrepeats(4,num),Polygonswithoutrepeats(6,num),Polygonswithoutrepeats(8,num),Polygonswithoutrepeats(0,num)]
      Ycs = [Polygonswithoutrepeats(1,num),Polygonswithoutrepeats(3,num),Polygonswithoutrepeats(5,num),Polygonswithoutrepeats(7,num),Polygonswithoutrepeats(9,num),Polygonswithoutrepeats(1,num)]
      polyfill, Xcs,Ycs, color=200,linestyle = 0,thick = 0.5 ;color = 200 - yellow
      oplot, Xcs,Ycs

    endfor
  endif else begin
    numPentagons = 0
  endelse
  out = numPentagons
end

;User defined procedure 9

pro HexagonCountandPlot,defectClassification,unt,out


  polygonCount = fltarr(12,unt+1)
  polygonCountTrack = fltarr(1,unt+1)
  polygonCount(0,*) = defectClassification(0,*)
  polygonCount(1,*) = defectClassification(1,*)
  polygonCount(2,*) = defectClassification(2,*)
  polygonCount(3,*) = defectClassification(3,*)
  polygonCount(4,*) = defectClassification(4,*)
  polygonCount(5,*) = defectClassification(5,*)
  polygonCount(6,*) = defectClassification(6,*)
  polygonCount(7,*) = defectClassification(7,*)
  polygonCount(8,*) = defectClassification(8,*)
  polygonCount(9,*) = defectClassification(9,*)
  polygonCount(10,*) = defectClassification(10,*)
  polygonCount(11,*) = defectClassification(11,*)

  polygonCountTrack(0,*) = defectClassification(14,*)

  row= where(polygonCountTrack(0,*) EQ 66)
  numPolygon = size(row,/N_ELEMENTS)

  if numPolygon GT 1 then begin
    polygonCount = polygonCount(*,row)
    numPolygon = size(polygonCount)

    numPolygon =  numPolygon(2,0)
    coordinates = fltarr(12,unt)

    for num=0, numPolygon-1 do begin

      Xcs = [polygonCount(0,num),polygonCount(2,num),polygonCount(4,num),polygonCount(6,num),polygonCount(8,num),polygonCount(10,num)]
      Ycs = [polygonCount(1,num),polygonCount(3,num),polygonCount(5,num),polygonCount(7,num),polygonCount(9,num),polygonCount(11,num)]

      Xcs = Xcs(sort(Xcs))
      Ycs = Ycs(sort(Ycs))

      coordinates(0,num) = Xcs(0,0)
      coordinates(1,num) = Ycs(0,0)
      coordinates(2,num) = Xcs(1,0)
      coordinates(3,num) = Ycs(1,0)
      coordinates(4,num) = Xcs(2,0)
      coordinates(5,num) = Ycs(2,0)
      coordinates(6,num) = Xcs(3,0)
      coordinates(7,num) = Ycs(3,0)
      coordinates(8,num) = Xcs(4,0)
      coordinates(9,num) = Ycs(4,0)
      coordinates(10,num) = Xcs(5,0)
      coordinates(11,num) = Ycs(5,0)

    endfor

    track = fltarr(1,numPolygon)

    for num=0,numPolygon-1 do begin
      if track(0,num) NE -5 then begin

        tempCoord = coordinates(*,num)

        row = where(coordinates(0,*) EQ tempCoord(0,0) and coordinates(1,*) EQ tempCoord(1,0) and coordinates(2,*) EQ tempCoord(2,0) and coordinates(3,*) EQ tempCoord(3,0) and coordinates(4,*) EQ tempCoord(4,0) and coordinates(5,*) EQ tempCoord(5,0)and coordinates(6,*) EQ tempCoord(6,0) and coordinates(7,*) EQ tempCoord(7,0)and coordinates(8,*) EQ tempCoord(8,0)and coordinates(9,*) EQ tempCoord(9,0)and coordinates(10,*) EQ tempCoord(10,0)and coordinates(11,*) EQ tempCoord(11,0))

        tempnum = size(row)

        for num2=0,tempnum(1,0)-1 do begin

          track (0,row(num2,0)) = -5

        endfor
        track(0,num) = 5

      endif
    endfor


    row = where(track EQ 5)
    Polygonswithoutrepeats = polygonCount(*,row)
    numHexagons= size(Polygonswithoutrepeats,/N_ELEMENTS)
    numHexagons = numHexagons/12



    for num=0, numHexagons-1 do begin

      Xcs = [Polygonswithoutrepeats(0,num),Polygonswithoutrepeats(2,num),Polygonswithoutrepeats(4,num),Polygonswithoutrepeats(6,num),Polygonswithoutrepeats(8,num),Polygonswithoutrepeats(10,num),Polygonswithoutrepeats(0,num)]
      Ycs = [Polygonswithoutrepeats(1,num),Polygonswithoutrepeats(3,num),Polygonswithoutrepeats(5,num),Polygonswithoutrepeats(7,num),Polygonswithoutrepeats(9,num),Polygonswithoutrepeats(11,num),Polygonswithoutrepeats(1,num)]
      polyfill, Xcs,Ycs, color=215,linestyle = 0,thick = 0.5 ;color = 215 - orange
      oplot, Xcs,Ycs

    endfor
  endif else begin
    numHexagons = 0
  endelse
  out = numHexagons





end

;User defined procedure 10

pro VertexClassification,defectClassification,numberOfParticles,Ulater,frameNumber,out

  tempdefectClass = [defectClassification(0,*),defectClassification(1,*),defectClassification(14,*),defectClassification(15,*)]

  vertexClass = FLTARR(11,numberOfParticles)
  tempVar =0
  stop = 0
  start = 0

  for particle=0,numberOfParticles-1 do begin

    stop = Ulater(0,particle)+stop
    vertexClass(9,particle)=Ulater(0,particle)
    tempVar = 2

    for num =start,stop-1 do begin

      vertexClass(tempVar,particle) = tempdefectClass(2,num)/11
      tempVar = tempVar+1

    endfor
    vertexClass(0,particle)=tempdefectClass(0,start)
    vertexClass(1,particle)=tempdefectClass(1,start)
    start=stop

  endfor


  for particle =0, numberOfParticles-1 do begin
    tempVar = vertexClass(9,particle)
    threes =0
    fours = 0
    fives = 0
    sixes = 0
    found = 0

    for num=2,tempVar+1 do begin

      tempVar2 = vertexClass(num,particle)
      if tempVar2 EQ 3 then begin
        threes = threes+1
      endif
      if tempVar2 EQ 4 then begin
        fours = fours+1
      endif
      if tempVar2 EQ 5 then begin
        fives = fives+1
      endif
      if tempVar2 EQ 6 then begin
        sixes = sixes+1
      endif
      if tempVar2 NE 3 && tempVar2 NE 4 && tempVar2 NE 5 && tempVar2 NE 6 then begin
        found = 1
      endif
    endfor


    if found EQ 1 then begin
      vertexClass(10,particle)=27
    endif
    found =0

    if threes EQ 7 && fours EQ 0 && fives EQ 0 && sixes EQ 0 && tempVar EQ 7 then begin
      vertexClass(10,particle)=11
    endif

    if threes EQ 6 && fours EQ 0 && fives EQ 0 && sixes EQ 0 && tempVar EQ 6 then begin
      vertexClass(10,particle)=1
    endif

    if threes EQ 5 && fours EQ 1 && fives EQ 0 && sixes EQ 0 && tempVar EQ 6 then begin
      vertexClass(10,particle)=6
    endif

    if threes EQ 5 && fours EQ 0 && fives EQ 1 && sixes EQ 0 && tempVar EQ 6 then begin
      vertexClass(10,particle)=17
    endif


    if threes EQ 4 && fours EQ 1 && fives EQ 0 && sixes EQ 0 && tempVar EQ 5 then begin
      vertexClass(10,particle)=5
    endif

    if threes EQ 5 && fours EQ 0 && fives EQ 0 && sixes EQ 0 && tempVar EQ 5 then begin
      vertexClass(10,particle)=10
    endif

    if threes EQ 4 && fours EQ 0 && fives EQ 1 && sixes EQ 0 && tempVar EQ 5 then begin
      vertexClass(10,particle)=12
    endif

    if threes EQ 4 && fours EQ 0 && fives EQ 0 && sixes EQ 1 && tempVar EQ 5 then begin
      vertexClass(10,particle)=20
    endif

    if threes EQ 2 && fours EQ 2 && fives EQ 1 && sixes EQ 0 && tempVar EQ 5 then begin
      vertexClass(10,particle)=23
    endif

    if threes EQ 0 && fours EQ 4 && fives EQ 0 && sixes EQ 0 && tempVar EQ 4 then begin
      vertexClass(10,particle)=4
    endif

    if threes EQ 1 && fours EQ 3 && fives EQ 0 && sixes EQ 0 && tempVar EQ 4 then begin
      vertexClass(10,particle)=7
    endif

    if threes EQ 2 && fours EQ 0 && fives EQ 2 && sixes EQ 0 && tempVar EQ 4 then begin
      vertexClass(10,particle)=21
    endif

    if threes EQ 2 && fours EQ 1 && fives EQ 0 && sixes EQ 1 && tempVar EQ 4 then begin
      vertexClass(10,particle)=25
    endif

    if threes EQ 3 && fours EQ 2 && fives EQ 0 && sixes EQ 0 && tempVar EQ 5 then begin
      vertexClass(10,particle)=2
    endif

    if threes EQ 2 && fours EQ 3 && fives EQ 0 && sixes EQ 0 && tempVar EQ 5 then begin
      vertexClass(10,particle)=8
    endif

    if threes EQ 3 && fours EQ 1 && fives EQ 1 && sixes EQ 0 && tempVar EQ 5 then begin
      vertexClass(10,particle)=13
    endif

    if threes EQ 3 && fours EQ 0 && fives EQ 2 && sixes EQ 0 && tempVar EQ 5 then begin
      vertexClass(10,particle)=22
    endif

    if threes EQ 2 && fours EQ 1 && fives EQ 1 && sixes EQ 0 && tempVar EQ 4 then begin
      vertexClass(10,particle)=15
    endif

    if threes EQ 1 && fours EQ 2 && fives EQ 1 && sixes EQ 0 && tempVar EQ 4 then begin
      vertexClass(10,particle)=16
    endif


  endfor

  vert5 = [[4,4,3,3,3],[4,3,4,3,3],[4,4,4,3,3],[4,4,3,4,3],[5,4,3,3,3],[5,3,4,3,3],[5,5,3,3,3],[5,3,5,3,3],[5,4,4,3,3]]
  vert5track =[2,3,8,9,13,14,22,24,23]

  arrangement = 0


  for particle = 0, numberOfParticles-1 do begin
    track = 0
    tempvar = vertexClass(10,particle)
    found = 0

    if tempvar EQ 2 or tempvar EQ 8 or tempvar EQ 13 or tempvar EQ 22 or tempvar EQ 23 then begin
      arrangement = [vertexClass(2,particle),vertexClass(3,particle),vertexClass(4,particle),vertexClass(5,particle),vertexClass(6,particle)]


      for num=0,8 do begin
        tempVert = vert5(*,num)

        for m=1,5 do begin



          temp = tempVert(0,0)
          tempVert(0,0) = tempVert(1,0)
          tempVert(1,0) = tempVert(2,0)
          tempVert(2,0) = tempVert(3,0)
          tempVert(3,0) = tempVert(4,0)
          tempVert(4,0) = temp


          if (arrangement(0,0) EQ tempVert(0,0) && arrangement(1,0) EQ tempVert(1,0)&& arrangement(2,0) EQ tempVert(2,0)&& arrangement(3,0) EQ tempVert(3,0) && arrangement(4,0) EQ tempVert(4,0)) or (arrangement(0,0) EQ tempVert(4,0) && arrangement(1,0) EQ tempVert(3,0)&& arrangement(2,0) EQ tempVert(2,0)&& arrangement(3,0) EQ tempVert(1,0) && arrangement(4,0) EQ tempVert(0,0)) then begin
            track = num
            found = 1
            vertexClass(10,particle)= vert5track(track,0)
          endif

        endfor

      endfor

      if found NE 1 then begin
        vertexClass(10,particle) = 0
      endif

    endif

  endfor


  vert4 = [[5,4,3,3],[5,3,4,3],[5,4,4,3],[4,5,4,3],[6,4,3,3],[5,5,3,3]]
  vert4track =[15,18,16,19,25,21]

  arrangement = 0


  for particle = 0, numberOfParticles-1 do begin
    track = 0
    tempvar = vertexClass(10,particle)
    found = 0
    if tempvar EQ 15 or tempvar EQ 16 or tempvar EQ 25 or tempvar EQ 21  then begin
      arrangement = [vertexClass(2,particle),vertexClass(3,particle),vertexClass(4,particle),vertexClass(5,particle)]


      for num=0,5 do begin
        tempVert = vert4(*,num)

        for m=1,4 do begin

          temp = tempVert(0,0)
          tempVert(0,0) = tempVert(1,0)
          tempVert(1,0) = tempVert(2,0)
          tempVert(2,0) = tempVert(3,0)
          tempVert(3,0) = temp

          if (arrangement(0,0) EQ tempVert(0,0) && arrangement(1,0) EQ tempVert(1,0)&& arrangement(2,0) EQ tempVert(2,0)&& arrangement(3,0) EQ tempVert(3,0)) or (arrangement(0,0) EQ tempVert(3,0) && arrangement(1,0) EQ tempVert(2,0)&& arrangement(2,0) EQ tempVert(1,0)&& arrangement(3,0) EQ tempVert(0,0)) then begin
            track = num
            found = 1
            vertexClass(10,particle)= vert4track(track,0)
          endif

        endfor

      endfor

      if found NE 1 then begin
        vertexClass(10,particle) = 0
      endif


    endif

  endfor


  vertexStat = FLTARR(2,27)

  for num=0,25 do begin

    temp = num+1
    vertexStat(0,num) = temp
    tempvar=0

    for particle=0,numberOfParticles-1 do begin

      if vertexClass(10,particle) EQ temp then begin
        tempvar = tempvar+1
      endif

    endfor

    vertexStat(1,num)= tempvar
  endfor

  tempvar=0
  tempvar2=0

  for particle=0,numberOfParticles-1 do begin

    if vertexClass(10,particle) EQ 27 then begin

      tempvar2 = tempvar2+1

    endif

    if vertexClass(10,particle) EQ 0 then begin
      tempvar = tempvar+1
      vertexClass(10,particle) = 26
    endif

  endfor


  vertexStat(1,25)= tempvar
  vertexStat(1,26)= tempvar2
  vertexStat(0,26)= 27

  vertexStatnew = fltarr(27,1)
  vertexStatnew(0,0) = frameNumber
  b= total(vertexStat(1,*))
  vertexStatnew(1,0) = b - tempvar2

  for num=2,26 do begin
    vertexStatnew(num,0) = vertexStat(1,num-2)
  endfor

  out = vertexStatnew
end