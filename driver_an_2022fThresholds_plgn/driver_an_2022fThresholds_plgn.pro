pro driver_an_2022fThresholds_plgn


  ;Supplementary material to W. D. Suranga Ruhunusiri, J. Goree, Yan Feng, and Bin Liu,
  ;   "Polygon construction to investigate melting in 2D strongly-coupled dusty plasma",
  ;   to be published in Physical Review E, 2011
  ;
  ;This is the source code that we used to identify geometrical defects in dusty plasma,
  ;but it will also be useful for other physical systems and simulated systems.
  ;
  ;About the language:
  ;   This code is written in IDL, and has been tested with IDL version 7.1;
  ;   Hints are provided below for readers interested in porting this code to another language.
  ;   We have also provided a source code in MATLAB in the supplementay material as an example of a ported code.
  ;
  ;The polygon construction method is based on M. A. Glaser and N. A. Clark,
  ;   Phys. Rev. A 41, 4585 (1990), Adv. Chem. Phys. 83, 543 (1993).
  ;
  ;Organization of the code is as follows:
  ;1.Read a particle position data file (INPUT1) in text format. An example input data file is provided with this distribution.
  ;
  ;2.Identify bonds between particle positions by Delaunay triangulation. Here, we use the TRIANGULATION routine in IDL, but
  ;   the reader who wishes to port to another language can use another Delaunay triangulation routine.
  ;
  ;3.Bonds are removed that are opposite to bond angles greater than a user-adjustable threshold, (INPUT3).
  ;   The value of this threshold is hardwired in the code; as written here it is 75 degrees, but the user may change this.
  ;   This method of removing bonds is one of two methods used by Glaser and Clark (1990,1993); their other method (removing bonds longer
  ;   than a certain threshold) is not implemented here.
  ;   Here and elsewhere in the code we use an IDL routine WHERE which searches the contents of an array for a particular value
  ;   and returns the corresponding indices; users wishing to port the code will need to find a suitable substitute.
  ;
  ;4.Diagnostics peculiar to IDL
  ;   We generate one image, for each set of input particle positions. This image is a color presentation of the polygons.
  ;   As an example of this image, see Fig. 2(b) and (c) in our paper.
  ;   Porting this graphics portion of the code to another language should not be difficult.
  ;
  ;5.Diagnostics not peculiar to IDL
  ;   A text file is output, and it contains:
  ;      - counts of the different kinds of polygons (triangles, quadrilaterals, pentagons, etc.)
  ;      - counts of vertex types (Type A, Type B, etc., as described in our paper and by Glaser and Clark.)                                                         }
  ;
  ;                                                       }
  ;ALL the INPUTS, other than the input text file for particle positions, are HARDWIRED into this code. This code has no user dialog.                                                                         }
  ;
  ;Before executing the code, the user should edit (as required) all INPUTS 1 through 8, as listed below:                                      }
  ;
  ;
  ;INPUT1: data_file : Enter the path of a txt file containing three columns
  ;                     (First column contains the X coordinates of particles,
  ;                      second column contains the Y coordinates of particles,
  ;                      third column contains the frame number
  ;                      where the columns are delimited by an arbitrary number of spaces.
  ;                      This txt file has no header line)
  ;
  ;                      as an example of a line in this text file:
  ;
  ;                      491.868         20.5              0
  ;
  ;                      where the first two columns are the floating-point X and Y coordinates of a particle, and the third column is an integer for the "frame number".
  ;                      Here, the term "frame number" corresponds to a snapshot in time. Typically, one would have several thousand particles for each frame number.                                                                              }
  ;
  ;INPUT2: length_data_file(integer): the number of rows in the above txt file
  ;
  ;
  ;INPUT3: begin_frame (integer): the first frame number in the above txt file
  ;        For example, in the sample particle position data file given (particle_positions.txt) has particle position
  ;        data for fames 0 through 1
  ;        So you should specify 0 for this input
  ;
  ;INPUT4: last_frame (integer):the last frame number in the above txt file
  ;        For example, in the sample particle position data file given (particle_positions.txt) has particle position
  ;        data for fames 0 through 1
  ;        So you should specify 1 for this input
  ;
  ;INPUT5: angleThreshold(float):the threshold angle for bond removal
  ;        Default value is set at 75.0 degrees as suggested by M. A. Glaser and N. A. Clark
  ;
  ;INPUT6: image_file(string): folder location to save polygon construction images.
  ;        The format of the images are TIF.
  ;
  ;INPUT7: fname_defect_stat(string): folder location and a name for a txt file to save geometrical defect counts for
  ;        polygon construction images
  ;
  ;        as an example of a line in this text file:
  ;
  ;        0.000000      1212.00         2108.00         75.0000              7.00000         0.000000
  ;
  ;        First column- frame number, second column- number of particles, third column- number of triangles,
  ;        fourth column- number of quadrilaterals, fifth column- number of pentagons
  ;        sixth column- number of hexagons.
  ;
  ;INPUT8: fname_vertex_stat(string): folder location and a name for a txt file to save vertex type counts for
  ;        polygon construction images
  ;
  ;        as an example of a line in this text file:
  ;
  ;        0.0        851.0        37.0        21.0   ....   0.0
  ;
  ;        First column- frame number, second column- total number of vertices, third column- number of A vertices,
  ;        fourth column- number of B vertices, etc...
  ;
  ;
  ;After reading the manuscript and with some programming knowledge the reader will be able to
  ;understand the self-explanatory variables.
  ;

  forceXlen = 1600
  forceYlen = 900
  screenWidth = forceXlen

  ;INPUT1
  data_file = '\\BRECKENRIDGE4\expAnalysisBackup\c_14226_vid56\20181201temp\04_code_polygon_construction_IDL\inputs\ff1-660_20181119positionTrimmedForVid056_6400solid_15427284684POLYGON.txt'

  ;INPUT2
  length_data_file = 1010884

  ;INPUT3
  startFrame = 550

  ;INPUT4
  endFrame = 551

  ;INPUT5
  angleThreshold = 75.0

  ;INPUT6
  image_file = '\\BRECKENRIDGE4\expAnalysisBackup\c_14226_vid56\20181201temp\04_code_polygon_construction_IDL\outputs'

  ;INPUT7
  fname_defect_stat = '\\BRECKENRIDGE4\expAnalysisBackup\c_14226_vid56\20181201temp\04_code_polygon_construction_IDL\outputs\Defect_stat.txt'

  ;INPUT8
  fname_vertex_stat = '\\BRECKENRIDGE4\expAnalysisBackup\c_14226_vid56\20181201temp\04_code_polygon_construction_IDL\outputs\Vertex_stat.txt'


  ;The following 'CLOSE' procedure is peculiar to IDL.
  ;The CLOSE procedure closes the file units specified as arguments.
  ;It is used here to close the files 'data_file','fname_defect_stat' and 'fname_vertex_stat', if they are already open.
  close,1 ;File unit 1 will be assigned to 'data_file'
  close,3 ;File unit 3 will be assigned to 'fname_defect_stat'
  close,4 ;File unit 4 will be assigned to 'fname_vertex_stat'

  length = length_data_file

  ;The following 'FLTARR' function is peculiar to IDL.
  ;The FLTARR function creates a floating-point vector or array of the specified dimensions.
  ;It is used here to create a floating-point array with 3 columns and 'length' of rows.
  tempPositions = fltarr(3,length)

  ;The following 'OPENR' procedure is peculiar to IDL.
  ;OPENR opens an existing file for input only.
  ;It is used here to open 'data_file' for input only.
  openr,1,data_file

  ;The followig 'READF' procedure is peculiar to IDL.
  ;The READ procedure performs formatted input into variables.
  ;It is used here to read 'data_file' to the array named 'tempPositions'.
  readf,1,tempPositions
  close,1

  ;The following 'OPENW' procedure is peculiar to IDL.
  ;OPENW opens a new file for output. If the file exists, its old contents are destroyed.
  ;It is used here to open 'fname_defect_stat' for output.
  openw,3,fname_defect_stat ;File unit 3 is assigned to 'Defect_stat.txt'

  openw,4,fname_vertex_stat ;File unit 4 is assigned to 'Vertex_stat.txt'

  arrayParticlePositions = tempPositions

  for frameNumber=startFrame,endFrame do begin

    ;The following 'WHERE' function is peculiar to IDL.
    ;WHERE selects elements of an array satisfying a given criteria.
    ;It is used here to find the matrix row indices of the second column in matrix 'arrayParticlePositions'
    ;where the matrix element has a value equal to 'frameNumber'. This is useful in selecting the X and Y coordinates
    ;of particles in the current frame.
    row = where(arrayParticlePositions(2,*) EQ frameNumber)

    Xtemp= arrayParticlePositions(0,row) ;X coordinates of particle positions in the current frame.
    Ytemp= arrayParticlePositions(1,row) ;Y coordinates of particle positions in the current frame.

    ;The following 'SIZE' function is peculiar to IDL.
    ;SIZE function returns size and type information for its argument.
    ;It is used here to find the number of rows in the matrix 'Xtemp'.
    ;The number of rows in matrix 'Xtemp' is equal to the total number of particles in the current frame.
    numberOfParticles = size(Xtemp,/N_ELEMENTS) ;Total number of particle coordinates in the current frame.

    X = fltarr(numberOfParticles,1)
    Y = fltarr(numberOfParticles,1)

    for particle=0,numberOfParticles-1 do begin
      X(particle,0) = Xtemp(0,particle)
      Y(particle,0) = Ytemp(0,particle)
    endfor

    ;The following 'MAX' function is peculiar to IDL.
    ;It is used below to find the maximum value of the specified matrix.
    ;The following 'MIN' function is peculiar to IDL.
    ;It is used below to find the minimum value of the specified matrix.
    ;   xlen = max(X) - min(X)
    ;   ylen = max(Y) - min(Y)

    xlen = forceXlen
    ylen = forceYlen

    ratio = DOUBLE(ylen)/DOUBLE(xlen)

    ;The following procedures 'DEVICE', 'LOADCT' and 'WINDOW' are used to set up graphic display options.
    device, retain=2, decomposed =0
    loadct, 39
    !p.color = 0
    !p.background = 255
    window, 2, xsize = screenWidth, ysize =screenWidth*ratio

    ;The following 'PLOT' routine is perculiar to IDL.
    ;The PLOT procedure draws graphs of vector arguments.
    ;Below we use 'PLOT' for plotting particle positions.

    plot,X,Y,psym=3, isotropic=1, xrange = [0,forceXlen], yrange = [0,forceYlen],xstyle = 1, ystyle=1,title = 'polygon construction for frame' + string(frameNumber)

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
    file = image_file+string(frameNumber)+'.polygon_construction'+string(frameNumber)+'.tif'

    Write_Tiff, file, image24, 1

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

      ;    polyfill, Xcs,Ycs, color=140,linestyle = 0,thick = 0.5 ;color = 140 - green
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
      ;    polyfill, Xcs,Ycs, color=250,linestyle = 0,thick = 0.5 ;color = 250 - red
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
      ;    polyfill, Xcs,Ycs, color=200,linestyle = 0,thick = 0.5 ;color = 200 - yellow
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
      ;    polyfill, Xcs,Ycs, color=215,linestyle = 0,thick = 0.5 ;color = 215 - orange
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