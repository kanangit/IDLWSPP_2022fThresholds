PRO driver_2022fThresholds_polyfill

coreName = 'triangle'

  datapath = 'e:\OneDrive - University of Iowa\bDocs\expAnalysisBackup\c_2022_polygons\c_20220615polygons'
  CD, datapath
  CD, 'outputs'
  
    seconds = STRING(SYSTIME(/seconds),FORMAT='(I18)')
  filename = STRCOMPRESS(coreName, /REMOVE_ALL)

  filename_eps = filename+ STRCOMPRESS('_' + seconds + '.eps', /REMOVE_ALL)
  
  originalDevice = !d.name
  set_plot, 'ps'
  device, filename=filename_eps
  ;device, xsize= 1.5 * (3 + 3/8), ysize=  ratio * 1.5 * ((3 + 3/8)), /inches
  device, FONT_SIZE = 12, /TIMES
  device, color=1, bits_per_pixel=24
  device, /encapsulated

  ; Create the vectors of X and Y values:
  X = [30, 100, 100, 30] & Y = [30, 30, 100, 100]

  ; Fill the polygon with color index 175:
  POLYFILL, X, Y, , COLOR = 140, linestyle = 0,thick = 0.5
  device, /close_file
  set_plot, originalDevice
END