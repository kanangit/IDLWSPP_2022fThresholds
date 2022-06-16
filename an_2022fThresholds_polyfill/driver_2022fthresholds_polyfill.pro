PRO driver_2022fThresholds_polyfill

  coreName = 'triangle'

  datapath = 'e:\OneDrive - University of Iowa\bDocs\expAnalysisBackup\c_2022_polygons\c_20220615polygons'
  CD, datapath ;khuj tebe local history
  CD, 'outputs'

  seconds = STRING(SYSTIME(/seconds),FORMAT='(I18)')
  filename = STRCOMPRESS(coreName, /REMOVE_ALL)

  filename_eps = filename+ STRCOMPRESS('_' + seconds + '.eps', /REMOVE_ALL)


  ;resol = 0.0254303 ;mm/px
  resol = 1

  plot_ybegin_mm = 0;
  plot_yend_mm = 120;
  plot_xbegin_mm = 0
  plot_xend_mm = 120

  screenWidth = 900
  ratio = DOUBLE(1)

  ; Create the vectors of X and Y values:
  X = [30, 100, 100, 30] & Y = [30, 30, 100, 100]

  ;--------------------
  ;--------------------
  originalDevice = !d.name
  set_plot, 'WIN'
  DEVICE, GET_DECOMPOSED=old_decomposed
  device, retain=2, decomposed =0
  loadct, 39
  !p.color = 0
  !p.background = 255
  window, 2, xsize = screenWidth, ysize =screenWidth*ratio

  plot,X*resol,Y*resol, isotropic=1, $
    xrange = [plot_xbegin_mm,plot_xend_mm], yrange = [plot_ybegin_mm,plot_yend_mm], $
    xstyle = 1, $
    ystyle=1,title = 'Voronoi map', charsize = 2.0, thick = 8.0, charthick = 2, $
    /NODATA

  POLYFILL, X*resol, Y*resol, COLOR = 140, linestyle = 0,thick = 0.5 & $
    oplot, X*resol, Y*resol
    oplot, SHIFT(X,1)*resol, SHIFT(Y,1)*resol

  set_plot, originalDevice


  ;--------------------
  ;--------------------
  originalDevice = !d.name
  set_plot, 'ps'
  loadct, 39
  !p.color = 0
  !p.background = 255
  device, filename=filename_eps
  ;device, xsize= 1.5 * (3 + 3/8), ysize=  ratio * 1.5 * ((3 + 3/8)), /inches
  device, FONT_SIZE = 12, /TIMES
  device, color=1, bits_per_pixel=24
  device, /encapsulated



  ; Fill the polygon with color index 175:

  ;POLYFILL, X, Y, COLOR = 140, /DEVICE

  plot,X*resol,Y*resol, isotropic=1, $
    xrange = [plot_xbegin_mm,plot_xend_mm], yrange = [plot_ybegin_mm,plot_yend_mm], $
    xstyle = 1, $
    ystyle=1,title = 'Voronoi map', charsize = 2.0, thick = 8.0, charthick = 2, $
    /NODATA

  POLYFILL, X*resol, Y*resol, COLOR = 140, linestyle = 0,thick = 0.5 & $
    oplot, SHIFT(X,1)*resol, SHIFT(Y,1)*resol


  device, /close_file
  set_plot, originalDevice

  DEVICE, DECOMPOSED=old_decomposed
END