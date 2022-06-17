FUNCTION plot_polygon, corename, N_vert, R, colorNo

  seconds = STRING(SYSTIME(/seconds),FORMAT='(I18)')
  filename = STRCOMPRESS(coreName, /REMOVE_ALL)

  filename_eps = filename+ STRCOMPRESS('_' + seconds + '.eps', /REMOVE_ALL)
  
  cam_resol = 0.0254303 ;mm/px
  plot_ybegin = 600;
  plot_yend = 1000;
  plot_xbegin = 0
  plot_xend = 1000
  screenWidth = 1600
  ratio = DOUBLE(plot_yend - plot_ybegin)/DOUBLE(plot_xend - plot_xbegin)
s = get_polygon(N_vert, R)

plot_x_offset = (plot_xend + plot_xbegin) / 2.0d
plot_y_offset = (plot_yend + plot_ybegin) / 2.0d

X = s.X + plot_x_offset
Y = s.Y + plot_y_offset

X_mm =  px_to_mm(X, plot_xbegin, cam_resol)
Y_mm =  px_to_mm(Y, plot_ybegin, cam_resol)
plot_xbegin_mm = px_to_mm(plot_xbegin, plot_xbegin, cam_resol)
plot_xend_mm = px_to_mm(plot_xend, plot_xbegin, cam_resol)
plot_ybegin_mm = px_to_mm(plot_ybegin, plot_ybegin, cam_resol)
plot_yend_mm = px_to_mm(plot_yend, plot_ybegin, cam_resol)


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

plot,X_mm,Y_mm, isotropic=1, $
  xrange = [plot_xbegin_mm,plot_xend_mm], yrange = [plot_ybegin_mm,plot_yend_mm], $
  xstyle = 1, $
  ystyle=1,title = 'Voronoi map', charsize = 2.0, thick = 8.0, charthick = 2, $
  /NODATA

POLYFILL, X_mm, Y_mm, COLOR = colorNo, linestyle = 0,thick = 0.5 & $
  oplot, X_mm, Y_mm
oplot, SHIFT(X_mm,1), SHIFT(Y_mm,1)

set_plot, originalDevice


;--------------------
;--------------------
originalDevice = !d.name
set_plot, 'ps'
loadct, 39
!p.color = 0
!p.background = 255
device, filename = filename_eps
;device, xsize= 1.5 * (3 + 3/8), ysize=  ratio * 1.5 * ((3 + 3/8)), /inches
device, FONT_SIZE = 12, /TIMES
device, color=1, bits_per_pixel=24
device, /encapsulated



plot,X_mm,Y_mm, isotropic=1, $
  xrange = [plot_xbegin_mm,plot_xend_mm], yrange = [plot_ybegin_mm,plot_yend_mm], $
  xstyle = 1, $
  ystyle=1,title = ' ', charsize = 2.0, thick = 8.0, charthick = 2, $
  /NODATA

POLYFILL, X_mm, Y_mm, COLOR = colorNo, linestyle = 0,thick = 0.5 & $
  oplot, X_mm, Y_mm
oplot, SHIFT(X_mm,1), SHIFT(Y_mm,1)

device, /close_file
set_plot, originalDevice

DEVICE, DECOMPOSED=old_decomposed

is_success = 1

RETURN, is_success

END