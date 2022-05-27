FUNCTION oldplot_defects, XmyFrame, YmyFrame, shockPos, filename

  DEVICE, GET_DECOMPOSED=old_decomposed
  device, retain=2, decomposed =0
  loadct, 39
  !p.color = 0
  !p.background = 255

  forceXlen = 1100
  forceYlen = 1200
  screenWidth = 1200
  ratio = DOUBLE(forceYlen)/DOUBLE(forceXlen)

  plot_ybegin = 600;
  plot_yend = 1000;
  plot_xbegin = 0
  plot_xend = 1100

  preshock_offset = 50.0
  postshock_offset = 50.0

  shockfront_x = [shockPos, shockPos]
  shockFront_y = [plot_ybegin, plot_yend]

  postshock_right_border = shockPos - preshock_offset
  preshock_left_border = shockPos + postshock_offset

  forceXlen = 1100
  forceYlen = 1200
  screenWidth = forceXlen
  ratio = DOUBLE(plot_yend - plot_ybegin)/DOUBLE(plot_xend - plot_xbegin)

  window, 2, xsize = screenWidth, ysize =screenWidth*ratio



  ; Triangulate it:
  TRIANGULATE, XmyFrame, YmyFrame, tr, CONN=C
  ;  stop

  N = N_ELEMENTS(XmyFrame);
  ;N = 1400




  ;  stop
  ;plot,XmyFrame,YmyFrame,psym=3, isotropic=1, xrange = [0,forceXlen], yrange = [0,forceYlen],xstyle = 1, ystyle=1,title = 'Voronoi map for unperturbed liquid'
  plot,XmyFrame,YmyFrame,psym=3, isotropic=1, xrange = [plot_xbegin,plot_xend], yrange = [plot_ybegin,plot_yend],xstyle = 1, $
    ystyle=1,title = 'Voronoi map', charsize = 2.0, thick = 8.0, charthick = 2, $
    /NODATA


  voronXLeftMargin = min(XmyFrame)
  voronXRightMargin = max(XmyFrame)
  voronYLeftMargin = min(YmyFrame)
  voronYRightMargin = max(YmyFrame)
  num_hexagons = LONG(0)
  num_pentagons = LONG(0)
  num_heptagons = LONG(0)
  num_othergons = LONG(0)
  preshock_num_hexagons = LONG(0)
  preshock_num_pentagons = LONG(0)
  preshock_num_heptagons = LONG(0)
  preshock_num_othergons = LONG(0)
  postshock_num_hexagons = LONG(0)
  postshock_num_pentagons = LONG(0)
  postshock_num_heptagons = LONG(0)
  postshock_num_othergons = LONG(0)
  FOR I=0, N-1 DO BEGIN & $
    ; Get the ith polygon:
    VORONOI, XmyFrame, YmyFrame, I, C, Xp, Yp & $
    ; Draw it:
    if (MAX(XP) LE voronXRightMargin AND MAX(YP) LE voronYRightMargin AND MIN(XP) GE voronXLeftMargin AND MIN(YP) GE voronYLeftMargin) then begin
    ;      POLYFILL, Xp, Yp, COLOR = 17 - LONG(DOUBLE((N_ELEMENTS(Xp)) - 10.0d/3.0d) * 3) & $
    coordNumber = N_ELEMENTS(XP)
    switch (coordNumber) of
      6: begin
        num_hexagons = num_hexagons + 1
        if (MIN(Xp) GE preshock_left_border) THEN $
          preshock_num_hexagons = preshock_num_hexagons + 1
        if (MAX(Xp) LE postshock_right_border) THEN $
          postshock_num_hexagons = postshock_num_hexagons + 1
        POLYFILL, Xp, Yp, COLOR = 140, linestyle = 0,thick = 0.5 & $
          oplot, Xp, Yp
        break
      end
      5: begin
        num_pentagons = num_pentagons + 1
        if (MIN(Xp) GE preshock_left_border) THEN $
          preshock_num_pentagons = preshock_num_pentagons + 1
        if (MAX(Xp) LE postshock_right_border) THEN $
          postshock_num_pentagons = postshock_num_pentagons + 1
        POLYFILL, Xp, Yp, COLOR = 40, linestyle = 0,thick = 0.5 & $
          oplot, Xp, Yp
        break
      end
      7: begin
        num_heptagons = num_heptagons + 1
        if (MIN(Xp) GE preshock_left_border) THEN $
          preshock_num_heptagons = preshock_num_heptagons + 1
        if (MAX(Xp) LE postshock_right_border) THEN $
          postshock_num_heptagons = postshock_num_heptagons + 1
        POLYFILL, Xp, Yp, COLOR = 215, linestyle = 0,thick = 0.5 & $
          oplot, Xp, Yp
        break
      end
      else: begin
        num_othergons = num_othergons + 1
        if (MIN(Xp) GE preshock_left_border) THEN $
          preshock_num_othergons = preshock_num_othergons + 1
        if (MAX(Xp) LE postshock_right_border) THEN $
          postshock_num_othergons = postshock_num_othergons + 1
        POLYFILL, Xp, Yp, COLOR = (coordNumber*30 + 20  ), linestyle = 0,thick = 0.5 & $
          oplot, Xp, Yp
      end
    endswitch



  endif
  ;    print, "I = ", I
  ;    wait, 0.01


endfor

oplot,XmyFrame,YmyFrame,psym=3

num_defects = num_othergons + num_pentagons + num_pentagons
num_total = num_hexagons + num_defects
defect_ratio = DOUBLE(num_defects) / DOUBLE(num_total)

preshock_num_defects = preshock_num_othergons + preshock_num_pentagons + preshock_num_pentagons
preshock_num_total = preshock_num_hexagons + preshock_num_defects
preshock_defect_ratio = DOUBLE(preshock_num_defects) / DOUBLE(preshock_num_total)

postshock_num_defects = postshock_num_othergons + postshock_num_pentagons + postshock_num_pentagons
postshock_num_total = postshock_num_hexagons + postshock_num_defects
postshock_defect_ratio = DOUBLE(postshock_num_defects) / DOUBLE(postshock_num_total)

XYOUTS, -100, 525, STRCOMPRESS('postshock defect ratio = ' + STRING(postshock_defect_ratio))

;plotting the shock front position:
plots, shockfront_x, shockFront_y
; Reset device paramters

;save the file as tiff:

scaled = TVRD()
TVLCT, r, g, b, /Get
s = Size(scaled, /Dimensions)
image24 = BytArr(3, s[0], s[1])
image24[0,*,*] = r[scaled]
image24[1,*,*] = g[scaled]
image24[2,*,*] = b[scaled]
image24 = TVRD(True=1)
image24 = Reverse(image24, 3)


Write_Tiff, filename, image24, 1;, XRESOL = 300, YRESOL = 300

DEVICE, DECOMPOSED=old_decomposed
;LOADCT, 0



return, defect_ratio
END