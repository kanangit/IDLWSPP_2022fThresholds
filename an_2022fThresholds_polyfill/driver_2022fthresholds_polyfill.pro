PRO driver_2022fThresholds_polyfill



  datapath = 'e:\OneDrive - University of Iowa\bDocs\expAnalysisBackup\c_2022_polygons\c_20220617polygons\'
  CD, datapath ;khuj tebe local history 
  CD, 'outputs'

  aWignerZeiss = 16.337134d; px
  R = aWignerZeiss


  coreName = 'vrn_hexagons'
  coord_number = 6
  COLOR = 140
  plot_suc = plot_polygon(coreName, coord_number, R, COLOR)
  
  coreName = 'vrn_pentagons'
  coord_number = 5
  COLOR = 40
  plot_suc = plot_polygon(coreName, coord_number, R, COLOR)
  stop
  
  coreName = 'vrn_heptagons'
  coord_number = 7
  COLOR = 215
  plot_suc = plot_polygon(coreName, coord_number, R, COLOR)
  
  
  coreName = 'vrn_quadrugons'
  coord_number = 4
  COLOR = coord_number*30 + 20
  plot_suc = plot_polygon(coreName, coord_number, R, COLOR)
 
  
  
  coreName = 'vrn_manygons'
  coord_number = 3
  COLOR = coord_number*30 + 20
  print, COLOR
  plot_suc = plot_polygon(coreName, coord_number, R, COLOR)
  stop

END