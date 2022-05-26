FUNCTION read_pulse_pos, path

  numParams = N_Params()
  IF (numParams EQ 0) THEN BEGIN
    filename = DIALOG_PICKFILE(/READ, FILTER = '*.csv')
  ENDIF
  IF (numParams EQ 1) THEN BEGIN
    filename = path
  ENDIF

  returnStructur = 0
  ;check if the template file exist
  existance = FILE_TEST('read_pulse_pos_test.sav')
  ; if it does not exist, ask user to create one:
  if (existance EQ 0) then begin
    rTemplate = ASCII_TEMPLATE(filename)
    SAVE, rTemplate, FILENAME='read_pulse_pos_test.sav'
    existance = 1
  endif
  ;obtain the file template rTemplate stored in the
  ;file 'read_pulse_pos_test.sav':
  RESTORE, 'read_pulse_pos_test.sav'
  ;input the data
  inputStructur = READ_ASCII(filename, template=rTemplate)

  RETURN, inputStructur
END