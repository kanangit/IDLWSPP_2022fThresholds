; v2.2
;
; NAME:
;       readImageJK
;
; PURPOSE:
;       This function reads the output file "position.txt" produced by
;       the ImageJ macros "Particle_Identification". "position.txt"
;       file must consist of the columns of data in the following
;       order: particle number, frame label, particle area,
;       paricle x-coordinate, particle y-coordinate, slice number.
;       It is assumbed that columns are separated by commas
;       and arbitrary number of spaces AND/OR tabs.
;
; CALLING SEQUENCE:
;       Result = readImageJK()
;
; INPUTS
;       The fuction doesn't have input variables.
;
; OUTPUTS
;       A structure with the following tags:
;         .iParticle  - an array of particle numbers
;         .iFrame     - an array of frame numbers
;         .X          - an array of x-coordinate of the particles
;         .Y          - an array of y-coordinate of the particles
;         .area       - an array with area of the particles
;         .error      - an array with errors in particle positions
;
;
; PROCEDURE:
;
;
; MODIFICATION HISTORY:
;           Written by:  Anton Kananovich, April 2016
;           Modified:    Anton Kananovich, May 2018.
;                         Optimized for memory economy. Now works
;                         correctly with framelables which contain file
;                         extension (for example, "fubar4834.tiff").
;                        Anton Kananovich, June 2018
;                         added a keyword lowmem. Set the keywork if your computer
;                         has insufficient memory
;                        Anton Kananovich, July 2018
;                         fixed the bug, when each element of the
;                         /lowmem version of the array
;                         returnStructur.iParticle was larger by unity
;                         than that of the usual version of the
;                         function.
;                         The iParticle and iFrame not return exact same numbers
;                         as in the file, without subtraction of 1 (unity)
;           Modified:    Anton Kananovich, July 2018
;                         added the 'path' input variable. Without the
;                         variable supplied, the function will work as
;                         before. With the variable supplied, the
;                         will open and process the file designated
;                         in the variable.
;-

function readImageJK, path, lowmem = lowmem
  numParams = N_Params()
  IF (numParams EQ 0) THEN BEGIN
    filename = DIALOG_PICKFILE(/READ, FILTER = '*.txt')
  ENDIF
  IF (numParams EQ 1) THEN BEGIN
    filename = path
  ENDIF

  if (KEYWORD_SET(lowmem)) then begin
    returnStructur = readImageJKSlow(filename)
  endif else begin
    returnStructur = readImageJKFast(filename)
  endelse
  RETURN, returnStructur
END
function readImageJKFast, filename

  nlines = FILE_LINES(filename) ;get the number of lines in the file
  sarr = STRARR(nlines) ;initializing a string array to store the data read from the file
  GET_LUN, lun
  OPENR, lun, filename
  READF, lun, sarr ;read all data at once. no cycles. in IDL it will work very fast
  FREE_LUN, lun
  sarr=STRCOMPRESS(temporary(sarr),/REMOVE_ALL); remove all white spaces and tabs in the strings. now the fields are separated only by commas
  slist = STRSPLIT(temporary(sarr), ',', /EXTRACT); split the strings into fields using comma as a delimiter
  nelements = N_ELEMENTS(slist)
  alist = slist.ToArray() ;convert the output of the strsplit() function to array
  slist = 0; free the memory


  ;create a structure in which we store all the data our fucntion will return:
  returnStructur = {iParticle:LONARR(nelements),iFrame:LONARR(nelements), $
    area:DBLARR(nelements), X:DBLARR(nelements),Y:DBLARR(nelements),error:DBLARR(nelements)}
  ;populating the structure fields with the fields values we read from file. again - no cycles. Cycles slow down IDL:
  returnStructur.iParticle = LONG(alist[*,0])
  returnStructur.area = DOUBLE(alist[*,2])
  returnStructur.X = DOUBLE(alist[*,3])
  returnStructur.Y = DOUBLE(alist[*,4])
  returnStructur.error = DOUBLE(alist[*,5])
  frameLabel = alist[*,1] ; we don't store the values of this field in the return structure,
  ;because we need to convert it first

  for i = 0L, nelements-1 do begin
    pos = STREGEX(frameLabel[i], '([0-9]+)(\.|$)', length=len, /SUBEXPR) ;extract the
    ;framenumber out of the frame label string using regular expressions
    ;The regular expression '[0-9]+$' means
    ;"give me the end of the string, which contains digits only"
    returnStructur.iFrame[i] = LONG(STRMID(frameLabel[i],pos[1],len[1])); now we are ready to pass the
    ;frame number to the corresponding field of the return structure
  endfor

  return, returnStructur ;returning the result
end


FUNCTION readImageJKSlow, filename
  returnStructur = 0
  ;check if the template file exist
  existance = FILE_TEST('readImageJKtemplate.sav')
  ; if it does not exist, ask user to create one:
  if (existance EQ 0) then begin
    rTemplate = ASCII_TEMPLATE(filename)
    SAVE, rTemplate, FILENAME='readImageJKtemplate.sav'
    existance = 1
  endif
  ;obtain the file template rTemplate stored in the
  ;file 'readImageJKtemplate.sav':
  RESTORE, 'readImageJKtemplate.sav'
  ;input the data
  inputStructur = READ_ASCII(filename, template=rTemplate)
  ;determine number of elements
  nelements = N_ELEMENTS(inputStructur.X)
  ;create a structure in which we store all the data our fucntion will return:
  returnStructur = {iParticle:LONARR(nelements),iFrame:LONARR(nelements), $
    area:DBLARR(nelements), X:DBLARR(nelements),Y:DBLARR(nelements),error:DBLARR(nelements)}

  ;prepare to convert framelabel (string) to frame number (integer)

  for i = 0L, nelements-1 do begin
    ;use the regexp syntax to determine the position of the frame number:
    pos = STREGEX(inputStructur.Framelabel[i], '([0-9]+)(\.|$)', length=len, /SUBEXPR) ;extract the
    ;framenumber out of the frame label string using regular expressions
    ;The regular expression '[0-9]+$' means
    ;"give me the end of the string, which contains digits only"
    returnStructur.iFrame[i] = LONG(STRMID(inputStructur.Framelabel[i],pos[1],len[1])); now we are ready to pass the
    ;frame number to the corresponding field of the return structure
  endfor
  returnStructur.iParticle = inputStructur.iParticle
  returnStructur.area = inputStructur.area
  returnStructur.X = inputStructur.X
  returnStructur.Y = inputStructur.Y
  returnStructur.error = inputStructur.error

  inputStructur = 0 ; free the memory
  RETURN, returnStructur
END