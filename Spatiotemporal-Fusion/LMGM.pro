;---------------------------------------------------------------------------------
;Rao Y, Zhu X, Chen J, et al. An improved method for producing high spatial-resolution NDVI time series datasets 
;with multi-temporal MODIS NDVI data and Landsat TM/ETM+ images[J]. Remote Sensing, 2015, 7(6): 7865-7891.
;
;This code is developed for remote sensing image spatiotemporal fusion.
;
;Input: Coarse image at t1 and , fine image at t1
;Output: fine image at t2
;
;parameter:
;     in the head of the main function
;
;Developed by Zhou Junxiong, email: zjxrs2018@mail.bnu.edu.cn
;
;Update history
;09/01/2018   first commit
;---------------------------------------------------------------------------------

; -------------------------------------------------------------------------------------------------------------------
; universal function: function for open the file
; -------------------------------------------------------------------------------------------------------------------
Pro GetData,ImgData = ImgData,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,$
  FileName = FileName,Map_info = map_Info, Fid = Fid, dims = dims
  Envi_Open_File,FileName,R_Fid = Fid
  Envi_File_Query,Fid,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type
  map_info = envi_get_map_info(fid=Fid)
  dims = [-1,0,ns - 1 ,0,nl - 1]
  case Data_Type Of
    1:ImgData = BytArr(ns,nl,nb)    ;  BYTE  Byte
    2:ImgData = IntArr(ns,nl,nb)    ;  INT  Integer
    3:ImgData = LonArr(ns,nl,nb)    ;  LONG  Longword integer
    4:ImgData = FltArr(ns,nl,nb)    ;  FLOAT  Floating point
    5:ImgData = DblArr(ns,nl,nb)    ;  DOUBLE  Double-precision floating
    6:ImgData = COMPLEXARR(ns,nl,nb); complex, single-precision, floating-point
    9:ImgData = DCOMPLEXARR(ns,nl,nb);complex, double-precision, floating-point
    12:ImgData = UINTARR(ns,nl,nb)   ; unsigned integer vector or array
    13:ImgData = ULONARR(ns,nl,nb)   ;  unsigned longword integer vector or array
    14:ImgData = LON64ARR(ns,nl,nb)   ;a 64-bit integer vector or array
    15:ImgData = ULON64ARR(ns,nl,nb)   ;an unsigned 64-bit integer vector or array
  EndCase
  For i = 0,nb-1 Do Begin
    Dt = Envi_Get_Data(Fid = Fid,dims = dims,pos=i)
    ImgData[*,*,i] = Dt[*,*]
  EndFor
End
;-------------------------------------------------------------------------------------------------------


pro LMGM
  t0=systime(1)                  ;the initial time of program running

  ;please set the following parameters
  ;----------------------------------------------------------------------
  winS=1                       ;set the half window size
  ;winO=1                       ;set the overlap of window
  min_class=4.0                ;set the estimated minimum and maximum number of classes
  max_class=6.0
  DN_min=0.0                   ;set the range of DN value of the image,If byte, 0 and 255
  DN_max=10000.0
  scale_factor=16.0              ;set the scale factor, it is integer=coarse resolution/fine resolution, e.g., 480/30=16
  
  Basename = 'G:\'
  ;------------------------------------------------------------------------

  ;open the first coarse NDVI
  CoarseNDVI1 = Basename + 'T3_480_NDVI'
  GetData,ImgData=coarse1,FileName = CoarseNDVI1, Fid=fid1
  envi_file_query,fid1,ns=ns_coarse,nl=nl_coarse,nb=nb_coarse,dims=dims_coarse
  
  ;open the second coarse NDVI
  CoarseNDVI2 = Basename + 'T2_480_NDVI'
  GetData,ImgData=coarse2,FileName = CoarseNDVI2, Fid=fid2

  ;open the fine img
  FineNDVI = Basename + 'T3_NDVI'
  GetData,ImgData=fine,FileName = FineNDVI, Fid=fid3
  envi_file_query,fid3,ns=ns_fine,nl=nl_fine,nb=nb_fine,dims=dims_fine
  
  ;open the fine img for classification
  FineImg = Basename + 'L7_2001_11_01_flaash_CIA_patch_img'
  GetData,ImgData=Img,FileName = FineImg, Fid=fid4
  
  
  ;get spectral classes from fine resolution image at fine by isodata
  ;parameter of isodata
  CHANGE_THRESH = .05
  NUM_CLASSES = max_class
  ITERATIONS = 20
  ISO_MERGE_DIST = 0.05*DN_max
  ISO_MERGE_PAIRS = 2
  ISO_MIN_PIXELS = 200
  ISO_SPLIT_SMULT = 1
  ISO_SPLIT_STD = 0.05*DN_max
  MIN_CLASSES = min_class
  out_bname = 'IsoData'
  out_name=FineImg+'class_ISODATA'
  ENVI_DOIT, 'class_doit', fid=fid4, pos=indgen(nb_fine), dims=dims_fine, $
    out_bname=out_bname, out_name=out_name, method=4, $
    r_fid=r_fid, $
    NUM_CLASSES = NUM_CLASSES, $
    ITERATIONS = ITERATIONS, $
    CHANGE_THRESH = CHANGE_THRESH, $
    ISO_MERGE_DIST = ISO_MERGE_DIST, $
    ISO_MERGE_PAIRS = ISO_MERGE_PAIRS, $
    ISO_MIN_PIXELS = ISO_MIN_PIXELS, $
    ISO_SPLIT_SMULT = ISO_SPLIT_SMULT, $
    ISO_SPLIT_STD = ISO_SPLIT_STD, $
    MIN_CLASSES = MIN_CLASSES

  ;Result of classification
  ClassR = Envi_Get_Data(Fid = r_fid,dims = dims_fine,pos=0)

  ;Number of class
  classNumber = max(classR) - min(classR) + 1

  ;Abundance of Coarse pixel
  coarseAbun = dblarr(ns_coarse, nl_coarse, classNumber)

  ;calculate the Abundance
  ;how to optimize?
  for i = 0, ns_coarse-1 do begin
    for j = 0, nl_coarse-1 do begin
      fineTmp = ClassR[i*scale_factor:(i+1)*scale_factor-1, j*scale_factor:(j+1)*scale_factor-1]

      for k = 0, classNumber-1 do begin
        tmp = where(fineTmp eq k+1, num_ic)
        coarseAbun[i,j,k] = num_ic / (scale_factor*scale_factor)
      endfor
    endfor
  endfor

  ;Unmixing result
  finePre = dblarr(ns_fine, nl_fine)

  ;constrained
  coarseChange = coarse2 - coarse1
  
  stdK_coarse=stddev(coarseChange)
  minK_coarse=min(coarseChange)
  maxK_coarse=max(coarseChange)

  ;Unmixing
  for i = 0, ns_coarse-1 do begin
    for j = 0, nl_coarse-1 do begin
      ai=max([0,i-winS])                         ; the MODIS window location
      bi=min([ns_coarse-1,i+winS])             
      aj=max([0,j-winS])
      bj=min([nl_coarse-1,j+winS])


      n_N = (bi-ai) * (bj-aj+1)
      b = reform(coarseChange[ai:bi, aj:bj], n_N)*1.0
      A = reform(coarseAbun[ai:bi, aj:bj, *], n_N, classNumber)

      xub=fltarr(classNumber,1) + maxK_coarse+stdK_coarse
      xlb=fltarr(classNumber,1) + minK_coarse-stdK_coarse
      c=fltarr(1,classNumber)+1
      bc=[n_N*(maxK_coarse+stdK_coarse)]
      contype=[1]
      result=IMSL_LINLSQ(b, A, c, bc, bc, contype, Xlb = xlb, Xub = xub)       ; Constrained Least Square
      num_nan = finite(result)
      result[where(num_nan eq 0, /null)] = 1.0         ;set the NAN to 10000.0

      tmpClass = classR[i*scale_factor:(i+1)*scale_factor-1, j*scale_factor:(j+1)*scale_factor-1]
      patches = dblarr(scale_factor, scale_factor)
      
      for m=1, classNumber do begin
        tmp = where(tmpClass eq m, num_ic)
        if num_ic ne 0 then begin
          patches[tmp] = result[m-1]
        endif
      endfor
      finePre[i*scale_factor:(i+1)*scale_factor-1, j*scale_factor:(j+1)*scale_factor-1] = patches

    endfor
  endfor
  
  finePre = fine + finePre
  
  ; revise the outlier value in prediction result
  ;open the second coarse NDVI (the same size with fine NDVI)
  CoarseNDVI2 = Basename + 'T2_480_NDVI_biliNear'
  GetData,ImgData=coarse3,FileName = CoarseNDVI2, Fid=fid5
  ind = where((finePre gt 1.0) or (finePre lt -1.0), num_i)
  if num_i gt 0 then begin
    finePre[ind] = coarse3[ind]
  endif

  map_info = envi_get_map_info(fid = fid3)
  OutName = Basename + 'Result2'
  Envi_Write_Envi_File, finePre, Out_Name = OutName, r_fid=fid_temp, ns = ns_fine, nl = nl_fine, nb = 1, MAP_INFO=map_info

  envi_file_mng, id = fid1, /remove
  envi_file_mng, id = fid2, /remove
  envi_file_mng, id = fid3, /remove
  envi_file_mng, id = fid4, /remove
  envi_file_mng, id = fid5, /remove
  envi_file_mng, id = r_fid, /remove
  ;envi_file_mng, id = fid_temp, /remove

  print, 'time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'

end
