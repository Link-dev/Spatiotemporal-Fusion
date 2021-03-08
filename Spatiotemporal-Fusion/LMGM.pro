;---------------------------------------------------------------------------------
;Rao Y, Zhu X, Chen J, et al. An improved method for producing high spatial-resolution NDVI time series datasets 
;with multi-temporal MODIS NDVI data and Landsat TM/ETM+ images[J]. Remote Sensing, 2015, 7(6): 7865-7891.
;
;This code is developed for remote sensing image spatiotemporal fusion.
;
;Input: Coarse image at t1 and fine image at t1
;Output: fine image at t2
;Note: this code is implemented based on the single pair of Landsat and MODIS at t1 with a MODIS image at t2,
;but he original NDVI_LMGM could deal with the multiple pairs of Landsat and MODIS images.
;
;
;parameter:
;     in the head of the main function
;
;Developed by Junxiong Zhou, email: zjxrs2018@mail.bnu.edu.cn
;
;Update history
;09/01/2018   first commit
;03/08/2021   revised by Junxiong
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
  winS=3                       ;set the half window size
  min_class=4.0                ;set the estimated minimum and maximum number of classes
  max_class=6.0
  DN_min=0.0                   ;set the range of DN value of the image,If byte, 0 and 255
  DN_max=10000.0
  sf=16.0                      ;set the scale factor, it is integer=coarse resolution/fine resolution, e.g., 480/30=16
  
  ;------------------------------------------------------------------------
  Basename = 'D:\'
  ;open the first coarse NDVI at t1
  coarseFile1 = Basename + 'MOD_2001_10_07'
  GetData,ImgData=coarse1,FileName = coarseFile1, Fid=fid1
  envi_file_query,fid1,ns=nsc, nl=nlc, nb=nbc, dims=dims_coarse
  
  ;open the second coarse NDVI at t2
  coarseFile2 = Basename + 'MOD_2002_01_04'
  GetData,ImgData=coarse2,FileName = coarseFile2, Fid=fid2

  ;open the fine img at t1
  fineFile = Basename + 'LT_2001_10_07'
  GetData,ImgData=fine,FileName = fineFile, Fid=fid3
  envi_file_query,fid3,ns=ns, nl=nl, nb=nb, dims=dims_fine
  
  OutName = fineFile + '_LMGM'; Results Name
  
  ;open the fine img for classification
  classMap = fineFile + '_ISODATA'
  if ~file_test(classMap) then begin
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
    out_name=fineFile+'_ISODATA'
    ENVI_DOIT, 'class_doit', fid=fid3, pos=indgen(nb), dims=dims_fine, $
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
  endif
  GetData,ImgData=classR,FileName = classMap, Fid=fid4
  nc = max(classR) - min(classR) + 1  ;Number of class

  ;calculate the Abundance, this step could be optimized by matrix calculation
  coarseAbun = dblarr(nsc, nlc, nc) ;Abundance of coarse pixels
  for i = 0, nsc-1 do begin
    for j = 0, nlc-1 do begin
      fineTmp = ClassR[i*sf:(i+1)*sf-1, j*sf:(j+1)*sf-1]

      for k = 0, nc-1 do begin
        tmp = where(fineTmp eq k+1, num_ic)
        coarseAbun[i,j,k] = num_ic/(sf*sf)
      endfor
    endfor
  endfor

  ;constraints
  coarseChange = coarse2 - coarse1 ;coarse temporal change
  stdK_coarse=stddev(coarseChange)
  minK_coarse=min(coarseChange)
  maxK_coarse=max(coarseChange)
  
  ;Unmix the coarse temporal change
  finePre = dblarr(ns, nl, nb)
  for bands = 0, nbc-1 do begin
    for i = 0, nsc-1 do begin
      for j = 0, nlc-1 do begin
        ai=max([0,i-winS])                         ; the MODIS window location
        bi=min([nsc-1,i+winS-1])
        aj=max([0,j-winS])
        bj=min([nlc-1,j+winS-1])

        n_N = (bi-ai+1) * (bj-aj+1)
        b = reform(coarseChange[ai:bi, aj:bj, bands], n_N)*1.0
        A = reform(coarseAbun[ai:bi, aj:bj, *], n_N, nc)

        xub=fltarr(nc,1) + min([maxK_coarse+stdK_coarse, DN_max-DN_min])
        xlb=fltarr(nc,1) + max([minK_coarse-stdK_coarse, DN_min-DN_max])
        c=fltarr(1,nc)+1
        bc=[n_N*(maxK_coarse+stdK_coarse)]
        contype=[1]
        result=IMSL_LINLSQ(b, A, c, bc, bc, contype, Xlb = xlb, Xub = xub)       ; Constrained Least Square
        num_nan = finite(result)
        result[where(num_nan eq 0, /null)] = 1.0*(DN_max-DN_min)         ;set the NAN to the range of the data

        tmpClass = classR[i*sf:(i+1)*sf-1, j*sf:(j+1)*sf-1]
        patches = dblarr(sf, sf)

        for m=1, nc do begin
          tmp = where(tmpClass eq m, num_ic)
          if num_ic ne 0 then begin
            patches[tmp] = result[m-1]
          endif
        endfor
        finePre[i*sf:(i+1)*sf-1, j*sf:(j+1)*sf-1, bands] = patches
      endfor
    endfor
  endfor
  finePre = fine + finePre

  map_info = envi_get_map_info(fid = fid3)
  Envi_Write_Envi_File, finePre, Out_Name = OutName, r_fid=fid_temp, ns = ns, nl = nl, nb = nb, MAP_INFO=map_info

;  envi_file_mng, id = fid1, /remove
;  envi_file_mng, id = fid2, /remove
;  envi_file_mng, id = fid3, /remove
;  envi_file_mng, id = fid4, /remove
;  envi_file_mng, id = fid_temp, /remove

  print, 'time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'

end
