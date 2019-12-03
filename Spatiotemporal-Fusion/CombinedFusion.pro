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

pro CombinedFusionByLMGM

  t0=systime(1)                  ;the initial time of program running

  ;please set the following parameters
  ;----------------------------------------------------------------------
  winS=3                       ;set the half window size
  scale_factor=16.0              ;set the scale factor, it is integer=coarse resolution/fine resolution, e.g., 480/30=16

  Basename = 'G:\'
  
  ;open the first coarse NDVI
  CoarseNDVI1 = Basename + 'T1_480_NDVI'
  GetData,ImgData=coarse1,FileName = CoarseNDVI1, Fid=fid1
  envi_file_query,fid1,ns=ns_coarse,nl=nl_coarse,nb=nb_coarse,dims=dims_coarse
  
  ;open the second coarse NDVI
  CoarseNDVI2 = Basename + 'T2_480_NDVI'
  GetData,ImgData=coarse2,FileName = CoarseNDVI2, Fid=fid2
  
  ;open the third coarse NDVI
  CoarseNDVI3 = Basename + 'T3_480_NDVI'
  GetData,ImgData=coarse3,FileName = CoarseNDVI3, Fid=fid3
  
  CoarseChange1 = abs(coarse1-coarse2)
  CoarseChange2 = abs(coarse3-coarse2)
  
  ;open the first result
  ResultNDVI1 = Basename + 'Result1'
  GetData,ImgData=result1,FileName = ResultNDVI1, Fid=fid4
  envi_file_query,fid4,ns=ns_fine,nl=nl_fine,nb=nb_fine,dims=dims_fine
  
  ;open the second result
  ResultNDVI2 = Basename + 'Result2'
  GetData,ImgData=result2,FileName = ResultNDVI2, Fid=fid5
  
  ;Final Result
  finePre = dblarr(ns_fine, nl_fine)

  ;Combined
  for i = 0, ns_coarse-1 do begin
    for j = 0, nl_coarse-1 do begin
      ai=max([0,i-winS])                         ; the MODIS window location
      bi=min([ns_coarse-1,i+winS])               ; 预设的窗口大小为3*3，w取值为1
      aj=max([0,j-winS])
      bj=min([nl_coarse-1,j+winS])
      
      
      tmp1 = 1.0/total(CoarseChange1[ai:bi, aj:bj])
      tmp2 = 1.0/total(CoarseChange2[ai:bi, aj:bj])
      
      w1 = tmp1 / (tmp1+tmp2)
      w2 = tmp2 / (tmp1+tmp2)
      
      finePre[i*scale_factor:(i+1)*scale_factor-1, j*scale_factor:(j+1)*scale_factor-1] = w1*result1[i*scale_factor:(i+1)*scale_factor-1, j*scale_factor:(j+1)*scale_factor-1] $
         + w2*result2[i*scale_factor:(i+1)*scale_factor-1, j*scale_factor:(j+1)*scale_factor-1]
      
    endfor
  endfor
  
  map_info = envi_get_map_info(fid = fid4)
  OutName = Basename + 'ResultCombined'
  Envi_Write_Envi_File, finePre, Out_Name = OutName, r_fid=fid_temp, ns = ns_fine, nl = nl_fine, nb = 1, MAP_INFO=map_info
  
  envi_file_mng, id = fid1, /remove
  envi_file_mng, id = fid2, /remove
  envi_file_mng, id = fid3, /remove
  envi_file_mng, id = fid4, /remove
  envi_file_mng, id = fid5, /remove
  
end