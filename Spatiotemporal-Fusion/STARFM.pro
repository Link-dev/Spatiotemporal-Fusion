;---------------------------------------------------------------------------------
;Gao F, Masek J, Schwaller M, et al. On the blending of the Landsat and MODIS surface reflectance: Predicting daily Landsat surface reflectance[J]. 
;IEEE Transactions on Geoscience and Remote sensing, 2006, 44(8): 2207-2218.
;
;This code try to reproduce the STARFM with IDL.
;
;Input: Coarse image in t1 and t2
;Output: fine image in t1
;
;Note that:
; The size of fine image and coarse must be equal
; 
;
;parameter:
;     in the head of the main function
;
;Developed by Zhou Junxiong, email: zjxrs2018@mail.bnu.edu.cn
;
;Update history
;31/12/2018   first version
;03/12/2019   optimize codes
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

pro STARFM
  t0=systime(1)                  ;the initial time of program running

  ;please set the following parameters
  ;----------------------------------------------------------------------
  scale = 10000.0                       ;10000 or 1
  DN_min=0.0                            ;set the range of DN value of the image,If byte, 0 and 255
  DN_max=1.0*scale
  ratio=16.0                            ;set the resolution ratio, it is integer=coarse resolution/fine resolution, e.g., 480/30=16
  fine_uncertainty = 0.005*scale        ;uncertainty of fine image
  coarse_uncertainty = 0.005*scale      ;uncertainty of coarse image
  combined_uncertainty = sqrt(coarse_uncertainty^2+fine_uncertainty^2); 

  num_class = 40
  win=floor(ratio*1.5/2)                ;set the half window size, if this parameter is so big, the process will be so time-consuming
  ;------------------------------------------------------------------------

  Basename = 'E:\test\'
    
  ;open the first coarse NDVI
  CoarseNDVI1 = Basename + 'MODIS20041126_SIMU'
  GetData,ImgData=coarse1,FileName = CoarseNDVI1, Fid=fid1
  envi_file_query,fid1,ns=ns_coarse,nl=nl_coarse,nb=nb_coarse,dims=dims_coarse
   
  CoarseNDVI1 = Basename + 'MODIS20041212_SIMU'
  GetData,ImgData=coarse2,FileName = CoarseNDVI1, Fid=fid2
  
  ;open the fine img
  FineNDVI = Basename + 'Landsat20041126'
  GetData,ImgData=fine,FileName = FineNDVI, Fid=fid3
  envi_file_query,fid3,ns=ns_fine,nl=nl_fine,nb=nb_fine,dims=dims_fine
  
  ;Assume that difference in coarse is comparable in fine
  DFusion = fine + coarse2 - coarse1
  
  ;Spectral difference between coarse and fine
  SDistance = abs(coarse1 - fine)
  SDistance_uncerntainty = sqrt((Fine_uncertainty^2+coarse_uncertainty^2))
  
  ;Temporal difference between coarse1 and coarse2
  TDistance = abs(coarse2 - coarse1)
  TDistance_uncerntainty = sqrt((2*coarse_uncertainty^2))
  
  ;Result of prediction
  finePre = dblarr(ns_fine,nl_fine,nb_fine)
  
  ;compute the distance of each pixel in the window with the target pixel (integrate window)
  D_D_all=((win-indgen(win*2+1)#(intarr(1,win*2+1)+1))^2+$
    (win-(intarr(win*2+1)+1)#indgen(1,win*2+1))^2)^0.5
  D_D_all=reform(D_D_all,(win*2+1)*(win*2+1))
  
  ;similar test
  similar_th=fltarr(nb_fine)
  for iband=0,nb_fine-1,1 do begin
    similar_th[iband,0]=stddev(fine[*,*,iband])*2.0/num_class
  endfor
  
  ;Prediction
  for i = 0, ns_fine-1 do begin
    for j = 0, nl_fine-1 do begin
      ai=max([0,i-win])                         
      bi=min([ns_fine-1,i+win])               
      aj=max([0,j-win])
      bj=min([nl_fine-1,j+win])
      
      n_N = (bi-ai+1)*(bj-aj+1)
      ci=i-ai      ;location of target pixel
      cj=j-aj
      col_wind=indgen(bi-ai+1)#(intarr(1,bj-aj+1)+1)*1.0
      row_wind=(intarr(bi-ai+1)+1)#indgen(1,bj-aj+1)*1.0
      if ((bi-ai+1)*(bj-aj+1) lt (win*2.0+1)*(win*2.0+1)) then begin   ;not an integrate window
        D_D_cand=((ci-col_wind)^2+(cj-row_wind)^2)^0.5
      endif else begin
        D_D_cand=D_D_all      ;integrate window
      endelse
      SpatialDis=(1.0+D_D_cand/win)

      for k = 0, nb_fine-1 do begin
        ;Combined distance
        ; if input only contains one pair of Landsat and MODIS, then only use difference between MODIS and Landsat   (STARFM, Feng Gao)
        CDistance = alog(SDistance[ai:bi,aj:bj,k]+1)*SpatialDis+0.0001   ; *alog(TDistance[ai:bi,aj:bj,k]+1)
        
        ; very close (less than uncertainty), gets maximum weight
        ind = where(SDistance[ai:bi,aj:bj,k] le combined_uncertainty)
        CDistance[ind] = 1.0
        
        ;exclude some worse and not similar neighbor pixels
        index = where(abs(fine[ai:bi,aj:bj,k]-fine[i,j,k]) ge similar_th[k] or $
          SDistance[ai:bi,aj:bj,k] ge (max(SDistance[ci,cj,k]+SDistance_uncerntainty)) or $
          TDistance[ai:bi,aj:bj,k] ge (max(TDistance[ci,cj,k]+TDistance_uncerntainty)))
        CDistance[index] = -9999
        
        indcand = where(CDistance ne -9999)
        CDistance = CDistance[indcand]
        
        ;use a nomalized reverse distance as the weight function
        CDistance = 1.0/CDistance
        
        weight = CDistance / total(CDistance)
        DValue = (DFusion[ai:bi, aj:bj, k])[indcand]
        
        finePre[i,j,k] = total(weight * DValue)
      endfor
      
    endfor
  endfor

  map_info = envi_get_map_info(fid = fid3)
  OutName = Basename + 'STARFM_Gwy'
  Envi_Write_Envi_File, finePre, Out_Name = OutName, r_fid=fid_temp,$
     ns = ns, nl = nl, nb = nb, MAP_INFO=map_info

  envi_file_mng, id = fid1, /remove
  envi_file_mng, id = fid2, /remove
  envi_file_mng, id = fid3, /remove
  ;envi_file_mng, id = fid_temp, /remove

  print, 'time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'
end