; fn

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  
begin
  
  f = addfile(in_fn,"r")
  time = f->time
  lat = f->lat
  lon = f->lon
  tas = f->tas
  LW_clear_monthly = f->LW_clear_monthly
  LW_cloud_monthly = f->LW_cloud_monthly
  SW_cloud_monthly = f->SW_cloud_monthly
  SW_clear_monthly = f->SW_clear_monthly

  newlat@units = "degrees_north"
  newlon@units = "degrees_east"
  
  ntime = dimsizes(time)
  nlat = dimsizes(newlat)
  nlon = dimsizes(newlon)
  
  newtas = area_conserve_remap(lon,lat,tas,newlon,newlat,False)
  newLW_cloud_monthly = area_conserve_remap(lon,lat,LW_cloud_monthly,newlon,newlat,False)
  newLW_clear_monthly = area_conserve_remap(lon,lat,LW_clear_monthly,newlon,newlat,False)
  newSW_clear_monthly = area_conserve_remap(lon,lat,SW_clear_monthly,newlon,newlat,False)
  newSW_cloud_monthly = area_conserve_remap(lon,lat,SW_cloud_monthly,newlon,newlat,False)
  
  dir_name = "./"
  new_file = addfile(out_fn, "c")
  
  setfileoption(new_file,"DefineMode",True)
  
  fAtt = True
  fAtt@title = "This file is the coolest"
  fAtt@source_file = in_fn
  fAtt@conventions = "None"
  fAtt@creation_date = systemfunc("date")
  fileattdef(new_file,fAtt)

  dim_names = (/"time", "lat", "lon"/)
  dim_sizes = (/ -1, nlat, nlon /)
  dim_unlim = (/ True, False, False /)
  filedimdef(new_file,dim_names,dim_sizes,dim_unlim)

  filevardef(new_file, "time" ,typeof(time),(/"time"/))
  filevardef(new_file, "lat"  ,typeof(newlat),(/"lat"/))
  filevardef(new_file, "lon"  ,typeof(newlon),(/"lon"/))
  filevardef(new_file, "tas"  ,typeof(newtas),(/"time","lat","lon"/))
  filevardef(new_file, "LW_clear_monthly"   ,typeof(LW_clear_monthly) ,(/"time","lat","lon"/))
  filevardef(new_file, "LW_cloud_monthly"   ,typeof(LW_cloud_monthly) ,(/"time","lat","lon"/))
  filevardef(new_file, "SW_cloud_monthly"   ,typeof(SW_cloud_monthly) ,(/"time","lat","lon"/))
  filevardef(new_file, "SW_clear_monthly"   ,typeof(SW_clear_monthly) ,(/"time","lat","lon"/))
  
  filevarattdef(new_file,"time",time)                    ; copy time attributes
  filevarattdef(new_file,"lat",lat)                     ; copy lat attributes
  filevarattdef(new_file,"lon",lon)                     ; copy lon attributes
  filevarattdef(new_file,"tas",tas)
  filevarattdef(new_file,"LW_clear_monthly",LW_clear_monthly)
  filevarattdef(new_file,"LW_cloud_monthly",LW_cloud_monthly)
  filevarattdef(new_file,"SW_cloud_monthly",SW_cloud_monthly)
  filevarattdef(new_file,"SW_clear_monthly",SW_clear_monthly)
  
  setfileoption(new_file,"DefineMode",False)
  
  new_file->time = (/time/)
  new_file->lat = (/newlat/)
  new_file->lon = (/newlon/)
  new_file->tas = (/newtas/)
  new_file->LW_clear_monthly = (/newLW_clear_monthly/)
  new_file->LW_cloud_monthly = (/newLW_cloud_monthly/)
  new_file->SW_cloud_monthly = (/newSW_cloud_monthly/)
  new_file->SW_clear_monthly = (/newSW_clear_monthly/)
end
