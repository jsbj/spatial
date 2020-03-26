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

  newlat@units = "degrees_north"
  newlon@units = "degrees_east"
  
  ntime = dimsizes(time)
  nlat = dimsizes(newlat)
  nlon = dimsizes(newlon)

  LW_clear_monthly = f->LW_clear_monthly
  LW_cloud_monthly = f->LW_cloud_monthly
  SW_cloud_monthly = f->SW_cloud_monthly
  SW_clear_monthly = f->SW_clear_monthly

  newLW_cloud_monthly = linint2(lon,lat,LW_cloud_monthly,True,newlon,newlat,0)
  newLW_cloud_monthly!0   = "time"
  newLW_cloud_monthly!1   = "lat"
  newLW_cloud_monthly!2   = "lon"
  newLW_cloud_monthly&time = time
  newLW_cloud_monthly&lat = newlat
  newLW_cloud_monthly&lon = newlon

  newLW_clear_monthly = linint2(lon,lat,LW_clear_monthly,True,newlon,newlat,0)
  newLW_clear_monthly!0   = "time"
  newLW_clear_monthly!1   = "lat"
  newLW_clear_monthly!2   = "lon"
  newLW_clear_monthly&time = time
  newLW_clear_monthly&lat = newlat
  newLW_clear_monthly&lon = newlon

  newSW_clear_monthly = linint2(lon,lat,SW_clear_monthly,True,newlon,newlat,0)
  newSW_clear_monthly!0   = "time"
  newSW_clear_monthly!1   = "lat"
  newSW_clear_monthly!2   = "lon"
  newSW_clear_monthly&time = time
  newSW_clear_monthly&lat = newlat
  newSW_clear_monthly&lon = newlon

  newSW_cloud_monthly = linint2(lon,lat,SW_cloud_monthly,True,newlon,newlat,0)
  newSW_cloud_monthly!0   = "time"
  newSW_cloud_monthly!1   = "lat"
  newSW_cloud_monthly!2   = "lon"
  newSW_cloud_monthly&time = time
  newSW_cloud_monthly&lat = newlat
  newSW_cloud_monthly&lon = newlon

  LW_clear_annual = f->LW_clear_annual
  LW_cloud_annual = f->LW_cloud_annual
  SW_cloud_annual = f->SW_cloud_annual
  SW_clear_annual = f->SW_clear_annual

  newLW_cloud_annual = linint2(lon,lat,LW_cloud_annual,True,newlon,newlat,0)
  newLW_cloud_annual!0   = "time"
  newLW_cloud_annual!1   = "lat"
  newLW_cloud_annual!2   = "lon"
  newLW_cloud_annual&time = time
  newLW_cloud_annual&lat = newlat
  newLW_cloud_annual&lon = newlon

  newLW_clear_annual = linint2(lon,lat,LW_clear_annual,True,newlon,newlat,0)
  newLW_clear_annual!0   = "time"
  newLW_clear_annual!1   = "lat"
  newLW_clear_annual!2   = "lon"
  newLW_clear_annual&time = time
  newLW_clear_annual&lat = newlat
  newLW_clear_annual&lon = newlon

  newSW_clear_annual = linint2(lon,lat,SW_clear_annual,True,newlon,newlat,0)
  newSW_clear_annual!0   = "time"
  newSW_clear_annual!1   = "lat"
  newSW_clear_annual!2   = "lon"
  newSW_clear_annual&time = time
  newSW_clear_annual&lat = newlat
  newSW_clear_annual&lon = newlon

  newSW_cloud_annual = linint2(lon,lat,SW_cloud_annual,True,newlon,newlat,0)
  newSW_cloud_annual!0   = "time"
  newSW_cloud_annual!1   = "lat"
  newSW_cloud_annual!2   = "lon"
  newSW_cloud_annual&time = time
  newSW_cloud_annual&lat = newlat
  newSW_cloud_annual&lon = newlon
  
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
  filevardef(new_file, "LW_clear_monthly"   ,typeof(LW_clear_monthly) ,(/"time","lat","lon"/))
  filevardef(new_file, "LW_cloud_monthly"   ,typeof(LW_cloud_monthly) ,(/"time","lat","lon"/))
  filevardef(new_file, "SW_cloud_monthly"   ,typeof(SW_cloud_monthly) ,(/"time","lat","lon"/))
  filevardef(new_file, "SW_clear_monthly"   ,typeof(SW_clear_monthly) ,(/"time","lat","lon"/))
  filevardef(new_file, "LW_clear_annual"   ,typeof(LW_clear_annual) ,(/"time","lat","lon"/))
  filevardef(new_file, "LW_cloud_annual"   ,typeof(LW_cloud_annual) ,(/"time","lat","lon"/))
  filevardef(new_file, "SW_cloud_annual"   ,typeof(SW_cloud_annual) ,(/"time","lat","lon"/))
  filevardef(new_file, "SW_clear_annual"   ,typeof(SW_clear_annual) ,(/"time","lat","lon"/))

  
  filevarattdef(new_file,"time",time)                    ; copy time attributes
  filevarattdef(new_file,"lat",lat)                     ; copy lat attributes
  filevarattdef(new_file,"lon",lon)                     ; copy lon attributes
  filevarattdef(new_file,"LW_clear_monthly",LW_clear_monthly)
  filevarattdef(new_file,"LW_cloud_monthly",LW_cloud_monthly)
  filevarattdef(new_file,"SW_cloud_monthly",SW_cloud_monthly)
  filevarattdef(new_file,"SW_clear_monthly",SW_clear_monthly)
  filevarattdef(new_file,"LW_clear_annual",LW_clear_annual)
  filevarattdef(new_file,"LW_cloud_annual",LW_cloud_annual)
  filevarattdef(new_file,"SW_cloud_annual",SW_cloud_annual)
  filevarattdef(new_file,"SW_clear_annual",SW_clear_annual)

  
  setfileoption(new_file,"DefineMode",False)
  
  new_file->time = (/time/)
  new_file->lat = (/newlat/)
  new_file->lon = (/newlon/)
  new_file->LW_clear_monthly = (/newLW_clear_monthly/)
  new_file->LW_cloud_monthly = (/newLW_cloud_monthly/)
  new_file->SW_cloud_monthly = (/newSW_cloud_monthly/)
  new_file->SW_clear_monthly = (/newSW_clear_monthly/)
  new_file->LW_clear_annual = (/newLW_clear_annual/)
  new_file->LW_cloud_annual = (/newLW_cloud_annual/)
  new_file->SW_cloud_annual = (/newSW_cloud_annual/)
  new_file->SW_clear_annual = (/newSW_clear_annual/)
end
