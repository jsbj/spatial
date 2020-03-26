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
  rlutcs = f->rlutcs
  rlut = f->rlut
  rsut = f->rsut
  rsutcs = f->rsutcs

  newlat@units = "degrees_north"
  newlon@units = "degrees_east"
  
  ntime = dimsizes(time)
  nlat = dimsizes(newlat)
  nlon = dimsizes(newlon)
  
  newtas = area_conserve_remap(lon,lat,tas,newlon,newlat,False)
  newrlut = area_conserve_remap(lon,lat,rlut,newlon,newlat,False)
  newrlutcs = area_conserve_remap(lon,lat,rlutcs,newlon,newlat,False)
  newrsutcs = area_conserve_remap(lon,lat,rsutcs,newlon,newlat,False)
  newrsut = area_conserve_remap(lon,lat,rsut,newlon,newlat,False)
  
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
  filevardef(new_file, "rlutcs"   ,typeof(rlutcs) ,(/"time","lat","lon"/))
  filevardef(new_file, "rlut"   ,typeof(rlut) ,(/"time","lat","lon"/))
  filevardef(new_file, "rsut"   ,typeof(rsut) ,(/"time","lat","lon"/))
  filevardef(new_file, "rsutcs"   ,typeof(rsutcs) ,(/"time","lat","lon"/))
  
  filevarattdef(new_file,"time",time)                    ; copy time attributes
  filevarattdef(new_file,"lat",lat)                     ; copy lat attributes
  filevarattdef(new_file,"lon",lon)                     ; copy lon attributes
  filevarattdef(new_file,"tas",tas)
  ; filevarattdef(new_file,"rlutcs",rlutcs)
  ; filevarattdef(new_file,"rlut",rlut)
  ; filevarattdef(new_file,"rsut",rsut)
  ; filevarattdef(new_file,"rsutcs",rsutcs)
  
  setfileoption(new_file,"DefineMode",False)
  
  new_file->time = (/time/)
  new_file->lat = (/newlat/)
  new_file->lon = (/newlon/)
  new_file->tas = (/newtas/)
  new_file->rlutcs = (/newrlutcs/)
  new_file->rlut = (/newrlut/)
  new_file->rsut = (/newrsut/)
  new_file->rsutcs = (/newrsutcs/)
end
