

function estimate_and_analyze(model_name,run_type,feedback_type,grid="orig",years_used=:all,adjacence=:all,tropics_width=:all,polar_mask=:none,feedback_model_name=nothing)
  estimate(model_name,run_type,feedback_type,grid,years_used,adjacence,tropics_width,polar_mask,feedback_model_name)
  feedbacks(model_name,run_type,:global_global_time_scales,grid,false,feedback_type,years_used,adjacence,tropics_width,polar_mask,feedback_model_name)
  # if years_used == :all
    if !(feedback_type in [:global_local_multiple,:global_local_multiple_regridded]) && finite_difference(model_name,run_type,grid,feedback_type,false,years_used,adjacence,tropics_width,polar_mask)
      finite_difference_difference(model_name,run_type,grid,feedback_type,false,years_used,adjacence,tropics_width,polar_mask)
    end
  # end
end

# function regrid_and_analyze(model_name,run_type,feedback_type,grid,lat_lon,years_used=:all,feedback_model_name=nothing)
#   regrid_estimate(model_name,run_type,feedback_type,grid,lat_lon,years_used,feedback_model_name)
#   feedbacks(model_name,run_type,:global_global_time_scales,grid,true,feedback_type,years_used,:all,:all,:none,feedback_model_name)
#   # if years_used == :all
#     if finite_difference(model_name,run_type,grid,feedback_type,true,years_used)
#       finite_difference_difference(model_name,run_type,grid,feedback_type,true,years_used)
#     end
#   # end
# end

# function reset_time_scales()
  # rm ../netcdfs/*/finite_difference/*
  # rm ../netcdfs/*/finite_difference_difference/*
#   rm ../jlds/orig/*time_scales*
#   rm ../jlds/orig/estimated/*time_scales*
#   rm ../jlds/12x24_equaldist/*time_scales*
#   rm ../jlds/12x24_equaldist/estimated/*time_scales*
#   rm ../jlds/scores/*/*
# end

# function reset_scores()
  # rm ../jlds/scores/*/*
# end

function create_model_dict()
  model_dict = Dict{Any,Any}()
  
  for control_fn in glob("netcdfs/orig/raw/tas/*control*")
    split_tas_fn = split(split(control_fn,"/")[end],"_")
    model_name = split_tas_fn[3]
    years = Int(parse(Float64,split(split_tas_fn[5],".")[1])) # Int(split(split_tas_fn[5],".")[1])
    model_dict[model_name] = Dict{Any,Any}("control" => years)
    
    abrupt_fns = glob("netcdfs/orig/raw/tas/*$(model_name)*abrupt*")
    if length(abrupt_fns) > 0
      model_dict[model_name]["abrupts"] = Dict{Any,Any}()
      for abrupt_fn in abrupt_fns
        model_dict[model_name]["abrupts"][split(abrupt_fn,"_")[4]] = Int(parse(Float64,split(split(abrupt_fn,"_")[5],".")[1]))
      end
    end
    
    ramp_fns = glob("netcdfs/orig/raw/tas/*$(model_name)*1pct*")
    if length(ramp_fns) > 0
      model_dict[model_name]["ramps"] = Dict{Any,Any}()
      for ramp_fn in ramp_fns
        model_dict[model_name]["ramps"][split(ramp_fn,"_")[4]] = Int(parse(Float64,split(split(ramp_fn,"_")[5],".")[1]))
      end
    end
    
    rcp_fns = glob("netcdfs/orig/raw/tas/*$(model_name)*rcp*")
    if length(rcp_fns) > 0
      model_dict[model_name]["rcps"] = Dict{Any,Any}()
      for rcp_fn in rcp_fns
        model_dict[model_name]["rcps"][split(rcp_fn,"_")[4]] = Int(parse(Float64,split(split(rcp_fn,"_")[5],".")[1]))
      end
    end
  end
  
  model_dict
end

function standardize(model_name,run_type)
  fn = to_standardized_fn(model_name,run_type)
  
  if !isfile(fn)
    println("making standardized $(fn)")
    println(to_raw_fn("tas",model_name,run_type))
    dset = to_dset(to_raw_fn("tas",model_name,run_type))
    if "latitude_1" in dset.keys()
      dset = dset.rename(Dict("latitude_1" => "lat","longitude_1" => "lon","longitude_1_bnds" => "lon_bnds","latitude_1_bnds" => "lat_bnds"))
    end
  
    for other_var in ["rlut","rlutcs","rsut","rsutcs"] #,"ts"]
      println(to_raw_fn(other_var,model_name,run_type))
      other_dset = to_dset(to_raw_fn(other_var,model_name,run_type))
      dset = dset.assign(; [(Symbol(other_var),(("time","lat","lon"),getproperty(other_dset,Symbol(other_var)).values))]...)
    end
    
    if dset.lat.values[1] > 0
      println("reversing latitudes")
      dset.lat.values = dset.lat.values[end:-1:1]
      for var_name in ["tas","rlut","rlutcs","rsut","rsutcs"]
        getproperty(dset,Symbol(var_name)).values = getproperty(dset,Symbol(var_name)).values[:,end:-1:1,:]
      end
    end
    
    if dset.time.values[end] == 0
      dset.time.values = collect(1:length(dset.time.values))
    end
    
    dset = dset.drop(setdiff(dset.keys(),["lat","lon","time","tas","rsut","rsutcs","rlut","rlutcs"]))
  
    println(fn)
    dset.to_netcdf(fn)
  else
    println("skipping standardized $(fn)")
  end
end

function crop(model_name,run_type,years)
  dset = to_dset(to_standardized_fn(model_name,run_type))
  dset = dset.isel(time=(years*12):(length(dset.time.values)-1))
  dset.to_netcdf(to_standardized_fn(model_name,run_type)*".tmp",engine="scipy")
  mv(to_standardized_fn(model_name,run_type)*".tmp",to_standardized_fn(model_name,run_type),remove_destination=true)
end

function crop_subset(model_name,run_type,indices)
  dset = to_dset(to_standardized_fn(model_name,run_type))
  dset = dset.isel(time=indices) # (years*12):(length(dset["time"].values)-1))
  dset.to_netcdf(to_standardized_fn(model_name,run_type)*".tmp",engine="scipy")
  mv(to_standardized_fn(model_name,run_type)*".tmp",to_standardized_fn(model_name,run_type),remove_destination=true)
end

function anomalize(model_name,run_type,grid="orig")
  if !isfile(to_anomalized_fn(model_name,run_type,grid))
    println("making anomalies for $(model_name)")
    if model_name in ["CCSM3","CNRMCM61","HadCM3L"]
      if !isfile(to_anomalized_fn(model_name,run_type,grid))
        run(`cp $(to_standardized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))`)
      end
      dset = to_dset(to_anomalized_fn(model_name,run_type,grid))
      
      for var_name in ["tas","rsutcs","rsut","rlutcs","rlut"]
        println("doing variable $(var_name)")
        control_climatology = to_climatology(to_dset(to_standardized_fn(model_name,"control",grid))[var_name].values)

        if dset.data_vars.keys() == var_names
          println("skipping " * to_anomalized_fn(model_name,run_type,grid))
          continue
        else
          data = getproperty(dset,Symbol(var_name)).values
          dset = dset.drop([var_name])
          to_anomaly!(data,control_climatology)
          if var_name == "tas"
            println("T_surf")
            dset = dset.assign(T_surf = (("time","lat","lon"),data))
          elseif var_name[end-1:end] == "cs"
            println(Symbol(uppercase(var_name[2]) * "W_clear"))
            dset = dset.assign(; [(Symbol(uppercase(var_name[2]) * "W_clear"), (("time","lat","lon"),-data))]...)
          else
            println(Symbol(uppercase(var_name[2]) * "W_cloud"))
            dset = dset.assign(; [(Symbol(uppercase(var_name[2]) * "W_cloud"), (("time","lat","lon"),-data-getproperty(dset,Symbol(uppercase(var_name[2]) * "W_clear")).values))]...)
          end
        end
        println("about to remove $(to_anomalized_fn(model_name,run_type,grid))")
        run(`rm $(to_anomalized_fn(model_name,run_type,grid))`)
        println("removed $(to_anomalized_fn(model_name,run_type,grid))")
        dset.to_netcdf(to_anomalized_fn(model_name,run_type,grid)) #,engine="scipy")
        println("saved new $(to_anomalized_fn(model_name,run_type,grid))")
      end
    else
      println("cdo ymonsub $(to_standardized_fn(model_name,run_type,grid)) -ymonmean $(to_standardized_fn(model_name,"control",grid)) $(to_anomalized_fn(model_name,run_type,grid))")
      run(`cdo ymonsub $(to_standardized_fn(model_name,run_type,grid)) -ymonmean $(to_standardized_fn(model_name,"control",grid)) $(to_anomalized_fn(model_name,run_type,grid))`)
      println("ncap2 -O -s 'LW_clear=-rlutcs' $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))")
      run(`ncap2 -O -s 'LW_clear=-rlutcs' $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))`)
      println("ncap2 -O -s 'LW_cloud=-(rlut-rlutcs)' $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))")
      run(`ncap2 -O -s 'LW_cloud=-(rlut-rlutcs)' $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))`)
      println("ncks -O -x -v rlut,rlutcs $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))")
      run(`ncks -O -x -v rlut,rlutcs $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))`)
      println("ncap2 -O -s 'SW_clear=-rsutcs' $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))")
      run(`ncap2 -O -s 'SW_clear=-rsutcs' $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))`)
      println("ncap2 -O -s 'SW_cloud=-(rsut-rsutcs)' $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))")
      run(`ncap2 -O -s 'SW_cloud=-(rsut-rsutcs)' $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))`)
      println("ncks -O -x -v rsut,rsutcs $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))")
      run(`ncks -O -x -v rsut,rsutcs $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))`)
      println("ncrename -O -v tas,T_surf $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))")
      run(`ncrename -O -v tas,T_surf $(to_anomalized_fn(model_name,run_type,grid)) $(to_anomalized_fn(model_name,run_type,grid))`)
    end
    
    globalize!(to_anomalized_fn(model_name,run_type,grid))
  else
    println("skipping $(to_anomalized_fn(model_name,run_type,grid))")
  end
  nothing
end


function globalize!(fn::String)
  println("globalizing $(fn)")
  dset = globalize(to_dset(fn))
  # if any([occursin("global",var_name) for var_name in vars(dset)])
  #   println("skipping globalizing $(fn)")
  # else
  # mv(to_anomalized_fn(model_name,run_type),to_anomalized_fn(model_name,run_type)*".tmp")
  # dset[:to_netcdf](to_anomalized_fn(model_name,run_type))
  # rm(to_anomalized_fn(model_name,run_type)*".tmp")
  mv(fn,fn*".tmp")  
  dset.to_netcdf(fn)
  rm(fn*".tmp")
  # end
end

function variance!(fn::String)
  println("variance $(fn)")
  dset = globalize(to_dset(fn))
  # if any([occursin("global",var_name) for var_name in vars(dset)])
  #   println("skipping globalizing $(fn)")
  # else
  # mv(to_anomalized_fn(model_name,run_type),to_anomalized_fn(model_name,run_type)*".tmp")
  # dset[:to_netcdf](to_anomalized_fn(model_name,run_type))
  # rm(to_anomalized_fn(model_name,run_type)*".tmp")
  mv(fn,fn*".tmp")  
  dset.to_netcdf(fn)
  rm(fn*".tmp")
  # end
end

function globalize(dset)
  for var_name in dset.data_vars.keys() # ["T_surf","LW_clear","LW_cloud","SW_clear","SW_cloud"]
    if !occursin("global",var_name)
      global_var_name = var_name * "_global"
      if !(global_var_name in vars(dset))
        println(global_var_name)
        dset = dset.assign(; [(Symbol(global_var_name), (("time"),to_spatial_average(dset,var_name)))]...)
      end
    end
  end
  dset
end

function globalize(dset)
  for var_name in dset.data_vars.keys() # ["T_surf","LW_clear","LW_cloud","SW_clear","SW_cloud"]
    if !occursin("global",var_name)
      global_var_name = var_name * "_global"
      if !(global_var_name in vars(dset))
        println(global_var_name)
        dset = dset.assign(; [(Symbol(global_var_name), (("time"),to_spatial_average(dset,var_name)))]...)
      end
    end
  end
  dset
end

function variance(dset)
  for var_name in dset.data_vars.keys() # ["T_surf","LW_clear","LW_cloud","SW_clear","SW_cloud"]
    if !occursin("global",var_name) && !occursin("variance",var_name)
      global_var_name = var_name * "_variance"
      if !(global_var_name in vars(dset))
        println(global_var_name)
        dset = dset.assign(; [(Symbol(global_var_name), (("time"),to_spatial_average(dset,var_name)))]...)
      end
    end
  end
  dset
end

function special_globalize(dset)
  for var_name in dset.data_vars.keys() # ["T_surf","LW_clear","LW_cloud","SW_clear","SW_cloud"]
    if !occursin("global",var_name)
      special_global_var_name = var_name * "_special_global"
      if !(special_global_var_name in vars(dset))
        dset = dset.assign(; [(Symbol(special_global_var_name), (("time"),to_special_spatial_average(dset,var_name)))]...)
      end
    end
  end
  dset
end


function feedbacks(model_name,run_type,feedback_type,grid="orig",regridded_estimate=false,estimate_type=:none,years_used=:all,adjacence=:all,tropics_width=:all,polar_mask=:none,feedback_model_name=nothing)
  fn = estimate_type == :none ? to_anomalized_fn(model_name,run_type,grid) : (regridded_estimate ? to_regridded_estimated_fn(model_name,run_type,estimate_type,grid,years_used=years_used,feedback_model_name=feedback_model_name) : to_estimated_fn(model_name,run_type,estimate_type,grid,years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=polar_mask,feedback_model_name=feedback_model_name))
  println(fn)
  dset = to_dset(fn)
  feedback_fn = to_feedback_fn(model_name,run_type,feedback_type,grid,estimate_type=estimate_type,years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=polar_mask,feedback_model_name=feedback_model_name)
  # println(feedback_fn)
  if !isfile(feedback_fn)
    println("making feedbacks for $(feedback_fn)")
    if estimate_type != :none
      # println(join(vcat(split(fn,"/")[(1:2).+length(split(fn_prefix,"/")).-1],"anomalized",join(split(split(fn,"/")[4+length(split(fn_prefix,"/")).-1],"_")[1:2],"_")*".nc"),"/"))
      T_surf = to_array(to_dset(fn_prefix * join(vcat(split(fn,"/")[(1:2).+length(split(fn_prefix,"/")).-1],"anomalized",join(split(split(fn,"/")[4+length(split(fn_prefix,"/")).-1],"_")[1:2],"_")*".nc"),"/")),"T_surf",:global)
    else
      if feedback_type in [:local_local,:global_local,:global_local_multiple,:global_local_multiple_first_half,:global_local_multiple_second_half,:global_local_multiple_regridded,:full,:full_first_half,:full_second_half]
        T_surf,nlat,nlon = to_vector(dset.T_surf.values)
      else
        T_surf = to_array(dset,"T_surf",:global)
      end

      if (years_used != :all) || (feedback_type in [:full_first_half,:full_second_half,:global_local_multiple_first_half,:global_local_multiple_second_half,])
        years = Int(size(T_surf)[1]/12)
        if years_used != :all
          if to_full_time_scale_boolean(T_surf,"$(years-years_used+1) to end","monthly")
            T_surf = to_time_scale(T_surf,"$(years-years_used+1) to end","monthly")
          else
            println("years used longer than run_length")
            return false
          end
        elseif feedback_type in [:full_first_half,:global_local_multiple_first_half]
          T_surf = to_time_scale(T_surf,"1 to $(Int(floor(years/2)))","monthly")
        elseif feedback_type in [:full_second_half,:global_local_multiple_second_half]
          T_surf = to_time_scale(T_surf,"$(Int(ceil((years+1)/2))) to end","monthly")
        end
      end
    end

    if (feedback_type == :global_global_time_scales) && (length(T_surf) <= 600)
      println("$(fn) not long enough for our analysis")
      return nothing
    end

    if feedback_type in [:global_global,:global_global_time_scales,:full,:full_first_half,:full_second_half]
      feedback_dict = Dict()
      for freq in freqs
        for sky in skies
          if feedback_type == :global_global
            fluxes = to_array(dset,freq*sky,:global)
            if years_used != :all
              fluxes = to_time_scale(fluxes,"$(years-years_used+1) to end","monthly")
            end
            feedback_dict[freq*sky] = to_annual_and_monthly_feedback_dict(T_surf,fluxes)
          elseif feedback_type == :global_global_time_scales
            if estimate_type != :none
              T_surf_ann = to_annual_average(T_surf)
              feedback_dict[freq*sky] = OrderedDict()
              for time_scale in time_scales_dict[run_type]
                if to_full_time_scale_boolean(T_surf_ann,time_scale)
                  # println("T_surf_ann: $(length(to_time_scale(to_fancy_average(T_surf_ann),time_scale)))")
                  # println("flux ann: $(length(to_time_scale(to_fancy_average(to_array(dset,freq*sky*"_annual",:global)),time_scale)))")
                  # println("flux mon: $(length(to_time_scale(to_fancy_average(to_array(dset,freq*sky*"_monthly",:global)),time_scale)))")
                  println("third time")
                  println(dset.keys())
                  feedback_dict[freq*sky][time_scale] = Dict(
                    "annual" => to_feedback_dict(to_time_scale(to_fancy_average(T_surf_ann),time_scale),to_time_scale(to_fancy_average(to_array(dset,freq*sky*"_annual",:global)),time_scale)),
                    "monthly" => to_feedback_dict(to_time_scale(to_fancy_average(T_surf_ann),time_scale),to_time_scale(to_fancy_average(to_array(dset,freq*sky*"_monthly",:global)),time_scale)),
                    "month_agnostic" => to_feedback_dict(to_time_scale(to_fancy_average(T_surf_ann),time_scale),to_time_scale(to_fancy_average(to_array(dset,freq*sky*"_month_agnostic",:global)),time_scale)),
                    "seasonal" => to_feedback_dict(to_time_scale(to_fancy_average(T_surf_ann),time_scale),to_time_scale(to_fancy_average(to_array(dset,freq*sky*"_seasonal",:global)),time_scale))
                  )
                  
                  if any([split(var,"_")[end-1] == "special" for var in vars(dset)])
                    feedback_dict[freq*sky][time_scale]["special_annual"] = to_feedback_dict(to_time_scale(to_fancy_average(T_surf_ann),time_scale),to_time_scale(to_fancy_average(to_array(dset,freq*sky*"_annual",:special_global)),time_scale))
                    feedback_dict[freq*sky][time_scale]["special_monthly"] = to_feedback_dict(to_time_scale(to_fancy_average(T_surf_ann),time_scale),to_time_scale(to_fancy_average(to_array(dset,freq*sky*"_monthly",:special_global)),time_scale))
                    feedback_dict[freq*sky][time_scale]["special_month_agnostic"] = to_feedback_dict(to_time_scale(to_fancy_average(T_surf_ann),time_scale),to_time_scale(to_fancy_average(to_array(dset,freq*sky*"_month_agnostic",:special_global)),time_scale))
                    feedback_dict[freq*sky][time_scale]["special_seasonal"] = to_feedback_dict(to_time_scale(to_fancy_average(T_surf_ann),time_scale),to_time_scale(to_fancy_average(to_array(dset,freq*sky*"_seasonal",:special_global)),time_scale))
                  end
                end
              end
            else
              fluxes = to_array(dset,freq*sky,:global)
              special_fluxes = to_array(dset,freq*sky,:special_global)
              feedback_dict[freq*sky] = OrderedDict()
              for time_scale in time_scales_dict[run_type]
                if to_full_time_scale_boolean(T_surf,time_scale,"monthly")
                  feedback_dict[freq*sky][time_scale] = to_feedback_dict(to_time_scale(to_fancy_average(to_annual_average(T_surf)),time_scale),to_time_scale(to_fancy_average(to_annual_average(fluxes)),time_scale))
                  # to_annual_seasonal_and_monthly_feedback_dict(to_time_scale(T_surf,time_scale,"monthly"),to_time_scale(fluxes,time_scale,"monthly"))
                  feedback_dict[freq*sky][time_scale]["special_annual"] = to_feedback_dict(to_time_scale(to_fancy_average(to_annual_average(T_surf)),time_scale),to_time_scale(to_fancy_average(to_annual_average(special_fluxes)),time_scale))
                end
              end
            end
          elseif feedback_type in [:full,:full_first_half,:full_second_half]
            fluxes,nlat,nlon = to_vector(to_array(dset,freq*sky,:full))
            if years_used != :all
              fluxes = to_time_scale(fluxes,"$(years-years_used+1) to end","monthly")
            elseif feedback_type == :full_first_half
              fluxes = to_time_scale(fluxes,"1 to $(Int(floor(years/2)))","monthly")
            elseif feedback_type == :full_second_half
              fluxes = to_time_scale(fluxes,"$(Int(ceil((years+1)/2))) to end","monthly")
            end
            if adjacence == :all
              feedback_dict[freq*sky] = to_annual_seasonal_and_monthly_feedback_dict(T_surf,fluxes,true)
            else
              feedback_dict[freq*sky] = to_annual_seasonal_and_monthly_feedback_dict(T_surf,fluxes,nlat,nlon,adjacence,tropics_width,polar_mask) #; adjacence=adjacence)
            end
          end
        end
      end
      save(feedback_fn,"feedback_dict",feedback_dict)
    else
      feedback_dset = xr.Dataset(Dict(),coords=Dict("lat" => dset.lat.values, "lon" => dset.lon.values))
      feedback_dset.lat.attrs = Dict("units"=>"degrees_north")
      feedback_dset.lon.attrs = Dict("units"=>"degrees_east")
      if feedback_type in [:global_local_multiple,:global_local_multiple_first_half,:global_local_multiple_second_half]
        orig_dset = to_dset(to_anomalized_fn(model_name,run_type,grid))
      end
      for freq in freqs
        for sky in skies
          if feedback_type in [:global_local_multiple,:global_local_multiple_first_half,:global_local_multiple_second_half]
            fluxes = to_array(orig_dset,freq*sky,:global)
            special_fluxes = to_array(orig_dset,freq*sky,:special_global)
          elseif feedback_type in [:global_local,:global_local_multiple_regridded]
            fluxes = to_array(dset,freq*sky,:global)
            special_fluxes = to_array(dset,freq*sky,:special_global)
          else
            fluxes,nlat,nlon = to_vector(to_array(dset,freq*sky,:full))
          end

          if (years_used != :all)
            fluxes = to_time_scale(fluxes,"$(years-years_used+1) to end","monthly")
            special_fluxes = to_time_scale(special_fluxes,"$(years-years_used+1) to end","monthly")
          elseif feedback_type == :global_local_multiple_first_half
            fluxes = to_time_scale(fluxes,"1 to $(Int(floor(years/2)))","monthly")
            special_fluxes = to_time_scale(special_fluxes,"1 to $(Int(floor(years/2)))","monthly")
          elseif feedback_type == :global_local_multiple_second_half
            fluxes = to_time_scale(fluxes,"$(Int(ceil((years+1)/2))) to end","monthly")
            special_fluxes = to_time_scale(special_fluxes,"$(Int(ceil((years+1)/2))) to end","monthly")
          end

          feedback_dict = to_annual_seasonal_and_monthly_feedback_dict(T_surf,fluxes,feedback_type in [:global_local_multiple,:global_local_multiple_first_half,:global_local_multiple_second_half,:global_local_multiple_regridded])
          if !(feedback_type in [:local_global,:local_local])
            special_feedback_dict = to_annual_seasonal_and_monthly_feedback_dict(T_surf,special_fluxes,feedback_type in [:global_local_multiple,:global_local_multiple_first_half,:global_local_multiple_second_half,:global_local_multiple_regridded])
          end
          for var_name in ["feedback"] #,"feedback_unct","R2","est_s2"]
            if var_name in ["R2","est_s2"] 
              if feedback_type in [:global_local_multiple,:global_local_multiple_first_half,:global_local_multiple_second_half,:global_local_multiple_regridded]
                feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_annual_"*var_name),(feedback_dict["annual"][var_name][1,1]))]...)
                feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_month_agnostic_"*var_name),(feedback_dict["month_agnostic"][var_name][1,1]))]...)
                feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_special_annual_"*var_name),(special_feedback_dict["annual"][var_name][1,1]))]...)
                feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_special_month_agnostic_"*var_name),(special_feedback_dict["month_agnostic"][var_name][1,1]))]...)
              end
            else
              feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_annual_"*var_name),(("lat","lon"),reshape(feedback_dict["annual"][var_name],(nlat,nlon))))]...)
              feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_month_agnostic_"*var_name),(("lat","lon"),reshape(feedback_dict["month_agnostic"][var_name],(nlat,nlon))))]...)
              if feedback_type == :global_local_multiple_regridded
                feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_special_annual_"*var_name),(("lat","lon"),reshape(special_feedback_dict["annual"][var_name],(nlat,nlon))))]...)
                feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_special_month_agnostic_"*var_name),(("lat","lon"),reshape(special_feedback_dict["month_agnostic"][var_name],(nlat,nlon))))]...)
              end
            end
            for i in 1:12
              if var_name in ["R2","est_s2"]
                if feedback_type in [:global_local_multiple,:global_local_multiple_first_half,:global_local_multiple_second_half,:global_local_multiple_regridded]
                  feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_monthly_$(lpad(i,2,"0"))_"*var_name),(feedback_dict["monthly"][i][var_name][1,1]))]...)
                  feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_special_monthly_$(lpad(i,2,"0"))_"*var_name),(special_feedback_dict["monthly"][i][var_name][1,1]))]...)
                end
              else
                feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_monthly_$(lpad(i,2,"0"))_"*var_name),(("lat","lon"),reshape(feedback_dict["monthly"][i][var_name],(nlat,nlon))))]...)
                if feedback_type == :global_local_multiple_regridded
                  feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_special_monthly_$(lpad(i,2,"0"))_"*var_name),(("lat","lon"),reshape(special_feedback_dict["monthly"][i][var_name],(nlat,nlon))))]...)
                end
              end
            end
            for i in 1:4 # MAM -> 1, JJA -> 2, etc.
              if var_name in ["R2","est_s2"]
                if feedback_type in [:global_local_multiple,:global_local_multiple_first_half,:global_local_multiple_second_half,:global_local_multiple_regridded]
                  feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_seasonal_$(lpad(i,2,"0"))_"*var_name),(feedback_dict["seasonal"][i][var_name][1,1]))]...)
                  feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_special_seasonal_$(lpad(i,2,"0"))_"*var_name),(special_feedback_dict["seasonal"][i][var_name][1,1]))]...)
                end
              else
                feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_seasonal_$(lpad(i,2,"0"))_"*var_name),(("lat","lon"),reshape(feedback_dict["seasonal"][i][var_name],(nlat,nlon))))]...)
                if feedback_type == :global_local_multiple_regridded
                  feedback_dset = feedback_dset.assign(; [(Symbol(freq*sky*"_special_seasonal_$(lpad(i,2,"0"))_"*var_name),(("lat","lon"),reshape(special_feedback_dict["seasonal"][i][var_name],(nlat,nlon))))]...)
                end
              end
            end
          end
        end
      end
      feedback_dset.to_netcdf(feedback_fn)
    end
  else
    println("skipping feedback $(feedback_fn)")
  end
end

# function nonlocal_and_local_from_full(model_name,grid;years_used=:all,adjacence=:all,feedback_model_name = nothing)
# end
#
# model_name = "CESM104"
# run_type = "abrupt4x"
# feedback_type = :full_tropics
# grid = "12x24_equaldist"
# years_used=:all
# adjacence=:all
# tropics_width=:all
# polar_mask=:none
# feedback_model_name = nothing
function estimate(model_name,run_type,feedback_type,grid="orig",years_used=:all,adjacence=:all,tropics_width=:all,polar_mask=:none,feedback_model_name = nothing;tropics=false)
  fn = to_estimated_fn(model_name,run_type,feedback_type,grid;years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=polar_mask,feedback_model_name=(feedback_model_name == model_name ? nothing : feedback_model_name),tropics=tropics)
  if !isfile(fn)
    println("making estimates for $(fn)")

    anomaly_dset = to_dset(to_anomalized_fn(model_name,run_type,grid))
    estimate_dset = xr.Dataset(Dict(),coords=Dict("time" => to_annual_average(anomaly_dset.time.values), "lat" => anomaly_dset.lat.values, "lon" => anomaly_dset.lon.values))
  
    if string(feedback_type)[1:4] == "full"
      if feedback_type in [:full_first_half,:full_second_half]
        feedback_dict = to_feedback_dict(to_feedback_fn((feedback_model_name == nothing ? model_name : feedback_model_name),"control",feedback_type,grid; years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=polar_mask))
      else
        feedback_dict = to_feedback_dict(to_feedback_fn((feedback_model_name == nothing ? model_name : feedback_model_name),"control",:full,grid; years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=polar_mask))
      end
    else
      feedback_dset = to_dset(to_feedback_fn((feedback_model_name == nothing ? model_name : feedback_model_name),"control",feedback_type,grid; years_used=years_used))
    end
  
    full_T_surf = (feedback_type in [:local_local,:global_local_multiple,:global_local_multiple_regridded]) || (string(feedback_type)[1:4] == "full")
    
    T_surf = getproperty(anomaly_dset,Symbol((full_T_surf ? "T_surf" : "T_surf_global"))).values
    if full_T_surf && tropics
      new_T_surf = zeros(size(T_surf))
      new_T_surf[:,6:7,:] = T_surf[:,6:7,:]
      T_surf = new_T_surf
    end
    T_surf_ann = to_annual_average(T_surf)
    ntime = div(length(anomaly_dset.time.values),12)
    nlat = length(anomaly_dset.lat.values)
    nlon = length(anomaly_dset.lon.values)
    # freq = "LW"
    # sky = "_clear"
    for freq in ["LW","SW"]
      for sky in ["_clear","_cloud"]
        if string(feedback_type)[1:4] == "full"
          T_surf_ann_vector,nlat,nlon = to_vector(T_surf_ann)
          feedback_matrix = to_filtered_feedback_matrix(feedback_dict[freq*sky]["annual"]["feedback"],feedback_type)
          ann_est = to_map(T_surf_ann_vector*feedback_matrix',nlat,nlon)

          T_surf_mon_agnostic_vector,nlat,nlon = to_vector(T_surf)
          mon_agnostic_feedback_matrix = to_filtered_feedback_matrix(feedback_dict[freq*sky]["month_agnostic"]["feedback"],feedback_type)
          mon_agnostic_est = to_annual_average(to_map(T_surf_mon_agnostic_vector*mon_agnostic_feedback_matrix',nlat,nlon))
        else
          ann_est = (full_T_surf ? T_surf_ann : reshape(T_surf_ann,(ntime,1,1))) .* reshape(getproperty(feedback_dset,Symbol("$(freq)$(sky)_annual_feedback")).values,(1,nlat,nlon))
          mon_agnostic_est = to_annual_average((full_T_surf ? T_surf : reshape(T_surf,(12ntime,1,1))) .* reshape(getproperty(feedback_dset,Symbol("$(freq)$(sky)_month_agnostic_feedback")).values,(1,nlat,nlon)))
        end

        estimate_dset = estimate_dset.assign(; [(Symbol(freq*sky*"_annual"),(("time","lat","lon"),ann_est))]...)
        estimate_dset = estimate_dset.assign(; [(Symbol(freq*sky*"_month_agnostic"),(("time","lat","lon"),mon_agnostic_est))]...)

        mon_est = zeros(size(ann_est))
        seas_est = zeros(size(ann_est))
        for i in 1:12
          if string(feedback_type)[1:4] == "full"
            T_surf_mon_vector,nlat,nlon = to_vector(T_surf[i:12:end,:,:])
            monthly_feedback_matrix = to_filtered_feedback_matrix(feedback_dict[freq*sky]["monthly"][i]["feedback"],feedback_type)
            seasonal_feedback_matrix = to_filtered_feedback_matrix(feedback_dict[freq*sky]["seasonal"][to_season_i(i)]["feedback"],feedback_type)
            mon_est += to_map(T_surf_mon_vector*monthly_feedback_matrix',nlat,nlon)
            seas_est += to_map(T_surf_mon_vector*seasonal_feedback_matrix',nlat,nlon)
          else
            mon_est += (full_T_surf ? T_surf[i:12:end,:,:] : reshape(T_surf[i:12:end],(ntime,1,1))) .* reshape(getproperty(feedback_dset,Symbol(freq*sky*"_monthly_$(lpad(i,2,"0"))_feedback")).values,(1,nlat,nlon))
            seas_est += (full_T_surf ? T_surf[i:12:end,:,:] : reshape(T_surf[i:12:end],(ntime,1,1))) .* reshape(getproperty(feedback_dset,Symbol(freq*sky*"_seasonal_$(lpad(to_season_i(i),2,"0"))_feedback")).values,(1,nlat,nlon))
          end
        end
        estimate_dset = estimate_dset.assign(; [(Symbol(freq*sky*"_monthly"),(("time","lat","lon"),mon_est/12))]...)
        estimate_dset = estimate_dset.assign(; [(Symbol(freq*sky*"_seasonal"),(("time","lat","lon"),seas_est/12))]...)
        if feedback_type in [:global_local_multiple,:global_local_multiple_regridded]
          estimate_dset = estimate_dset.assign(; [(Symbol(freq*sky*"_annual_global"),(("time"),dropdims(sum(getproperty(estimate_dset,Symbol(freq*sky*"_annual")).values,dims=(2,3)),dims=(2,3))))]...)
          estimate_dset = estimate_dset.assign(; [(Symbol(freq*sky*"_month_agnostic_global"),(("time"),dropdims(sum(getproperty(estimate_dset,Symbol(freq*sky*"_month_agnostic")).values,dims=(2,3)),dims=(2,3))))]...)
          estimate_dset = estimate_dset.assign(; [(Symbol(freq*sky*"_monthly_global"),(("time"),dropdims(sum(getproperty(estimate_dset,Symbol(freq*sky*"_monthly")).values,dims=(2,3)),dims=(2,3))))]...)
          estimate_dset = estimate_dset.assign(; [(Symbol(freq*sky*"_seasonal_global"),(("time"),dropdims(sum(getproperty(estimate_dset,Symbol(freq*sky*"_seasonal")).values,dims=(2,3)),dims=(2,3))))]...)
        end
      end
    end 
    special_globalize(globalize(estimate_dset)).to_netcdf(fn)
  else
    println("skipping estimates for $(fn)")
  end
end

function regrid(model_name,run_type,grid,lat_lon)
  if !isdir(fn_prefix * "netcdfs/$(grid)")
    mkdir(fn_prefix * "netcdfs/$(grid)")
    mkdir(fn_prefix * "netcdfs/$(grid)/standardized")
    mkdir(fn_prefix * "netcdfs/$(grid)/anomalized")
    mkdir(fn_prefix * "netcdfs/$(grid)/feedbacks")
    mkdir(fn_prefix * "netcdfs/$(grid)/estimated")
    mkdir(fn_prefix * "netcdfs/$(grid)/regridded_estimated")
    mkdir(fn_prefix * "netcdfs/$(grid)/finite_difference")
    mkdir(fn_prefix * "netcdfs/$(grid)/finite_difference_difference")
    mkdir(fn_prefix * "netcdfs/$(grid)/normalized_finite_difference")
    mkdir(fn_prefix * "netcdfs/$(grid)/contributions")
    mkdir(fn_prefix * "jlds/$(grid)")
    mkdir(fn_prefix * "jlds/$(grid)/estimated")
  end

  in_fn = to_standardized_fn(model_name,run_type)
  out_fn = to_standardized_fn(model_name,run_type,grid)
  
  regrid_cmd(in_fn,out_fn,lat_lon)
end

function regrid_cmd(in_fn,out_fn,lat_lon,is_est=false)
  if !isfile(out_fn)  
    println("regridding $(in_fn) to $(out_fn)")
    # cd("scripts")
    if is_est == true
      # Shell.run("ncl 'in_fn=\"$("../"*in_fn)\"' 'out_fn=\"$("../"*out_fn)\"' 'newlat=(/$(join(lat_lon["lats"],","))/)' 'newlon=(/$(join(lat_lon["lons"],","))/)' regrid_est.ncl")
      # if is_est == :had
        println("here")
        Shell.run("ncl 'in_fn=\"$("../"*in_fn)\"' 'out_fn=\"$("../"*out_fn)\"' 'newlat=(/$(join(lat_lon["lats"],","))/)' 'newlon=(/$(join(lat_lon["lons"],","))/)' regrid_had.ncl")        
      # end
    else
      println("ncl 'in_fn=\"$("../"*in_fn)\"' 'out_fn=\"$("../"*out_fn)\"' 'newlat=(/$(join(lat_lon["lats"],","))/)' 'newlon=(/$(join(lat_lon["lons"],","))/)' regrid.ncl")
      # Shell.run("ncl 'in_fn=\"$("../"*in_fn)\"' 'out_fn=\"$("../"*out_fn)\"' 'newlat=(/$(join(lat_lon["lats"],","))/)' 'newlon=(/$(join(lat_lon["lons"],","))/)' regrid.ncl")
    end
    cd("..")
  else
    println("skipping $(out_fn)")
  end  
end

function regrid_estimate(model_name,run_type,feedback_type,grid,lat_lon,years_used=:all,feedback_model_name=nothing)
  fn = to_regridded_estimated_fn(model_name,run_type,feedback_type,grid;years_used=years_used,feedback_model_name=feedback_model_name)
  
  if !isfile(fn)
    println("making $(fn)")
    regrid_cmd(to_estimated_fn(model_name,run_type,feedback_type;years_used=years_used,feedback_model_name=feedback_model_name),fn,lat_lon,true)
    globalize!(fn)
  else
    println("skipping $(fn)")
  end
end

function finite_difference(model_name,run_type,grid="orig",estimate_type=:none,regridded_estimate=false,years_used=:all,adjacence=:all,tropics_width=:all,polar_mask=:none;tropics=false)
  fn = to_finite_difference_fn(model_name,run_type,grid,estimate_type,regridded_estimate;years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=polar_mask,tropics=tropics)
  if !isfile(fn)
    println("calculating finite difference for $(fn)")
    if estimate_type == :none
      starting_dset = to_dset(to_anomalized_fn(model_name,run_type,grid))
      factor = 12
    else
      if regridded_estimate
        starting_dset = to_dset(to_regridded_estimated_fn(model_name,run_type,estimate_type,grid;years_used=years_used))
      else
        starting_dset = to_dset(to_estimated_fn(model_name,run_type,estimate_type,grid;years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=polar_mask,tropics=tropics))
      end
      factor = 1
    end
    new_dset = xr.Dataset(Dict(),coords=Dict("lat" => starting_dset.lat.values, "lon" => starting_dset.lon.values))
    
    time_scales = time_scales_dict[run_type]
    
    if length(starting_dset.time) < factor * Int(parse(Float64,split(time_scales[end]," ")[1]))
      println("too short for finite diff")
      return false
    end
    
    nts = length(time_scales)
    # println("nts: ",nts)
    for var_name in vars(starting_dset)


      for i in 1:nts
        # println("evaluate ----")
        # println(size(starting_dset[var_name].values))
        # println(time_scales[i])
        # println(estimate_type == :none ? "monthly" : "annual")
        # println(to_full_time_scale_boolean(starting_dset[var_name].values,time_scales[i],estimate_type == :none ? "monthly" : "annual"))
        if to_full_time_scale_boolean(getproperty(starting_dset,Symbol(var_name)).values,time_scales[i],estimate_type == :none ? "monthly" : "annual")
          time_scale = to_time_scale(getproperty(starting_dset,Symbol(var_name)).values,time_scales[i],estimate_type == :none ? "monthly" : "annual")
          len = size(time_scale)[1]
          
          # println("iteration ",i)
          
          if nts == 2
            if i == 1
              range_1 = estimate_type == :none ? (1:60) : (1:5)
              range_2 = estimate_type == :none ? (61:228) : (6:19)
            elseif i == 2
              range_1 = estimate_type == :none ? (1:1800) : (1:150)
              range_2 = estimate_type == :none ? (1801:len) : (151:len)
            end
          elseif nts == 1
            range_1 = 1:Int(floor(len/2))
            range_2 = Int(ceil(len/2)):len
          end
          new_var_name = split(var_name,"_")[end] == "global" ? join(vcat(split(var_name,"_")[1:(end-1)],[i],["global"]),"_") : var_name * "_$i"
          start_part = dropdims(mean(time_scale[vcat([range_1],fill(Colon(),length(size(time_scale))-1))...],dims=1),dims=1)
          end_part = dropdims(mean(time_scale[vcat([range_2],fill(Colon(),length(size(time_scale))-1))...],dims=1),dims=1)
          new_dset = to_added_var(new_dset,new_var_name,("lat","lon")[1:(length(size(time_scale))-1)],end_part - start_part)
          
          if var_name == "T_surf"
            for m in 1:12
              time_scale_by_month = time_scale[m:12:end,:,:]
              len = size(time_scale_by_month)[1]
              if nts == 2
                if i == 1
                  range_1 = 1:5
                  range_2 = 6:19
                elseif i == 2
                  range_1 = 1:150
                  range_2 = 151:len
                end
              else
                range_1 = 1:Int(floor(len/2))
                range_2 = Int(ceil(len/2)):len
              end
              
              new_var_name = split(var_name,"_")[end] == "global" ? join(vcat(split(var_name,"_")[1:(end-1)],[i],["global"]),"_") : var_name * "_$(i)_m$(lpad(m,2,"0"))"
              start_part = dropdims(mean(time_scale_by_month[vcat([range_1],fill(Colon(),length(size(time_scale_by_month))-1))...],dims=1),dims=1)
              end_part = dropdims(mean(time_scale_by_month[vcat([range_2],fill(Colon(),length(size(time_scale_by_month))-1))...],dims=1),dims=1)
              new_dset = to_added_var(new_dset,new_var_name,("lat","lon")[1:(length(size(time_scale))-1)],end_part - start_part)
            end
          end
        end
      end
    end
    new_dset.to_netcdf(fn)
  else
    println("skipping finite difference for $(fn)")
  end
  return true
end

function normalized_finite_difference(model_name,run_type,grid="orig",estimate_type=:none,regridded_estimate=false,years_used=:all,adjacence=:all,tropics_width=:all,polar_mask=:none;tropics=false)
  fn = to_normalized_finite_difference_fn(model_name,run_type,grid,estimate_type,regridded_estimate;years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=polar_mask,tropics=tropics)
  
  if !isfile(fn)
    println("making $(fn)")
    true_dset = to_dset(to_finite_difference_fn(model_name,run_type,grid))
    new_dset = xr.Dataset(Dict(),coords=Dict("lat" => true_dset.lat.values, "lon" => true_dset.lon.values))
    fd_dset = to_dset(to_finite_difference_fn(model_name,run_type,grid,estimate_type,regridded_estimate;years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=polar_mask))
    
    if estimate_type == :none
      for flux in ["LW","SW"], sky in ["clear","cloud"], t in 1:maximum([(split(x,"_")[3] == "special") ? 0 : Int(parse(Float64,split(x,"_")[3])) for x in vars(true_dset)]), spatial_avg in ["","_global"]
        est_var_name = join([flux,sky,t],"_")*spatial_avg
        new_dset = new_dset.assign(; [(Symbol(est_var_name),(("lat","lon")[1:(isempty(spatial_avg) ? 2 : 0)],getproperty(fd_dset,Symbol(est_var_name)).values./getproperty(true_dset,Symbol(:T_surf_,t,:_global)).values))]...)
      end
    else
      for flux in ["LW","SW"], sky in ["clear","cloud"], t in 1:maximum([(split(x,"_")[3] == "special") ? 0 : Int(parse(Float64,split(x,"_")[3])) for x in vars(true_dset)]), average_type in ["monthly","annual","month_agnostic","seasonal"], spatial_avg in ["","_global"]
        est_var_name = join([flux,sky,average_type,t],"_")*spatial_avg
        new_dset = new_dset.assign(; [(Symbol(est_var_name),(("lat","lon")[1:(isempty(spatial_avg) ? 2 : 0)],getproperty(fd_dset,Symbol(est_var_name)).values./getproperty(true_dset,Symbol(:T_surf_,t,:_global)).values))]...)
      end
    end
    new_dset.to_netcdf(fn)
  else
    println("skipping $(fn)")
  end
end

function finite_difference_difference(model_name,run_type,grid="orig",estimate_type=:none,regridded_estimate=false,years_used=:all,adjacence=:all,tropics_width=:all,polar_mask=:none)
  true_dset = to_dset(to_finite_difference_fn(model_name,run_type,grid))
  est_dset = to_dset(to_finite_difference_fn(model_name,run_type,grid,estimate_type,regridded_estimate;years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=polar_mask))
  new_dset = xr.Dataset(Dict(),coords=Dict("lat" => true_dset.lat.values, "lon" => true_dset.lon.values))
  
  fn = to_finite_difference_difference_fn(model_name,run_type,grid,estimate_type,regridded_estimate;years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=polar_mask)
  if !isfile(fn)
    println("making $(fn)")
    # make difference
    for flux in ["LW","SW"], sky in ["clear","cloud"], t in 1:maximum([(split(x,"_")[3] == "special") ? 0 : Int(parse(Float64,split(x,"_")[3])) for x in vars(true_dset)]), average_type in ["monthly","annual","month_agnostic","seasonal"], spatial_avg in ["","_global"]
      est_var_name = join([flux,sky,average_type,t],"_")*spatial_avg
      var_name = join([flux,sky,t],"_")*spatial_avg
      # println(est_dset)
      # println(true_dset)
      new_dset = new_dset.assign(; [(Symbol(est_var_name),(("lat","lon")[1:(isempty(spatial_avg) ? 2 : 0)],getproperty(est_dset,Symbol(est_var_name)).values - getproperty(true_dset,Symbol(var_name)).values))]...)
      if spatial_avg == ""      
        new_dset = new_dset.assign(; [(Symbol(est_var_name*"_error"),((),sqrt(to_spatial_average(getproperty(new_dset,Symbol(est_var_name)).values.^2,new_dset.lat.values))))]...)
      end
      # new_dset = new_dset[:assign](; [(Symbol(est_var_name),(("lat","lon")[1:(isempty(spatial_avg) ? 2 : 0)],abs../(new_dset[est_var_name].values)))]...)
    end
    
    for t in 1:maximum([(split(x,"_")[3] == "special") ? 0 : Int(parse(Float64,split(x,"_")[3])) for x in vars(true_dset)]), average_type in ["monthly","annual","month_agnostic","seasonal"], spatial_avg in ["","_global"]
      est_var_name = join(["",average_type,t],"_")*spatial_avg
      answer = zeros(size(getproperty(new_dset,Symbol("LW_clear"*est_var_name)).values))
      for flux in ["LW_clear","LW_cloud","SW_clear","SW_cloud"]
        answer .+= getproperty(new_dset,Symbol(flux*est_var_name)).values
      end
      new_dset = new_dset.assign(; [(Symbol("N"*est_var_name),(("lat","lon")[1:(isempty(spatial_avg) ? 2 : 0)],answer))]...)
      if spatial_avg == ""      
        new_dset = new_dset.assign(; [(Symbol("N"*est_var_name*"_error"),((),sqrt(to_spatial_average(getproperty(new_dset,Symbol("N"*est_var_name)).values.^2,new_dset.lat.values))))]...)
      end
    end

    new_dset.to_netcdf(fn)
  else
    println("skipping $(fn)")
  end
end

# function score_dict(score_i,avg_type,grid,run_fn,full_type=:none,years_used=:all,adjacence=:all,tropics_width=:all)
#   if tropics_width == :all
#     adjacence = :all
#   end
#   fn = to_score_fn(score_i,avg_type,grid,run_fn,full_type,years_used,adjacence,tropics_width)
#   if !isfile(fn)
#     println("making $(fn)")
#
#     score_dict = Dict()
#     for var_name in ["LW_clear","LW_cloud","SW_clear","SW_cloud","N"]
#       score_dict[var_name] = Dict(
#         "global" => getfield(Main, Symbol("to_score_$(score_i)"))(:local_global,run_fn,avg_type,var_name,score_i == 1 ? "orig" : grid,years_used),
#         "local" => getfield(Main, Symbol("to_score_$(score_i)"))(:local_local,run_fn,avg_type,var_name,score_i == 1 ? "orig" : grid,years_used)
#       )
#
#       if (years_used == :all) || (years_used > to_min_years(adjacence,tropics_width,length(grid_dict[grid]["lats"]),length(grid_dict[grid]["lons"])))
#         score_dict[var_name]["full"] = getfield(Main, Symbol("to_score_$(score_i)"))(full_type == :none ? (score_i == 2 ? :full : :global_local_multiple) : full_type,run_fn,avg_type,var_name,grid,years_used,adjacence,tropics_width)
#       end
#     end
#
#     save(fn,"score_dict",score_dict)
#   else
#     println("skipping $(fn)")
#   end
# end


function score_dict(score_i,avg_type,run_fn,full_type=:none,years_used=:all,adjacence=:all,tropics_width=:all)
  if tropics_width == :all
    adjacence = :all
  end
  fn = to_score_fn(score_i,avg_type,run_fn,full_type,years_used,adjacence,tropics_width)
  if !isfile(fn)
    println("making $(fn)")
    
    score_dict = Dict()
    for var_name in ["LW_clear","LW_cloud","SW_clear","SW_cloud","N"]
      score_dict[var_name] = Dict(
        "global" => getfield(Main, Symbol("to_score_$(score_i)"))(:local_global,run_fn,avg_type,var_name,years_used),
        "local" => getfield(Main, Symbol("to_score_$(score_i)"))(:local_local,run_fn,avg_type,var_name,years_used)
      )
        
      score_dict[var_name]["full"] = getfield(Main, Symbol("to_score_$(score_i)"))(full_type == :none ? (score_i == 2 ? :full : :global_local_multiple) : full_type,run_fn,avg_type,var_name,years_used,adjacence,tropics_width)
    end
    
    save(fn,"score_dict",score_dict)
  else
    println("skipping $(fn)")
  end
end

function local_nonlocal(model_name,grid,lat_lon,half=:none)
  fn = to_feedback_fn(model_name,"control",Symbol("full_local_nonlocal$(half == :none ? "" : "_$(half)")"),grid)  

  if !isfile(fn)
    println("making $(fn)")
    dset = xr.Dataset(Dict(),coords=Dict("lat" => lat_lon["lats"], "lon" => lat_lon["lons"]))
    dset.lat[:attrs] = Dict("units"=>"degrees_north")
    dset.lon[:attrs] = Dict("units"=>"degrees_east")
  
    nlat = length(dset.lat)
    nlon = length(dset.lon)

    for (var_name,var_dict) in to_feedback_dict(to_feedback_fn(model_name,"control",Symbol("full$(half == :none ? "" : "_$(half)")"),grid))
      local_mat = Diagonal(diag(var_dict["annual"]["feedback"]))
      dset = to_added_var(dset,var_name*"_annual_feedback_local",("lat","lon"),reshape([to_spatial_average(reshape(local_mat[:,i],nlat,nlon),lat_lon["lats"]) for i in 1:(nlat * nlon)],nlat,nlon))
      nonlocal_mat = var_dict["annual"]["feedback"]-Diagonal(diag(var_dict["annual"]["feedback"]))
      dset = to_added_var(dset,var_name*"_annual_feedback_nonlocal",("lat","lon"),reshape([to_spatial_average(reshape(nonlocal_mat[:,i],nlat,nlon),lat_lon["lats"]) for i in 1:(nlat * nlon)],nlat,nlon))

      local_mat = Diagonal(diag(var_dict["month_agnostic"]["feedback"]))
      dset = to_added_var(dset,var_name*"_month_agnostic_feedback_local",("lat","lon"),reshape([to_spatial_average(reshape(local_mat[:,i],nlat,nlon),lat_lon["lats"]) for i in 1:(nlat * nlon)],nlat,nlon))
      nonlocal_mat = var_dict["month_agnostic"]["feedback"]-Diagonal(diag(var_dict["month_agnostic"]["feedback"]))
      dset = to_added_var(dset,var_name*"_month_agnostic_feedback_nonlocal",("lat","lon"),reshape([to_spatial_average(reshape(nonlocal_mat[:,i],nlat,nlon),lat_lon["lats"]) for i in 1:(nlat * nlon)],nlat,nlon))      

      for mon_i in 1:12
        local_mat = Diagonal(diag(var_dict["monthly"][mon_i]["feedback"]))
        dset = to_added_var(dset,var_name*"_monthly_$(lpad(mon_i,2,"0"))_feedback_local",("lat","lon"),reshape([to_spatial_average(reshape(local_mat[:,i],nlat,nlon),lat_lon["lats"]) for i in 1:(nlat * nlon)],nlat,nlon))
        nonlocal_mat = var_dict["monthly"][mon_i]["feedback"]-Diagonal(diag(var_dict["monthly"][mon_i]["feedback"]))
        dset = to_added_var(dset,var_name*"_monthly_$(lpad(mon_i,2,"0"))_feedback_nonlocal",("lat","lon"),reshape([to_spatial_average(reshape(nonlocal_mat[:,i],nlat,nlon),lat_lon["lats"]) for i in 1:(nlat * nlon)],nlat,nlon))
      end
      
      for seas_i in 1:4
        local_mat = Diagonal(diag(var_dict["seasonal"][seas_i]["feedback"]))
        dset = to_added_var(dset,var_name*"_seasonal_$(lpad(seas_i,2,"0"))_feedback_local",("lat","lon"),reshape([to_spatial_average(reshape(local_mat[:,i],nlat,nlon),lat_lon["lats"]) for i in 1:(nlat * nlon)],nlat,nlon))
        nonlocal_mat = var_dict["seasonal"][seas_i]["feedback"]-Diagonal(diag(var_dict["seasonal"][seas_i]["feedback"]))
        dset = to_added_var(dset,var_name*"_seasonal_$(lpad(seas_i,2,"0"))_feedback_nonlocal",("lat","lon"),reshape([to_spatial_average(reshape(nonlocal_mat[:,i],nlat,nlon),lat_lon["lats"]) for i in 1:(nlat * nlon)],nlat,nlon))
      end
    end
    dset.to_netcdf(fn)
  else
    println("skipping $(fn)")
  end
end

# contribution of each location to W/m2 per degree of warming after - contribution of each location to W/m2 per degree of warming before
# ((local warming * global local map) / global warming) 2 - ((local warming * global local map) / global warming) 1
function contribution(model_name,run_type,grid,lat_lon)
  fn = to_contribution_fn(model_name,run_type,grid)  

  if !isfile(fn)
    println("making contribution $(fn)")
    dset = xr.Dataset(Dict(),coords=Dict("lat" => lat_lon["lats"], "lon" => lat_lon["lons"]))
    dset.lat[:attrs] = Dict("units"=>"degrees_north")
    dset.lon[:attrs] = Dict("units"=>"degrees_east")
    nlat = length(dset.lat)
    nlon = length(dset.lon)

    f_dset = to_dset(to_feedback_fn(model_name,"control",:global_local_multiple_regridded,grid))
    fd_dset = to_dset(to_finite_difference_fn(model_name,run_type,grid))

    for var_name in var_names
      # dset = to_added_var(dset,var_name*"_annual",("lat","lon"),f_dset[var_name*"_annual_feedback"].values .* ΔΔT)

      before_part = (fd_dset.T_surf_1.values / fd_dset.T_surf_1_global.values[1])
      after_part = (fd_dset.T_surf_2.values / fd_dset.T_surf_2_global.values[1])
      ann = getproperty(f_dset,Symbol(var_name*"_annual_feedback")).values .* (after_part - before_part)
      special_ann = getproperty(f_dset,Symbol(var_name*"_special_annual_feedback")).values .* (after_part - before_part)

      before_part = mean([getproperty(fd_dset,Symbol("T_surf_1_m$(lpad(m,2,"0"))")).values .* getproperty(f_dset,Symbol(var_name * "_monthly_" * lpad(m,2,"0") * "_feedback")).values for m in 1:12]) / fd_dset.T_surf_1_global.values[1]
      after_part = mean([getproperty(fd_dset,Symbol("T_surf_2_m$(lpad(m,2,"0"))")).values .* getproperty(f_dset,Symbol(var_name * "_monthly_" * lpad(m,2,"0") * "_feedback")).values for m in 1:12]) / fd_dset.T_surf_2_global.values[1]
      mon = after_part - before_part        
      
      special_before_part = mean([getproperty(fd_dset,Symbol("T_surf_1_m$(lpad(m,2,"0"))")).values .* getproperty(f_dset,Symbol(var_name * "_special_monthly_" * lpad(m,2,"0") * "_feedback")).values for m in 1:12]) / fd_dset.T_surf_1_global.values[1]
      special_after_part = mean([getproperty(fd_dset,Symbol("T_surf_2_m$(lpad(m,2,"0"))")).values .* getproperty(f_dset,Symbol(var_name * "_special_monthly_" * lpad(m,2,"0") * "_feedback")).values for m in 1:12]) / fd_dset.T_surf_2_global.values[1]
      special_mon = special_after_part - special_before_part        
      
      before_part = mean([getproperty(fd_dset,Symbol("T_surf_1_m$(lpad(m,2,"0"))")).values .* getproperty(f_dset,Symbol(var_name * "_seasonal_" * lpad(m,2,"0") * "_feedback")).values for m in 1:4]) / fd_dset.T_surf_1_global.values[1]
      after_part = mean([getproperty(fd_dset,Symbol("T_surf_2_m$(lpad(m,2,"0"))")).values .* getproperty(f_dset,Symbol(var_name * "_seasonal_" * lpad(m,2,"0") * "_feedback")).values for m in 1:4]) / fd_dset.T_surf_2_global.values[1]
      seas = after_part - before_part        
      
      special_before_part = mean([getproperty(fd_dset,Symbol("T_surf_1_m$(lpad(m,2,"0"))")).values .* getproperty(f_dset,Symbol(var_name * "_special_seasonal_" * lpad(m,2,"0") * "_feedback")).values for m in 1:4]) / fd_dset.T_surf_1_global.values[1]
      special_after_part = mean([getproperty(fd_dset,Symbol("T_surf_2_m$(lpad(m,2,"0"))")).values .* getproperty(f_dset,Symbol(var_name * "_special_seasonal_" * lpad(m,2,"0") * "_feedback")).values for m in 1:4]) / fd_dset.T_surf_2_global.values[1]
      special_seas = special_after_part - special_before_part        

      dset = to_added_var(dset,var_name*"_annual",("lat","lon"),ann)
      dset = to_added_var(dset,var_name*"_monthly",("lat","lon"),mon)
      dset = to_added_var(dset,var_name*"_seasonal",("lat","lon"),seas)
      dset = to_added_var(dset,var_name*"_special_annual",("lat","lon"),special_ann)
      dset = to_added_var(dset,var_name*"_special_monthly",("lat","lon"),special_mon)
      dset = to_added_var(dset,var_name*"_special_seasonal",("lat","lon"),special_seas)
    end

    dset.to_netcdf(fn)
  else
    println("skipping $(fn)")
  end
end


# function add_feedback_dicts(feedback_dict_1,feedback_dict_2,T_surf,flux_1,flux_2)
#   flux_est1 = reshape(T_surf * feedback_dict_1["feedback"],length(T_surf))
#   flux_est2 = reshape(T_surf * feedback_dict_2["feedback"],length(T_surf))
#   f = sum((flux_1+flux_2).^2)
#
#   a = sum(flux_1.^2)
#   b = sum(flux_2.^2)
#   c = sum(flux_est1.^2)
#   d = sum(flux_est2.^2)
#   g = sum((flux_est1+flux_est2).^2)
#   R_fact = sum((flux_est1+flux_est2).^2)*sum(flux_1.^2)*sum(flux_2.^2) / (sum((flux_1+flux_2).^2)*(sum(flux_est1.^2)*sum(flux_2.^2)+sum(flux_est2.^2)*sum(flux_1.^2)))
#   s_fact = sum((flux_1+flux_2).^2-(flux_est1+flux_est2).^2)/sum(flux_1.^2+flux_2.^2-(flux_est1.^2+flux_est2.^2))
#
#   R_fact = g*a*b/(f*(c*b+a*d))
#   s_fact = (f-g)/(a+b-(c+d))
#   Dict(
#     "feedback" => feedback_dict_1["feedback"] + feedback_dict_2["feedback"],
#     "R2" => (feedback_dict_1["R2"] + feedback_dict_2["R2"]) * R_fact,
#     "est_s2" => s_fact*(feedback_dict_1["est_s2"]+feedback_dict_2["est_s2"]),
#     "feedback_est" => sqrt.(s_fact*(feedback_dict_1["feedback_unct"].^2+feedback_dict_2["feedback_unct"].^2))
#   )
# end

# model_name = "HadCM3L"
# run_type = "control"
# feedback_type = :full
# grid = "12x24_equaldist"
# regridded_estimate=false
# estimate_type=:none
# years_used=:all
# adjacence=:all
# feedback_model_name=nothing
# freq = "SW"
# sky = "_cloud"

# function mean_state(model_name,run_type,grid="orig",estimate_type=:none,regridded_estimate=false)
#   fn = to_mean_state_fn(model_name,run_type,grid,estimate_type,regridded_estimate)
#   if !isfile(fn)
#     println("calculating mean state for $(fn)")
#     if estimate_type == :none
#       starting_dset = to_dset(to_anomalized_fn(model_name,run_type,grid))
#     else
#       if regridded_estimate
#         starting_dset = to_dset(to_regridded_estimated_fn(model_name,run_type,estimate_type,grid))
#       else
#         starting_dset = to_dset(to_estimated_fn(model_name,run_type,estimate_type,grid))
#       end
#     end
#     new_dset = xr.Dataset(Dict(),coords=Dict("lat" => starting_dset["lat"].values, "lon" => starting_dset["lon"].values))
#     annual_time_scales = ["1 to 20","21 to 150","151 to end"]
#     monthly_time_scales = ["1 to 240","241 to 1800","1801 to end"]
#
#     for var_name in vars(starting_dset)
#       for i in 1:3
#         time_scale = to_time_scale(starting_dset[var_name].values,(estimate == :none ? monthly_time_scales[i] : annual_time_scales[i]))
#         len = size(time_scale)[1]
#         half = Int(floor(len/2))
#         new_var_name = split(var_name,"_")[end] == "global" ? join(vcat(split(var_name,"_")[1:(end-1)],[i],["global"]),"_"): var_name * "_$i"
#         new_dset = new_dset[:assign](; [(Symbol(new_var_name),(("lat","lon")[1:(length(size(time_scale))-1)],squeeze(mean(time_scale[vcat([(len-half+1):len],fill(Colon(),length(size(time_scale))-1))...],1),1) - squeeze(mean(time_scale[vcat([1:half],fill(Colon(),length(size(time_scale))-1))...],1),1)))]...)
#       end
#     end
#     new_dset[:to_netcdf](fn)
#   else
#     println("skipping finite difference for $(fn)")
#   end
# end


function area_weighted_map(grid)
  ds = xr.Dataset(Dict(),coords=Dict("lat" => grid_dict[grid]["lats"], "lon" => grid_dict[grid]["lons"]))
  ds["lat"][:attrs] = Dict("units"=>"degrees_north")
  ds["lon"][:attrs] = Dict("units"=>"degrees_east")
  nlon = length(grid_dict[grid]["lons"])
  nlat = length(grid_dict[grid]["lats"])
  lat_borders = -90:(180/nlat):90
  ds = ds[:assign](area_weights = (("lat","lon"),[(sin(lat_borders[lat_i+1]*π/180) - sin(lat_borders[lat_i]*π/180))/(2*nlon) for lat_i in 1:nlat, lon_i in 1:nlon]))
  ds[:to_netcdf]("netcdfs/$(grid)/area_weights.nc")
end

