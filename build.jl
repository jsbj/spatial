# future things to add
# - do ts instead of tas
include("load.jl")

# for (run_types,model_names) in model_name_dict
#   for model_name in model_names
#     standardize(model_name,"control")
#     for run_type in forced_runs(model_name,run_types)
#       standardize(model_name,run_type)
#     end
#   end
# end

# crop("CCSM3","control",530)
# cropped one of the GFDL models
# crop_subset("MIROC","control",vcat(collect(1:(517*12)),collect((1+518*12):(681*12)))-1)
#
# reset_time_scales()
function build()
  for model_name in all_model_names # ["GISSE2R","GFDLCM3"]# all_model_names
    # for run_type in ["control","abrupt4x"] # all_runs(model_name)
    #   standardize(model_name,run_type)
    #   anomalize(model_name,run_type)
    #   for (grid,lat_lon) in grid_dict
    #     regrid(model_name,run_type,grid,lat_lon)
    #     anomalize(model_name,run_type,grid)
    #   end
    # end

    # # # control
    # for feedback_type in [:global_global,:local_global,:local_local]
    #   feedbacks(model_name,"control",feedback_type,"orig",false,:none,:all)
    # end

    for (grid,lat_lon) in grid_dict
      for feedback_type in [:local_global,:local_local] # ,:global_global,:global_local_multiple] #:global_global_regridded, :global_global_regridded,
        println(model_name * " " * "$(feedback_type)")
        feedbacks(model_name,"control",feedback_type,grid)
      end

      for years_used in years_useds
        for feedback_type in [:global_local_multiple_regridded,:global_local_multiple_first_half,:global_local_multiple_second_half]
          if (years_used == :all) || (years_used > to_min_years(:all,:all,length(lat_lon["lats"]),length(lat_lon["lons"])))
            feedbacks(model_name,"control",feedback_type,grid,false,:none,years_used)
          end
        end

        for reduced in reduceds
          if (years_used == :all) || years_used > to_min_years(reduced[1],reduced[2],length(lat_lon["lats"]),length(lat_lon["lons"])) && years_used < length(to_dset(to_anomalized_fn(model_name,"control"))["time"])/12
            feedback_fn = to_feedback_fn(model_name,"control",:full,grid,years_used=years_used,adjacence=reduced[1],tropics_width=reduced[2],polar_mask=reduced[3],feedback_model_name=model_name)

            # midway([:feedbacks,model_name,"control",:full,grid,false,:none,years_used,reduced[1],reduced[2],reduced[3],model_name])
            for feedback_type in [:full,:full_first_half,:full_second_half]
              feedbacks(model_name,"control",feedback_type,grid,false,:none,years_used,reduced[1],reduced[2],reduced[3],model_name)
            end
          end
        end
      end
      for feedback_type in [:local_global,:local_local,:full]
        estimate(model_name,"control",feedback_type,grid)
        finite_difference(model_name,"control",grid,feedback_type)
      end
    end

    for (grid,lat_lon) in grid_dict
      local_nonlocal(model_name,grid,lat_lon)
    end

    for run_type in ["abrupt4x"] #forced_runs(model_name) # "abrupt4x",
      # # # time_scaled
      # special_globalize!(to_anomalized_fn(model_name,run_type))
      # feedbacks(model_name,run_type,:global_global_time_scales)
      # finite_difference(model_name,run_type)
      # # # end

      # for feedback_type in [:local_global,:local_local] # ,:global_local]
      #   estimate_and_analyze(model_name,run_type,feedback_type,"orig",:all)
      #   # estimate(model_name,run_type,feedback_type,"orig",years_used)
      #   feedbacks(model_name,run_type,:global_global_time_scales,"orig",false,feedback_type,:all)
      # end

      for (grid,lat_lon) in grid_dict
        globalize!(to_anomalized_fn(model_name,run_type,grid))
        special_globalize!(to_anomalized_fn(model_name,run_type,grid))
        feedbacks(model_name,run_type,:global_global_time_scales,grid)
        finite_difference(model_name,run_type,grid)


        estimate_and_analyze(model_name,run_type,:local_local,grid,:all)
        estimate_and_analyze(model_name,run_type,:local_global,grid,:all)

        for feedback_type in [:global_local_multiple_regridded] # :global_local_multiple] #,:global_local_multiple_regridded]
          estimate(model_name,run_type,feedback_type,grid,:all)
          feedbacks(model_name,run_type,:global_global_time_scales,grid,false,feedback_type,:all)
        end

        for years_used in years_useds
          # for feedback_type in [:local_global,:local_local] # :global_local] # ,:global_local_multiple,:global_local_multiple_regridded]
          #   regrid_and_analyze(model_name,run_type,feedback_type,grid,lat_lon,years_used)
          # end



          # for feedback_type in [:local_global_regridded,:local_local_regridded,:global_local_regridded]
          #   estimate_and_analyze(model_name,run_type,feedback_type,years_used)
          # end

          for reduced in reduceds
            # for feedback_model_name in all_model_names
            if (years_used == :all) || years_used > to_min_years(reduced[1],reduced[2],length(lat_lon["lats"]),length(lat_lon["lons"])) && years_used < length(to_dset(to_anomalized_fn(model_name,"control"))["time"])/12
              for feedback in [:full,:full_local,:full_nonlocal] # ,:full_tropics,:full_extratropics,:full_tropics_local,:full_extratropics_local,:full_tropics_nonlocal,:full_extratropics_nonlocal]
                estimate_and_analyze(model_name,run_type,feedback,grid,years_used,reduced[1],reduced[2],reduced[3],model_name)
                # efn = to_estimated_fn(model_name,run_type,:full,grid,years_used=years_used)
                # println(efn)
                # globalize!(efn)
                # special_globalize!(efn)
              end
              for feedback in [:full_first_half,:full_second_half]
                estimate_and_analyze(model_name,run_type,feedback,grid,years_used,reduced[1],reduced[2],reduced[3],model_name)
              end
            end
            # end
          end
        end

        contribution(model_name,run_type,grid,lat_lon)
      end
    end

    # 247: tropical ascent
    # 9: extratropics
    # 32: subsidence
    # to_individual_netcdf(model_name,"12x24_equaldist",[247,9,32])
    to_individual_netcdf(model_name,"12x24_equaldist",[258,9,32]) # what we use now
  end
end



grid = collect(keys(grid_dict))[1]
# []
[globalize!(to_anomalized_fn(model_name,"control",grid)) for model_name in all_model_names]
[special_globalize!(to_anomalized_fn(model_name,"control",grid)) for model_name in all_model_names]
# [globalize!(to_anomalized_fn(model_name,"abrupt4x",grid)) for model_name in all_model_names]
# [special_globalize!(to_anomalized_fn(model_name,"abrupt4x",grid)) for model_name in all_model_names]
build()








# [special_globalize!(to_estimated_fn(model_name,"abrupt4x",:local_local,grid)) for model_name in all_model_names]
# [special_globalize!(to_estimated_fn(model_name,"abrupt4x",:local_global,grid)) for model_name in all_model_names]
# [special_globalize!(to_estimated_fn(model_name,"abrupt4x",:full_local,grid)) for model_name in all_model_names]
# [special_globalize!(to_estimated_fn(model_name,"abrupt4x",:full_nonlocal,grid)) for model_name in all_model_names]
[[local_nonlocal(model_name,grid,grid_dict[grid],half) for model_name in all_model_names] for half in [:none,:first_half, :second_half,:none]]
# [estimate(model_name,"abrupt4x",:full,"12x24_equaldist",tropics=true) for model_name in all_model_names]
# [finite_difference(model_name,"abrupt4x","12x24_equaldist",:full,false,tropics=true) for model_name in all_model_names]
# [globalize!(to_estimated_fn(model_name,"abrupt4x",:full,"12x24_equaldist",tropics=true)) for model_name in all_model_names]
# [special_globalize!(to_estimated_fn(model_name,"abrupt4x",:full,"12x24_equaldist",tropics=true)) for model_name in all_model_names]
# [contribution(model_name,"abrupt4x","12x24_equaldist",grid_dict["12x24_equaldist"]) for model_name in all_model_names]
[normalized_finite_difference(model_name,"abrupt4x","12x24_equaldist") for model_name in all_model_names]
[normalized_finite_difference(model_name,"abrupt4x","12x24_equaldist",:full,false) for model_name in all_model_names]
[normalized_finite_difference(model_name,"abrupt4x","12x24_equaldist",:local_global,false) for model_name in all_model_names]
[normalized_finite_difference(model_name,"abrupt4x","12x24_equaldist",:local_local,false) for model_name in all_model_names]
[normalized_finite_difference(model_name,"abrupt4x","12x24_equaldist",:full_first_half,false) for model_name in all_model_names]
[normalized_finite_difference(model_name,"abrupt4x","12x24_equaldist",:full_second_half,false) for model_name in all_model_names]
#
# dset = to_dset(to_anomalized_fn("CESM104","control","12x24_equaldist"))
# variance_dset = xr.Dataset(Dict(),coords=Dict("lat" => dset.lat.values, "lon" => dset.lon.values))
# variance_dset.lat.attrs = Dict("units"=>"degrees_north")
# variance_dset.lon.attrs = Dict("units"=>"degrees_east")
# var_mean = mean([dropdims(var(to_dset(to_anomalized_fn(model_name,"control","12x24_equaldist")).T_surf.values,dims=1),dims=1) for model_name in all_model_names])
# variance_dset = variance_dset.assign(; [(:variance,(("lat","lon"),var_mean))]...)
# variance_dset.to_netcdf(to_variance_fn("12x24_equaldist"))