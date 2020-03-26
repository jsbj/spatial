
# file names
to_variance_fn(grid) = fn_prefix * "netcdfs/$(grid)/variance.nc"
to_raw_fn(var_name,model_name,run_type) = glob(fn_prefix * "netcdfs/orig/raw/$(var_name)/" * join([var_name,"mon",model_name,run_type,"*.nc"],"_"))[1]
to_standardized_fn(model_name,run_type,grid="orig") = fn_prefix * "netcdfs/$(grid)/standardized/" * join([model_name,run_type*".nc"],"_")
to_anomalized_fn(model_name,run_type,grid="orig") = fn_prefix * "netcdfs/$(grid)/anomalized/" * join([model_name,run_type*".nc"],"_")
to_full_suffix(years_used,adjacence,tropics_width,polar_mask,feedback_model_name) = (years_used==:all ? "" : "_$(years_used)used") * (adjacence==:all ? "" : "_$(adjacence)adjacent") * (tropics_width==:all ? "" : "_$(tropics_width)tropics") * (polar_mask==:none ? "" : "_$(polar_mask)polar") * (feedback_model_name == nothing ? "" : "_$(feedback_model_name)")
to_estimated_fn(model_name,run_type,feedback_type,grid="orig";years_used=:all,adjacence=:all,tropics_width=:all,polar_mask=:none,feedback_model_name=nothing,tropics=false) = fn_prefix * "netcdfs/$(grid)/estimated/$(model_name)_$(run_type)_$(feedback_type)" * to_full_suffix(years_used,adjacence,tropics_width,polar_mask,feedback_model_name == model_name ? nothing : feedback_model_name) * (tropics ? "_tropics" : "") * ".nc"
to_regridded_estimated_fn(model_name,run_type,feedback_type,grid="orig";years_used=:all,feedback_model_name=nothing) = fn_prefix * "netcdfs/$(grid)/regridded_estimated/$(model_name)_$(run_type)_$(feedback_type)" * to_full_suffix(years_used,:all,:all,:none,feedback_model_name == model_name ? nothing : feedback_model_name) * ".nc"
function to_feedback_fn(model_name,run_type,feedback_type,grid="orig";estimate_type=:none,years_used=:all,adjacence=:all,tropics_width=:all,polar_mask=:none,feedback_model_name=nothing)
  fn = if feedback_type in [:global_global,:global_global_time_scales,:full,:full_first_half,:full_second_half]
    if estimate_type == :none
      fn_prefix * "jlds/$(grid)/$(model_name)_$(run_type)_$(feedback_type)" * to_full_suffix(years_used,adjacence,tropics_width,polar_mask,feedback_model_name == model_name ? nothing : feedback_model_name) * ".jld"
    else
      fn_prefix * "jlds/$(grid)/estimated/$(model_name)_$(run_type)_$(estimate_type)_$(feedback_type)" * to_full_suffix(years_used,adjacence,tropics_width,polar_mask,feedback_model_name == model_name ? nothing : feedback_model_name) * ".jld"
    end
  else
    fn_prefix * "netcdfs/$(grid)/feedbacks/$(model_name)_$(run_type)_$(feedback_type)" * to_full_suffix(years_used,adjacence,tropics_width,polar_mask,feedback_model_name == model_name ? nothing : feedback_model_name) * ".nc"
  end
  # println(fn)
  fn
end
to_finite_difference_fn(model_name,run_type,grid="orig",estimate_type=:none,regridded_estimate=false;years_used=:all,adjacence=:all,tropics_width=:all,polar_mask=:none,tropics=false) = fn_prefix * "netcdfs/$(grid)/finite_difference/$(model_name)_$(run_type)$(estimate_type == :none ? "" : "_$(estimate_type)")$(regridded_estimate ? "_regridded" : "")" * to_full_suffix(years_used,adjacence,tropics_width,polar_mask,nothing) * (tropics ? "_tropics" : "") * ".nc"
to_normalized_finite_difference_fn(model_name,run_type,grid="orig",estimate_type=:none,regridded_estimate=false;years_used=:all,adjacence=:all,tropics_width=:all,polar_mask=:none,tropics=false) = fn_prefix * "netcdfs/$(grid)/normalized_finite_difference/$(model_name)_$(run_type)$(estimate_type == :none ? "" : "_$(estimate_type)")$(regridded_estimate ? "_regridded" : "")" * to_full_suffix(years_used,adjacence,tropics_width,polar_mask,nothing) * (tropics ? "_tropics" : "") * ".nc"
to_finite_difference_difference_fn(model_name,run_type,grid="orig",estimate_type=:none,regridded_estimate=false;years_used=:all,adjacence=:all,tropics_width=:all,polar_mask=:none) = fn_prefix * "netcdfs/$(grid)/finite_difference_difference/$(model_name)_$(run_type)$(estimate_type == :none ? "" : "_$(estimate_type)")$(regridded_estimate ? "_regridded" : "")" * to_full_suffix(years_used,adjacence,tropics_width,polar_mask,nothing) * ".nc"
to_contribution_fn(model_name,run_type,grid="orig") = fn_prefix * "netcdfs/$(grid)/contributions/$(model_name)_$(run_type).nc"
to_uncertainty_fn(model_name,time_scale,flux,sample,version,grid,area) = fn_prefix * "jlds/uncertainties/$(join([model_name,replace(time_scale," "=>"_"),flux,sample,version,grid,area],"_")).jld"

# to_mean_by(ar,i) = squeeze(mean(reshape(ar,(i,Int(length(ar)/i))),1)',2)
to_mean_by(ar,i) = dropdims(mean(reshape(ar,(i,Int(size(ar)[1]/i),size(ar)[2:length(size(ar))]...)),dims=1),dims=1)
to_fancy_average(list) = vcat(list[1:20],to_mean_by(list[21:50],5),to_mean_by(list[51:200],10),to_mean_by(list[201:Int(floor(length(list)/100)*100)],100))
to_fancy_years(model_name) = Int(floor(length(to_dset(to_anomalized_fn(model_name,"abrupt4x"))["time"])/1200)) + 38
to_years(model_name,run_name) = Int(floor(length(to_dset(to_anomalized_fn(model_name,run_name))["time"])/12))

to_dset(fn) = xr.open_dataset(fn,decode_times=false) #,engine="pynio")
# function to_dset(fn)
#   println(fn)
#   xr.open_dataset(fn,decode_times=false) #,engine="pynio")
# end
to_feedback_dict(fn) = load(fn)["feedback_dict"]
# to_score_fn(score_i,avg_type,grid,fn,full_type=:none,years_used=:all,adjacence=:all,tropics_width=:all) = "jlds/scores/$(score_i)/$(fn)_$(grid)_$(avg_type)$((full_type == :none) ? "" : "_$(full_type)")$(years_used == :all ? "" : "_$(years_used)")$(adjacence == :all ? "" : "_$(adjacence)")$(tropics_width == :all ? "" : "_$(tropics_width)").jld"
to_score_fn(score_i,avg_type,grid,fn,full_type=:none,years_used=:all,adjacence=:all,tropics_width=:all) = fn_prefix * "jlds/scores/$(score_i)/$(fn)_$(grid)_$(avg_type)$((full_type == :none) ? "" : "_$(full_type)")$(years_used == :all ? "" : "_$(years_used)")$(adjacence == :all ? "" : "_$(adjacence)")$(tropics_width == :all ? "" : "_$(tropics_width)").jld"
to_score_dict(fn) = load(fn)["score_dict"]

to_spatial_average(dset::PyCall.PyObject, field) = to_spatial_average(getproperty(dset,Symbol(field)).values,sort(dset.lat.values))

to_season_i(month_i) = Int(floor(mod(month_i - 3, 12) / 3))+1
# function to_spatial_average(vals,lats)
#   weights = sin.(vcat([-90],(lats[2:end]+lats[1:end-1])/2,[90])*π/180)
#   weights = (weights[2:end] - weights[1:end-1])/2
#   if length(size(vals)) == 2
#     squeeze(mean(vals,length(size(vals))),length(size(vals)))' * weights
#   else
#     squeeze(mean(vals,length(size(vals))),length(size(vals))) * weights
#   end
# end

function to_spatial_average(vals,lats)
  weights = cos.(lats*π/180)
  weights /= sum(weights)
  if length(size(vals)) == 2
    dropdims(mean(vals,dims=length(size(vals))),dims=length(size(vals)))' * weights
  else
    dropdims(mean(vals,dims=length(size(vals))),dims=length(size(vals))) * weights
  end
end

# just 30ºS and north
to_special_spatial_average(dset::PyCall.PyObject, field) = to_special_spatial_average(getproperty(dset,Symbol(field)).sel(lat = dset.lat.values[Int(1+length(dset.lat)/3):end]).values,sort(dset.lat.values[Int(1+length(dset.lat)/3):end]))

function to_special_spatial_average(vals,lats)
  weights = sin.(vcat([-30],(lats[2:end]+lats[1:end-1])/2,[90])*π/180)
  weights = (weights[2:end] - weights[1:end-1])/1.5
  if length(size(vals)) == 2
    dropdims(mean(vals,dims=length(size(vals))),dims=length(size(vals)))' * weights
  else
    dropdims(mean(vals,dims=length(size(vals))),dims=length(size(vals))) * weights
  end
end

to_time_scale(arr,phrase,avg_type="annual",return_type="annual") = arr[vcat([to_time_scale_range(arr,phrase,avg_type)],fill(Colon(),length(size(arr))-1))...]

function to_time_scale_range(arr,phrase,avg_type="annual")
  a = split(phrase," to ")
  true_lst = size(arr)[1]
  fst = avg_type != "annual" ? ((avg_type == "monthly" ? 12 : 4)*(Int(parse(Float64,a[1]))-1)+1) : Int(parse(Float64,a[1]))
  lst = (a[end] == "end" ? true_lst : (avg_type == "monthly" ? 12 : (avg_type == "seasonal" ? 4 : 1))*Int(parse(Float64,a[2])))
  fst:min(lst,true_lst)
end

function to_full_time_scale_boolean(arr,phrase,avg_type="annual")
  rng = to_time_scale_range(arr,phrase,avg_type)
  a = split(phrase," to ")
  if a[end] == "end"
    length(rng) > 0
  else
    length(rng) == (avg_type == "monthly" ? 12 : (avg_type == "seasonal" ? 4 : 1))*(parse(Float64,a[2])-(parse(Float64,a[1])-1))
  end
end

# to_annual_average(time_series) = squeeze(mean(reshape(time_series,(12,div(size(time_series)[1],12),size(time_series)[2:end]...)),1),1)
to_annual_average(time_series) = dropdims(mean(reshape(time_series,(12,div(size(time_series)[1],12),size(time_series)[2:end]...)),dims=1),dims=1)

function to_array(dset,var_name,averaging=:full)
  if var_name[1] == 'T'
    getproperty(dset,Symbol(var_name * (averaging == :full ? "" : "_$(averaging)"))).values
  else
    var_name_parts = split(var_name,"_")
    if (length(var_name_parts) == 1) || !(var_name_parts[2] in ["clear","cloud"])
      to_array(dset,join(insert!(copy(var_name_parts),2,"clear"),"_"),averaging)+to_array(dset,join(insert!(copy(var_name_parts),2,"cloud"),"_"),averaging)
    else
      if var_name_parts[1] == "N"
        to_array(dset,"LW_"*join(var_name_parts[2:end],"_"),averaging)+to_array(dset,"SW_"*join(var_name_parts[2:end],"_"),averaging)
      else
        getproperty(dset,Symbol(var_name * (averaging == :full ? "" : "_$(averaging)"))).values
      end
    end
  end
end

function to_filtered_feedback_matrix(feedback_matrix,feedback_type)
  nvec = size(feedback_matrix)[1]
  feedback_list = split(string(feedback_type),"_")
  if "local" in feedback_list
    feedback_matrix = Diagonal(diag(feedback_matrix))
  elseif "nonlocal" in feedback_list
    [feedback_matrix[i,i] = 0 for i in 1:nvec]
  end

  if ("tropics" in feedback_list) || ("extratropics" in feedback_list)
    if nvec == 288
      tropics = sort(union((6 + ((1:24)-1)*12),(7 + ((1:24)-1)*12)))
    else
      tropics = sort(union((11 + ((1:48)-1)*24),(12 + ((1:48)-1)*24),(13 + ((1:48)-1)*24),(14 + ((1:48)-1)*24)))
    end
    if "tropics" in feedback_list
      feedback_matrix[:,setdiff(1:nvec,tropics)] = 0
    elseif "extratropics" in feedback_list
      feedback_matrix[:,tropics] = 0
    end
  end
  feedback_matrix
end

function to_dset_patch(dset,lat_i,lon_i,adjacence)
  nlat = length(dset.lat.values)
  nlon = length(dset.lon.values)
  dset.isel(lon = [j%nlon for j in (lon_i-adjacence-1):(lon_i+adjacence-1)],lat = [j%nlat for j in max(lat_i-adjacence-1,0):min(lat_i+adjacence-1,nlat-1)])
end
to_patch(dset,var_name,lat_i,lon_i,adjacence) = to_dset_patch(dset,lat_i,lon_i,adjacence)[var_name].values


function to_vector(lat_lon)
  old_size = size(lat_lon)
  nlat,nlon = old_size[end-1:end]
  if length(old_size) == 2
    reshape(lat_lon,(old_size[1]*old_size[2])),nlat,nlon
  else length(old_size) == 3
    nlat,nlon = old_size[2:3]
    reshape(lat_lon,(old_size[1],old_size[2]*old_size[3])),nlat,nlon
  end
end

to_map(vector,nlat,nlon) = reshape(vector,(size(vector)[1],nlat,nlon))

to_gregory_dict(dset,flux_name) = Dict("T_surf"=>to_annual_average(to_array(dset,"T_surf",:global)),"flux"=>to_annual_average(to_array(dset,flux_name,:global)))
to_climatology(time_series) = dropdims(mean(reshape(time_series,(12,div(size(time_series)[1],12),size(time_series)[2:end]...)),dims=2),dims=2)
function to_anomaly!(data,climatology)
  for month in 1:12
    time_range = month:12:(size(data)[1])
    data[time_range,:,:] .-= reshape(climatology[month,:,:],(1,size(climatology[month,:,:])...))
  end
end

to_index(lat_i,lon_i) = lat_i + (lon_i-1)*12

to_min_years(adjacence,tropics_width,nlat,nlon) = ((adjacence == :all) || (tropics_width == :all)) ? nlat*nlon : (2adjacence+1)^2 + 2tropics_width * nlon

function to_feedback_dict(xs,ys,nlat,nlon,adjacence,tropics_width,polar_mask)
  xs .-= mean(xs,dims=1)
  ys .-= mean(ys,dims=1)

  nvec = nlat*nlon
  feedback = zeros(nvec,nvec)

  for lat_i in 1:nlat, lon_i in 1:nlon
    i_s = [to_index(([mod(adj_lat_i,nlat),mod(adj_lon_i,nlon)]+1)...) for adj_lon_i in ((lon_i-1)-adjacence):((lon_i-1)+adjacence), adj_lat_i in max((lat_i-1)-adjacence,0):min((lat_i-1)+adjacence,nlat-1)]
    if abs(lat_i-6.5) < 6-polar_mask
      band_lat_is = tropics_width == :special ? [3,8] : (6.5+[-.5,.5]*(2*tropics_width-1))
      i_s = union([collect(band_lat_i+(0:23)*12) for band_lat_i in (band_lat_is[1]):(band_lat_is[2])]...,i_s)
    else
      i_s = union(i_s)
    end
    i_s = [Int(i) for i in i_s]
    patch_xs = xs[:,i_s]
    patch_ys = ys[:,to_index(lat_i,lon_i)]
    # if month != nothing
    #   patch_xs = patch_xs[month:12:end,:]
    #   patch_ys = patch_ys[month:12:end,:]
    # end
    var_x = ((patch_xs')*patch_xs)^-1
    patch_feedback = (var_x*patch_xs')*patch_ys

    feedback[to_index(lat_i,lon_i),i_s] = patch_feedback
  end
  # feedback = reshape(feedback,(nvec,nvec))
  y_hats = xs*feedback

  # if full
  #   feedback = feedback'
  # end

  Dict("feedback"=>feedback) #,"R2"=>R2,"est_s2"=>est_s2,"feedback_unct"=>feedback_unct,"err"=>err)
end

function to_feedback_dict(xs,ys,multiple=false)
  d = -1
  # feedback = ((xs'*ys)/(xs'*xs))
  if (length(size(xs)) > 1) && (size(xs)[2] != 1)
    nvec = size(xs)[2]
    if (length(size(ys)) > 1) && (size(ys)[2] != 1) && !multiple
      feedback = [xs[:,i]'*ys[:,i]/(xs[:,i]'*xs[:,i]) for i in 1:nvec]
      y_hats = reshape(feedback,(1,nvec)).*xs

      R2 = (sum(y_hats.^2,dims=1)./sum(ys.^2,dims=1))' # ((y_hats'*y_hats)/(ys'*ys))
      err = ys-y_hats
      df = (size(xs)[end] == 1) ? size(xs)[end-1]-2 : size(xs)[end]-2
      est_s2 = (sum(err.^2,dims=1)/df)'
    else
      if multiple
        var_x = ((xs')*xs)^-1
        feedback = (var_x*xs')*ys
        if multiple
          feedback = feedback'
        end
        y_hats = multiple ? xs*feedback' : sum(reshape(feedback,(1,nvec)).*xs,dims=2)
        err = ys-y_hats
        df = (size(xs)[end] == 1) ? size(xs)[end-1]-2 : size(xs)[1]-(size(xs)[2]+1)
        est_s2 = (sum(err.^2,dims=1)/df)'
        feedback_unct = zeros(nvec,nvec)
        R2 = (sum(y_hats.^2,dims=1)./sum(ys.^2,dims=1))' # ((y_hats'*y_hats)/(ys'*ys))
        d = sum((err[2:end,:] - err[1:(end-1),:]).^2,dims=1) ./ sum(err .^2,dims=1)
       else
        feedback = [xs[:,i]'*ys/(xs[:,i]'xs[:,i]) for i in 1:nvec]
        y_hats = reshape(feedback,(1,nvec)).*xs
        R2 = (sum(y_hats.^2,dims=1)./sum(ys.^2,dims=1))' # ((y_hats'*y_hats)/(ys'*ys))
        err = ys-y_hats
        df = size(xs)[1] - 2
        est_s2 = (sum(err.^2,dims=1)/df)'
      end
    end
  else
    feedback = xs'*ys./(xs'*xs)
    y_hats = feedback.*xs

    R2 = (sum(y_hats.^2,dims=1)./sum(ys.^2,dims=1))' # ((y_hats'*y_hats)/(ys'*ys))
    err = ys-y_hats
    df = (size(xs)[end] == 1) ? size(xs)[end-1]-2 : size(xs)[end]-2
    est_s2 = (sum(err.^2,dims=1)/df)'
  end

  if (length(size(xs)) > 1) && (size(xs)[2] != 1)
    if (length(size(ys)) > 1) && (size(ys)[2] != 1)
      if !multiple
        feedback_unct = [quantile(TDist(df),0.975) * ((est_s2[i]./(xs[:,i]'*xs[:,i])).^.5) for i in 1:nvec]
      end
    else
      feedback_unct = quantile(TDist(df),0.975) * sqrt.(diag(var_x) * est_s2)
    end
  else
    # println(df)
    feedback_unct = 0.0 #quantile(TDist(df),0.975) * ((est_s2./(xs'*xs)).^.5)
  end
  Dict("feedback" => copy(feedback),"R2" => R2,"est_s2"=>est_s2,"feedback_unct"=>feedback_unct,"err"=>err,"d"=>d)
end

function to_annual_seasonal_and_monthly_feedback_dict(xs,ys,multiple=false)
  length_xs = size(xs)[1]
  seasonal_xs = to_mean_by(xs[3:(length_xs-1),:],3)
  seasonal_ys = to_mean_by(ys[3:(length_xs-1),:],3)
  Dict(
    "annual" => to_feedback_dict(to_annual_average(xs),to_annual_average(ys),multiple),
    "seasonal" => [to_feedback_dict(seasonal_xs[i:4:end,:],seasonal_ys[i:4:end,:],multiple) for i in 1:4],
    "monthly" => [to_feedback_dict(xs[i:12:end,:],ys[i:12:end,:],multiple) for i in 1:12],
    "month_agnostic" => to_feedback_dict(xs,ys,multiple)
  )
end

function to_annual_seasonal_and_monthly_feedback_dict(xs,ys,nlat,nlon,adjacence,tropics_width,polar_mask)
  Dict(
    "annual" => to_feedback_dict(to_annual_average(xs),to_annual_average(ys),nlat,nlon,adjacence,tropics_width,polar_mask),
    "monthly" => [to_feedback_dict(xs[i:12:end,:],ys[i:12:end,:],nlat,nlon,adjacence,tropics_width,polar_mask) for i in 1:12],
    "month_agnostic" => to_feedback_dict(xs,ys,nlat,nlon,adjacence,tropics_width,polar_mask)
  )
end

# 247: tropical ascent
# 9: extratropics
# 32: subsidence
function to_individual_netcdf(model_name,grid,fn_is)
  feedback_dict = to_feedback_dict(to_feedback_fn(model_name,"control",:full,grid))
  new_dset = xr.Dataset(Dict(),coords=Dict("lat" => grid_dict[grid]["lats"], "lon" => grid_dict[grid]["lons"].-180))
  new_dset.lat.attrs = Dict("units"=>"degrees_north")
  new_dset.lon.attrs = Dict("units"=>"degrees_east")
  for (var_name,var_dict) in feedback_dict
    for fn_i in fn_is
      # new_dset = new_dset.assign(; [(Symbol(var_name * "_$(fn_i)"),(("lat","lon"),reshape(var_dict["annual"]["feedback"][:,fn_i],(length(grid_dict[grid]["lats"]),length(grid_dict[grid]["lons"])))))]...)
      new_dset = new_dset.assign(; [(Symbol(var_name * "_$(fn_i)"),(("lat","lon"),reshape(sum([mean(x -> x["feedback"],var_dict["seasonal"])[:,mod((fn_i+12i),288)] for i in -1:1]),(length(grid_dict[grid]["lats"]),length(grid_dict[grid]["lons"])))))]...)
      
    end
  end
  new_dset.to_netcdf(to_feedback_fn(model_name,"control","full_selected",grid))
end

# function to_score_1(feedback_type,runs,averaging_type,var_name,grid="orig",years_used=:all,adjacence=:all,tropics_width=:all)
#   if runs == "long_run"
#     trues = [collect(to_feedback_dict(to_feedback_fn(model_name,get_run(model_name),:global_global_time_scales))[var_name])[end][2]["annual"]["feedback"][1] for model_name in all_model_names]
#     ests = [collect(to_feedback_dict(to_feedback_fn(model_name,get_run(model_name),:global_global_time_scales,grid,estimate_type=feedback_type,years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=(tropics_width == :all ? :none : 0)))[var_name])[end][2][averaging_type]["feedback"][1] for model_name in all_model_names]
#   # elseif runs == "abrupt4x_sans_MIROC"
#   #   trues = [collect(to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales))[var_name])[runs != "abrupt4x_starts" ? 2 : 1][2]["annual"]["feedback"][1] for model_name in model_names_abrupt4x_sans_MIROC]
#   #   ests = [collect(to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid,estimate_type=feedback_type))[var_name])[runs != "abrupt4x_starts" ? 2 : 1 ? 2 : 1][2][averaging_type]["feedback"][1] for model_name in model_names_abrupt4x_sans_MIROC]
#   else
#     trues = [collect(to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales))[var_name])[runs != "abrupt4x_starts" ? 2 : 1][2]["annual"]["feedback"][1] for model_name in model_names_abrupt4x]
#     spec_i = (runs != "abrupt4x_starts") ? 2 : 1
#     # println(spec_i)
#     # println(averaging_type)
    # ests = [collect(to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid,estimate_type=feedback_type,years_used=to_true_years_end(model_name,years_used),adjacence=adjacence,tropics_width=tropics_width,polar_mask=(tropics_width == :all ? :none : 0)))[var_name])[spec_i][2][averaging_type]["feedback"][1] for model_name in model_names_abrupt4x]
#   end
#   sqrt(sum((ests-trues).^2)/length(ests))
# end

# function to_score_1(feedback_type,runs,averaging_type,var_name,years_used=:all,adjacence=:all,tropics_width=:all)
#   if feedback_type == :global_local_multiple
#     grid_tmp = (model_name) -> (model_name == "GISSE2R") ? "24x48_equaldist" : "12x24_equaldist"
#   else
#     grid_tmp = (model_name) -> "orig"
#   end
#
#   if runs == "long_run"
#     trues = [collect(to_feedback_dict(to_feedback_fn(model_name,get_run(model_name),:global_global_time_scales))[var_name])[end][2]["annual"]["feedback"][1] for model_name in all_model_names]
#     ests = [collect(to_feedback_dict(to_feedback_fn(model_name,get_run(model_name),:global_global_time_scales,grid_tmp(model_name),estimate_type=feedback_type,years_used=years_used,adjacence=adjacence,tropics_width=tropics_width,polar_mask=(tropics_width == :all ? :none : 0)))[var_name])[end][2][averaging_type]["feedback"][1] for model_name in all_model_names]
#   else
#     trues = [collect(to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales))[var_name])[runs != "abrupt4x_starts" ? 2 : 1][2]["annual"]["feedback"][1] for model_name in model_names_abrupt4x]
#     spec_i = (runs != "abrupt4x_starts") ? 2 : 1
#     ests = [collect(to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid_tmp(model_name),estimate_type=feedback_type,years_used=to_true_years_end(model_name,years_used),adjacence=adjacence,tropics_width=tropics_width,polar_mask=(tropics_width == :all ? :none : 0)))[var_name])[spec_i][2][averaging_type]["feedback"][1] for model_name in model_names_abrupt4x]
#   end
#   sqrt(sum((ests-trues).^2)/length(ests))
# end


function to_special_score_1(flux,period,years_used)
  grid_tmp = (model_name) -> "12x24_equaldist" # (model_name == "GISSE2R") ? "24x48_equaldist" : "12x24_equaldist"
  years_used_tmp = (model_name) -> to_years(model_name,"control") == years_used ? :all : years_used

  if period == :change
    early_trues = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid_tmp(model_name)))[flux]["2 to 20"]["special_annual"]["feedback"][1] for model_name in model_names_abrupt4x]
    early_ests = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid_tmp(model_name),estimate_type=:full,years_used=years_used_tmp(model_name)))[flux]["2 to 20"]["special_monthly"]["feedback"][1] for model_name in model_names_abrupt4x]
    late_trues = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid_tmp(model_name)))[flux]["21 to end"]["special_annual"]["feedback"][1] for model_name in model_names_abrupt4x]
    late_ests = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid_tmp(model_name),estimate_type=:full,years_used=years_used_tmp(model_name)))[flux]["21 to end"]["special_monthly"]["feedback"][1] for model_name in model_names_abrupt4x]
    trues = late_trues - early_trues
    ests = late_ests - early_ests
  else
    time_scale = period == :early ? "2 to 20" : "21 to end"
    trues = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid_tmp(model_name)))[flux][time_scale]["special_annual"]["feedback"][1] for model_name in model_names_abrupt4x]
    ests = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid_tmp(model_name),estimate_type=:full,years_used=years_used_tmp(model_name)))[flux][time_scale]["special_monthly"]["feedback"][1] for model_name in model_names_abrupt4x]
  end
  sqrt(sum((ests-trues).^2)/length(ests))
end

function to_score_1(flux,period,years_used)
  grid_tmp = (model_name) -> "12x24_equaldist" # (model_name == "GISSE2R") ? "24x48_equaldist" : "12x24_equaldist"
  years_used_tmp = (model_name) -> to_years(model_name,"control") == years_used ? :all : years_used

  if period == :change
    early_trues = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales))[flux]["2 to 20"]["annual"]["feedback"][1] for model_name in model_names_abrupt4x]
    early_ests = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid_tmp(model_name),estimate_type=:full,years_used=years_used_tmp(model_name)))[flux]["2 to 20"]["monthly"]["feedback"][1] for model_name in model_names_abrupt4x]
    late_trues = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales))[flux]["21 to end"]["annual"]["feedback"][1] for model_name in model_names_abrupt4x]
    late_ests = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid_tmp(model_name),estimate_type=:full,years_used=years_used_tmp(model_name)))[flux]["21 to end"]["monthly"]["feedback"][1] for model_name in model_names_abrupt4x]
    trues = late_trues - early_trues
    ests = late_ests - early_ests
  else
    time_scale = period == :early ? "2 to 20" : "21 to end"
    trues = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales))[flux][time_scale]["annual"]["feedback"][1] for model_name in model_names_abrupt4x]
    ests = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid_tmp(model_name),estimate_type=:full,years_used=years_used_tmp(model_name)))[flux][time_scale]["monthly"]["feedback"][1] for model_name in model_names_abrupt4x]
  end
  sqrt(sum((ests-trues).^2)/length(ests))
end

function to_score_2_per_model(model_name,run_type,var_name,averaging_type,grid,feedback_type,regridded_estimate,i,years_used=:all,adjacence=:all,tropics_width=:all)
  fd_fn = to_finite_difference_fn(model_name,run_type,grid)
  fdd_fn = to_finite_difference_difference_fn(model_name,run_type,grid,feedback_type,regridded_estimate,years_used=years_used,adjacence=adjacence,tropics_width=tropics_width)
  # println(fd_fn)
  # println(fdd_fn)
  fd_dset = to_dset(fd_fn)
  fdd_dset = to_dset(fdd_fn)
  global_value = getproperty(fd_dset,Symbol(join(["T_surf",i,"global"],"_"))).values[1]
  diff_values = getproperty(fdd_dset,Symbol(join([var_name,averaging_type,i],"_"))).values
  to_score_2_per_model(diff_values./global_value,grid)
end

function to_score_2_per_model_change(model_name,run_type,var_name,averaging_type,grid,feedback_type,regridded_estimate,years_used=:all,adjacence=:all,tropics_width=:all)
  fd_fn = to_finite_difference_fn(model_name,run_type,grid)
  fdd_fn = to_finite_difference_difference_fn(model_name,run_type,grid,feedback_type,regridded_estimate,years_used=years_used,adjacence=adjacence,tropics_width=tropics_width)
  # println(fd_fn)
  # println(fdd_fn)
  fd_dset = to_dset(fd_fn)
  fdd_dset = to_dset(fdd_fn)
  global_value_1 = getproperty(fd_dset,Symbol(join(["T_surf",1,"global"],"_"))).values[1]
  diff_values_1 = getproperty(fdd_dset,Symbol(join([var_name,averaging_type,1],"_"))).values
  global_value_2 = getproperty(fd_dset,Symbol(join(["T_surf",2,"global"],"_"))).values[1]
  diff_values_2 = getproperty(fdd_dset,Symbol(join([var_name,averaging_type,2],"_"))).values
  to_score_2_per_model(diff_values_2./global_value_2 - diff_values_1./global_value_1,grid)
end


function to_special_score_2_per_model_change(model_name,run_type,var_name,averaging_type,grid,feedback_type,regridded_estimate,years_used=:all,adjacence=:all,tropics_width=:all)
  fd_fn = to_finite_difference_fn(model_name,run_type,grid)
  fdd_fn = to_finite_difference_difference_fn(model_name,run_type,grid,feedback_type,regridded_estimate,years_used=years_used,adjacence=adjacence,tropics_width=tropics_width)
  # println(fd_fn)
  # println(fdd_fn)
  fd_dset = to_dset(fd_fn)
  fdd_dset = to_dset(fdd_fn)
  global_value_1 = getproperty(fd_dset,Symbol(join(["T_surf",1,"global"],"_"))).values[1]
  diff_values_1 = getproperty(fdd_dset,Symbol(join([var_name,averaging_type,1],"_"))).values
  global_value_2 = getproperty(fd_dset,Symbol(join(["T_surf",2,"global"],"_"))).values[1]
  diff_values_2 = getproperty(fdd_dset,Symbol(join([var_name,averaging_type,2],"_"))).values
  to_special_score_2_per_model(diff_values_2[5:end,:]./global_value_2 - diff_values_1[5:end,:]./global_value_1,grid)
end


# function to

# to_score_2_per_model(to_dset(fdd_fn)[join([var_name,averaging_type,i],"_")].values,to_dset(fd_fn)[join(["T_surf",i,"global"],"_")].values[1])

to_score_2_per_model(normalized_diff_values,grid) = sqrt(to_spatial_average((normalized_diff_values.^2),grid_dict[grid]["lats"]))

function to_special_score_2_per_model(model_name,run_type,var_name,averaging_type,grid,feedback_type,regridded_estimate,i,years_used=:all,adjacence=:all,tropics_width=:all)
  fd_fn = to_finite_difference_fn(model_name,run_type,grid)
  fdd_fn = to_finite_difference_difference_fn(model_name,run_type,grid,feedback_type,regridded_estimate,years_used=years_used,adjacence=adjacence,tropics_width=tropics_width)

  fd_dset = to_dset(fd_fn)
  fdd_dset = to_dset(fdd_fn)
  global_value = getproperty(fd_dset,Symbol(join(["T_surf",i,"global"],"_"))).values[1]
  diff_values = getproperty(fdd_dset,Symbol(join([var_name,averaging_type,i],"_"))).values
  to_special_score_2_per_model(diff_values[5:end,:]./global_value,grid)
end

to_special_score_2_per_model(normalized_diff_values,grid) = sqrt(to_special_spatial_average((normalized_diff_values.^2),grid_dict[grid]["lats"][5:end]))

# function to_score_2(feedback_type,runs,averaging_type,var_name,grid,years_used=:all,adjacence=:all,tropics_width=:all)
#   if runs == "long_run"
#     mean([to_score_2_per_model(model_name,get_run(model_name),var_name,averaging_type,grid,feedback_type,length(time_scales_dict[get_run(model_name)]),years_used,adjacence,tropics_width) for model_name in all_model_names])
#   else
#     mean([to_score_2_per_model(model_name,"abrupt4x",var_name,averaging_type,grid,feedback_type,runs != "abrupt4x_starts" ? 2 : 1,to_true_years_end(model_name,years_used),adjacence,tropics_width) for model_name in model_names_abrupt4x])
#   end
# end

# gridless
function to_score_2(feedback_type,runs,averaging_type,var_name,years_used=:all,adjacence=:all,tropics_width=:all)

  if runs == "long_run"
    mean([to_score_2_per_model(model_name,get_run(model_name),var_name,averaging_type,"12x24_equaldist",feedback_type,length(time_scales_dict[get_run(model_name)]),years_used,adjacence,tropics_width) for model_name in all_model_names])
  else
    mean([to_score_2_per_model(model_name,"abrupt4x",var_name,averaging_type,"12x24_equaldist",feedback_type,runs != "abrupt4x_starts" ? 2 : 1,to_true_years_end(model_name,years_used),adjacence,tropics_width) for model_name in model_names_abrupt4x])
  end
end

# function to_score_3(feedback_type,runs,averaging_type,var_name,grid,years_used=:all,adjacence=:all,tropics_width=:all)
#   trues = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales))[var_name] for model_name in model_names_abrupt4x]
#   true_diffs = [t[time_scales_dict["abrupt4x"][2]]["annual"]["feedback"][1] - t[time_scales_dict["abrupt4x"][1]]["annual"]["feedback"][1] for t in trues]
#   ests = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid,estimate_type=feedback_type))[var_name] for model_name in model_names_abrupt4x]
#   est_diffs = [t[time_scales_dict["abrupt4x"][2]][averaging_type]["feedback"][1] - t[time_scales_dict["abrupt4x"][1]][averaging_type]["feedback"][1] for t in ests]
#   sqrt(sum((est_diffs - true_diffs).^2)/length(est_diffs))
# end

# gridless
function to_score_3(feedback_type,runs,averaging_type,var_name,years_used=:all,adjacence=:all,tropics_width=:all)
  if feedback_type == :global_local_multiple
    grid_tmp = (model_name) -> (model_name == "GISSE2R") ? "24x48_equaldist" : "12x24_equaldist"
  else
    grid_tmp = (model_name) -> "orig"
  end

  trues = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales))[var_name] for model_name in model_names_abrupt4x]
  true_diffs = [t[time_scales_dict["abrupt4x"][2]]["annual"]["feedback"][1] - t[time_scales_dict["abrupt4x"][1]]["annual"]["feedback"][1] for t in trues]
  ests = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid_tmp(model_name),estimate_type=feedback_type))[var_name] for model_name in model_names_abrupt4x]
  est_diffs = [t[time_scales_dict["abrupt4x"][2]][averaging_type]["feedback"][1] - t[time_scales_dict["abrupt4x"][1]][averaging_type]["feedback"][1] for t in ests]
  sqrt(sum((est_diffs - true_diffs).^2)/length(est_diffs))
end

to_added_var(dset,var_name,dims,values) = dset.assign(; [(Symbol(var_name),(dims,values))]...)

function to_true_years_end(model_name,years_used)
  true_length = length(to_dset(to_anomalized_fn(model_name,"control"))["time"])/12
  years_used == :all ? :all : (true_length > years_used ? years_used : (true_length == years_used ? :all : nothing))
end

# function polyfit_exps(x,y,n,flag=false)
#   A = parse(Float64,[ parse(Float64,x[i])^p for i = 1:length(x), p = (flag == :no_constant ? 1 : 0):n ])
#   A \ y
# end
#
# function uncertainty(x,y,y_est,n)
#   SSE = sum((y .- y_est).^2)
#   MSE = SSE / (n - 2)
#   quantile(TDist(n-2),.975) * sqrt(MSE / sum((x-mean(x)).^2))
# end


# time_average

# climatology

# grid_to_vector

# vector_to_grid
