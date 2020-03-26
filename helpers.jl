forced_runs(model_name,forcing) = sort(collect(keys(models[model_name][forcing])),by = x -> Int(float(match(r"\d+",x).match)))
control_and_forced_runs(model_name,forcing) = vcat(["control"],forced_runs(model_name,forcing))
forced_runs(model_name) = intersect(vcat([forced_runs(model_name,forcing) for forcing in intersect(collect(keys(models[model_name])),["abrupts","ramps","rcps"])]...),["abrupt4x","rcp85","1pct2x"])
# forced_runs(model_name) = intersect(vcat([forced_runs(model_name,forcing) for forcing in intersect(collect(keys(models[model_name])),["abrupts","ramps","rcps"])]...),["abrupt2x","abrupt4x","abrupt8x","abrupt16x","abrupt32x","rcp85","1pct2x"])
all_runs(model_name) = vcat(["control"], forced_runs(model_name))
equal_lons(nlon) = collect((180/nlon):(360/nlon):(360-180/nlon))
equal_distance_lats(nlat) = collect(-(90-90/nlat):(180/nlat):(90-90/nlat))
equal_area_lats(nlat) = asin.(-(1-1/nlat):(2/nlat):(1-1/nlat))*180/π
vars(dset) = setdiff(dset.keys(),["lat","lon","time"])
get_run(model_name) = (model_name == "ECEARTH" ? "rcp85" : (model_name in model_names_1pct2x ? "1pct2x" : "abrupt4x"))

function midway(var_list)
  args_list = join([string(a)[1] == ':' ? string(a)[2:end] : a for a in var_list],":")
  
  cd("scripts")
  println(`sbatch --export=ARGS_LIST="$(args_list)" code.run`)
  # b = run(`sbatch --export=ARGS_LIST="$(args_list)" code.run`)
  cd("..")
  # b
end

# R2(x,y) =   /sum((y.-mean(y)).^2)