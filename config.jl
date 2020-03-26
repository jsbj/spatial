fn_prefix = "../" #"../../../../../project2/abbot/jsbj/spatial/"
annual_time_scales = ["1 to 20","21 to 150","151 to end"]
freqs = ["LW","SW","N"]
skies = ["_clear","_cloud",""]
var_names = ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]


model_names_abrupt4x = ["CCSM3","CESM104","ECEARTH","GISSE2R","HadCM3L","IPSLCM5A","MIROC","MPIESM12"] # "GFDLESM2M",
model_names_abrupt4x_sans_MIROC = setdiff(model_names_abrupt4x,["MIROC"])
model_names_1pct2x = ["GFDLCM3","GFDLESM2M","MIROC"]
model_names_rcp85 = ["ECEARTH"]
new_abrupt4x = ["ECEARTH","GFDLESM2M","MIROC"]
# all_model_names = ["CCSM3","CESM104","ECEARTH","GFDLCM3","GFDLESM2M","GISSE2R","HadCM3L","IPSLCM5A","MIROC","MPIESM12"] # ["HadCM3L"]# ["CCSM3","CESM104","ECEARTH","GFDLCM3","GFDLESM2M","GISSE2R","HadCM3L","IPSLCM5A","MIROC","MPIESM12"] # ["HadCM3L"] #,"HadCM3L","IPSLCM5A","MIROC","MPIESM12"] # ["CCSM3","CESM104","ECEARTH","GFDLCM3","GFDLESM2M","GISSE2R","HadCM3L","IPSLCM5A","MIROC","MPIESM12"]

model_name_dict = OrderedDict(
  "abrupts" => setdiff(setdiff(model_names_abrupt4x,model_names_1pct2x),model_names_rcp85),
  "ramps" => model_names_1pct2x,
  "rcps" => model_names_rcp85
)

# time_scales_dict = Dict(
#   "abrupt2x" => ["1 to 20","21 to 150","151 to end"],
#   "abrupt4x" => ["1 to 20","21 to 150","151 to end"],
#   "abrupt6x" => ["1 to 20","21 to 150","151 to end"],
#   "abrupt8x" => ["1 to 20","21 to 150","151 to end"],
#   "abrupt16x" => ["1 to 20","21 to 150","151 to end"],
#   "abrupt32x" => ["1 to 20","21 to 150","151 to end"],
#   "1pct2x" => ["71 to end"],
#   "1pct4x" => ["141 to end"],
#   "rcp85" => ["301 to end"]
# )

model_names_abrupt4x = ["CESM104","CNRMCM61","GISSE2R","HadCM3L","IPSLCM5A","MPIESM12"] # ["CNRMCM61"] # ["CCSM3"] #["CCSM3","CESM104","GISSE2R","HadCM3L","IPSLCM5A","MPIESM12"]
# model_names_abrupt4x = ["GISSE2R"]
all_model_names = model_names_abrupt4x
# all_model_names = ["MPIESM12"]
# model_names_1pct2x = ["GFDLCM3","GFDLESM2M"]
# all_model_names = ["CCSM3","CESM104","GFDLCM3","GFDLESM2M","GISSE2R","HadCM3L","IPSLCM5A","MPIESM12"]
# model_name_dict = Dict(
#   "abrupts" => model_names_abrupt4x,
#   "ramps" => model_names_1pct2x
# )
#

time_scales_dict = Dict(
  "control" => ["1 to end"],
  "abrupt2x" => ["2 to 20","21 to end"], # ["6 to 50","51 to end"], #["1 to 20","21 to 150","151 to end"],
  "abrupt4x" => ["2 to 20","21 to end"], # ["6 to 50","51 to end"], # ["3 to 50","51 to end"], # ["1 to 50","51 to end"], #["1 to 20","21 to 150","151 to end"],
  "abrupt6x" => ["2 to 20","21 to end"], # ["6 to 50","51 to end"], #["1 to 20","21 to 150","151 to end"],
  "abrupt8x" => ["2 to 20","21 to end"], # ["6 to 50","51 to end"], #["1 to 20","21 to 150","151 to end"],
  "abrupt16x" => ["2 to 20","21 to end"], # ["6 to 50","51 to end"], #["1 to 20","21 to 150","151 to end"],
  "abrupt32x" => ["2 to 20","21 to end"], # ["6 to 50","51 to end"], #["1 to 20","21 to 150","151 to end"],
  "1pct2x" => ["71 to end"], # ["71 to 120","121 to end"], #
  "1pct4x" => ["141 to end"], # ["141 to 190","191 to end"], #
  "rcp85" => ["301 to end"] # ["301 to 350","351 to end"] #
)

grid_dict = OrderedDict(
  "12x24_equaldist" => Dict("lats"=>equal_distance_lats(12),"lons"=>equal_lons(24)) #,
  # "24x48_equaldist" => Dict("lats"=>equal_distance_lats(24),"lons"=>equal_lons(48)) #,
  # "16x32_equaldist" => Dict("lats"=>equal_distance_lats(16),"lons"=>equal_lons(32)) #,
  # "12x24_equalarea" => Dict("lats"=>equal_area_lats(12),"lons"=>equal_lons(24)),
  # "24x12_equaldist" => Dict("lats"=>equal_distance_lats(24),"lons"=>equal_lons(12)),
  # "24x24_equaldist" => Dict("lats"=>equal_distance_lats(24),"lons"=>equal_lons(24)),
  # "24x36_equaldist" => Dict("lats"=>equal_distance_lats(24),"lons"=>equal_lons(36)),
  # "12x12_equaldist" => Dict("lats"=>equal_distance_lats(12),"lons"=>equal_lons(12)),
  # "6x12_equaldist" => Dict("lats"=>equal_distance_lats(6),"lons"=>equal_lons(12)),
  # "12x36_equaldist" => Dict("lats"=>equal_distance_lats(12),"lons"=>equal_lons(36)),
  # "8x16_equaldist" => Dict("lats"=>equal_distance_lats(8),"lons"=>equal_lons(16)),
  # "16x16_equaldist" => Dict("lats"=>equal_distance_lats(16),"lons"=>equal_lons(16)),
  # "8x32_equaldist" => Dict("lats"=>equal_distance_lats(8),"lons"=>equal_lons(32)),
  # "6x36_equaldist" => Dict("lats"=>equal_distance_lats(6),"lons"=>equal_lons(36)),
  # "6x24_equaldist" => Dict("lats"=>equal_distance_lats(6),"lons"=>equal_lons(24)),
  # "12x16_equaldist" => Dict("lats"=>equal_distance_lats(12),"lons"=>equal_lons(16)),
  # "8x24_equaldist" => Dict("lats"=>equal_distance_lats(8),"lons"=>equal_lons(24)),
  # "12x32_equaldist" => Dict("lats"=>equal_distance_lats(12),"lons"=>equal_lons(32)),
  # "16x24_equaldist" => Dict("lats"=>equal_distance_lats(16),"lons"=>equal_lons(24)),
  # "10x28_equaldist" => Dict("lats"=>equal_distance_lats(10),"lons"=>equal_lons(28)),
  # "12x28_equaldist" => Dict("lats"=>equal_distance_lats(12),"lons"=>equal_lons(28)),
  # "10x24_equaldist" => Dict("lats"=>equal_distance_lats(10),"lons"=>equal_lons(24))
)

reduceds = [(:all,:all,:none)] #,(2,1,0),(2,2,0),(2,3,0),(2,4,0),(2,5,0)] # [(2,4,0),(1,4,0)] # [(4,3,0)] :special #[(:all,:all,:none),(2,1,0),(2,2,0),(2,3,0),(2,4,0),(2,5,0),(1,1,0),(1,2,0),(1,3,0),(1,4,0),(1,5,0)] # [1,2,4,8,:all] [(4,3,0)] #

years_useds = [:all] # [5,10,20,50,100,200,350,500,750,1000,:all]
# years_useds = [350,500,750,1000,:all]

# function equal_distance1(nlat)
#   new_range = asin.(-1:(2/nlat):1)*180/π
#   (new_range[1:(end-1)] + new_range[2:end])/2
# end

# function equal_distance2(nlat)
#   new_range = asin.(-1:(2/nlat):1)*180/π
#   lats = zeros(nlat)
#   diff_tmp = (new_range[2] - new_range[1])/2
#   lats[1] = new_range[1]+diff_tmp
#   for i in 2:Int(nlat/2)
#     lats[i] = new_range[i]+diff_tmp
#     diff_tmp = (new_range[i+1] - new_range[i]) - diff_tmp
#   end
#
#   lats[(Int(nlat/2)+1):end] = -reverse(lats[1:(Int(nlat/2))])
#   lats
# end

function reset_time_scales(which_way)
  if which_way == "to three"
    # mv netcdfs/orig/finite_difference/*abrupt4x* to_two/netcdfs/orig/finite_difference/
    # mv netcdfs/12x24_equaldist/finite_difference/*abrupt4x* to_two/netcdfs/12x24_equaldist/finite_difference/
    # mv netcdfs/12x24_equaldist/finite_difference_difference/*abrupt4x* to_two/netcdfs/12x24_equaldist/finite_difference_difference/
    # mv jlds/orig/*abrupt4x* to_two/jlds/orig/
    # mv jlds/orig/estimated/*abrupt4x* to_two/jlds/orig/estimated/
    # mv jlds/12x24_equaldist/estimated/*abrupt4x* to_two/jlds/12x24_equaldist/estimated/
    # mv to_three/netcdfs/orig/finite_difference/* netcdfs/orig/finite_difference/
    # mv to_three/netcdfs/12x24_equaldist/finite_difference/* netcdfs/12x24_equaldist/finite_difference/
    # mv to_three/netcdfs/12x24_equaldist/finite_difference_difference/* netcdfs/12x24_equaldist/finite_difference_difference/
    # mv to_three/jlds/orig/*.nc jlds/orig/
    # mv to_three/jlds/orig/estimated/* jlds/orig/estimated/
    # mv to_three/jlds/12x24_equaldist/estimated/* jlds/12x24_equaldist/estimated/
  elseif which_way == "to two"
    # mv netcdfs/orig/finite_difference/*abrupt4x* to_three/netcdfs/orig/finite_difference/
    # mv netcdfs/12x24_equaldist/finite_difference/*abrupt4x* to_three/netcdfs/12x24_equaldist/finite_difference/
    # mv netcdfs/12x24_equaldist/finite_difference_difference/*abrupt4x* to_three/netcdfs/12x24_equaldist/finite_difference_difference/
    # mv jlds/orig/*abrupt4x* to_three/jlds/orig/
    # mv jlds/orig/estimated/*abrupt4x* to_three/jlds/orig/estimated/
    # mv jlds/12x24_equaldist/estimated/*abrupt4x* to_three/jlds/12x24_equaldist/estimated/
    # mv to_two/netcdfs/orig/finite_difference/* netcdfs/orig/finite_difference/
    # mv to_two/netcdfs/12x24_equaldist/finite_difference/* netcdfs/12x24_equaldist/finite_difference/
    # mv to_two/netcdfs/12x24_equaldist/finite_difference_difference/* netcdfs/12x24_equaldist/finite_difference_difference/
    # mv to_two/jlds/orig/* jlds/orig/
    # mv to_two/jlds/orig/estimated/* jlds/orig/estimated/
    # mv to_two/jlds/12x24_equaldist/estimated/* jlds/12x24_equaldist/estimated/
  end
end


  # if !isempty(glob("netcdfs/*/finite_difference{,_difference}/*"))
  #   Shell.run("rm netcdfs/*/finite_difference{,_difference}/*")
  # end
  # if !isempty(glob("jlds/*/*_time_scales.jld"))
  #   Shell.run("rm jlds/*/*_time_scales.jld")
  # end
  # if !isempty(glob("jlds/*/estimated/*abrupt4x*_time_scales.jld"))
  #   Shell.run("rm jlds/*/estimated/*_time_scales.jld")
  # end
