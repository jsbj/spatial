# future things to add
# - do ts instead of tas

include("load.jl")
rc("font", size=20)
rc("xtick.major", pad=8)
rc("ytick.major", pad=8)
rc("axes", linewidth=1)

# ion()
ioff()
mean_by(ar,i) = dropdims(mean(reshape(ar,(i,Int(length(ar)/i))),dims=1)',dims=2)

# to_annual_average(time_series) = squeeze(mean(reshape(time_series,(12,div(size(time_series)[1],12),size(time_series)[2:end]...)),1),1)
to_annual_average(time_series) = dropdims(mean(reshape(time_series,(12,div(size(time_series)[1],12),size(time_series)[2:end]...)),dims=1),dims=1)

function fancy_average_maker(lngth,stps)
  list_is = unique(Int.(ceil.(exp.(0:(log(lngth)/(stps-1)):log(lngth)))))
  ranges = [(i == 1 ? 1 : (list_is[i-1]+1)):list_is[i] for i in 1:length(list_is)]
  x(list) = vcat([mean(list[range]) for range in ranges]...)
end

function feedback_scale(model_name,method,flux,scale_i,averaging="seasonal",area=:all)
  dset = to_dset(to_anomalized_fn(model_name,"abrupt4x","12x24_equaldist"))
  Ts = to_gregory_dict(dset,"N")["T_surf"]
  fa = fancy_average_maker(1000,50)
  fancy_T = fa(Ts) # fancy_average(Ts)

  if method == :true
    monthly_flux = to_array(to_dset(to_anomalized_fn(model_name,"abrupt4x","12x24_equaldist")),flux,"$(area == :special ? "special_" : "")global")
    fluxes = to_annual_average(monthly_flux)
  else
    if method != :full_tropics
      efn = to_estimated_fn(model_name,"abrupt4x",(method in [:global,:local] ? Symbol("local_$(method)") : method),"12x24_equaldist")
    else
      efn = to_estimated_fn(model_name,"abrupt4x",:full,"12x24_equaldist",tropics=true)
    end
    # println(efn)
    fluxes = to_array(to_dset(efn),"$(flux)_$(averaging)","$(area == :special ? "special_" : "")global")
  end
  fancy_flux = fa(fluxes)

  if scale_i == 1
    F,l = polyfit_exps(fancy_T[2:15],fancy_flux[2:15],1)
  elseif scale_i == 2
    F,l = polyfit_exps(fancy_T[16:end],fancy_flux[16:end],1)
  else
    F,l_1 = polyfit_exps(fancy_T[2:15],fancy_flux[2:15],1)
    F,l_2 = polyfit_exps(fancy_T[16:end],fancy_flux[16:end],1)
    l = l_2 - l_1
  end
  l
end

# function fancy_average(list)
#
#
#
#   annual_list = to_annual_average(list)
#   years = Int(floor(length(annual_list)/100)*100)
#   vcat([mean(list[i:i+11]) for i in 1:24],[mean(list[i:i+11]) for i in 25:6:120],annual_list[10:20],mean_by(annual_list[21:50],5),mean_by(annual_list[51:200],10),mean_by(annual_list[201:years],100))
# end

function fancy_average(list)
  years = Int(floor(length(list)/200)*200)
  vcat(list[2:10],mean_by(list[11:30],2),mean_by(list[31:40],5),mean_by(list[41:200],10),mean_by(list[201:600],100),mean_by(list[601:years],200))
end

function polyfit_exps(x,y,n,flag=false)
  A = float([ float(x[i])^p for i = 1:length(x), p = (flag == :no_constant ? 1 : 0):n ])
  A \ y
end
rms(list) = round(sqrt(sum((1/length(list))*list.^2)),digits=2)
rms_line(feedback_dict,true_feedback_dict,area,period,averaging) = join([join([rpad(rms(feedback_dict[period][method][averaging][area][flux]-true_feedback_dict[period][area][flux]),4,"0") for method in [:global, :local, :full]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"
rms_line(feedback_dict,true_feedback_dict,area,period) = join([join([rpad(rms(feedback_dict[period][method][averaging][area][flux]-true_feedback_dict[period][area][flux]),4,"0") for (method,averaging) in [(:global,"month_agnostic"), (:local,"annual"), (:full,"monthly")]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"

# error_line(time_i::Number) = join([join([rpad(round(to_score_2_per_model("mmm","abrupt4x",flux,((method == :full) ? "monthly" : "annual"),grid,method,method == :local_global,time_i,45),2),4,0) for method in [:local_global, :local_local, :full]]," | ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," | ")

# function error_line(time_i::String)
#   join(
#     [
#       join(
#         [
#           rpad(
#             round(
#               to_score_2_per_model(diff_values,global_value),
#               2
#             ),
#             4,
#             0
#           )
#           for method in [:local_global, :local_local, :full]
#         ],
#         " | "
#       )
#       for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]
#     ],
#     " | "
#   )
# end
error_line(time_i,averaging,grid) = join([join([rpad(round(to_score_2_per_model("mmm","abrupt4x",flux,averaging,grid,method,false,time_i),digits=2),4,"0") for method in [:local_global, :local_local, :full]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"
error_line(averaging,grid) = join([join([rpad(round(to_score_2_per_model_change("mmm","abrupt4x",flux,averaging,grid,method,false),digits=2),4,"0") for method in [:local_global, :local_local, :full]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"
special_error_line(time_i,averaging,grid) = join([join([rpad(round(to_special_score_2_per_model("mmm","abrupt4x",flux,averaging,grid,method,false,time_i),digits=2),4,"0") for method in [:local_global, :local_local, :full]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"
special_error_line(averaging,grid) = join([join([rpad(round(to_special_score_2_per_model_change("mmm","abrupt4x",flux,averaging,grid,method,false),digits=2),4,"0") for method in [:local_global, :local_local, :full]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"

averaged_error_line(time_i,grid) = join([join([rpad(round(to_score_2_per_model("mmm","abrupt4x",flux,averaging,grid,method,false,time_i),digits=2),4,"0") for (method,averaging) in [(:full,"seasonal"),(:local_global,"seasonal"),(:local_local,"seasonal")]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"
averaged_error_change_line(grid) = join([join([rpad(round(to_score_2_per_model_change("mmm","abrupt4x",flux,averaging,grid,method,false),digits=2),4,"0") for (method,averaging) in [(:full,"seasonal"),(:local_global,"seasonal"),(:local_local,"seasonal")]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"
averaged_special_error_line(time_i,grid) = join([join([rpad(round(to_special_score_2_per_model("mmm","abrupt4x",flux,averaging,grid,method,false,time_i),digits=2),4,"0") for (method,averaging) in [(:full,"seasonal"),(:local_global,"seasonal"),(:local_local,"seasonal")]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"
averaged_special_error_line(grid) = join([join([rpad(round(to_special_score_2_per_model_change("mmm","abrupt4x",flux,averaging,grid,method,false),digits=2),4,"0") for (method,averaging) in [(:full,"seasonal"),(:local_global,"seasonal"),(:local_local,"seasonal")]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"


averaged_error_line(time_i,grid) = join([join([rpad(round(to_score_2_per_model(model,"abrupt4x","N",averaging,grid,method,false,time_i),digits=2),4,"0") for (method,averaging) in [(:full,"seasonal"),(:local_global,"seasonal"),(:local_local,"seasonal")]]," & ") for model in ["CESM104","CNRMCM61","GISSE2R","HadCM3L","IPSLCM5A","MPIESM12"]]," & ") * " \\\\"
averaged_error_change_line(grid) = join([join([rpad(round(to_score_2_per_model_change(model,"abrupt4x","N",averaging,grid,method,false),digits=2),4,"0") for (method,averaging) in [(:full,"seasonal"),(:local_global,"seasonal"),(:local_local,"seasonal")]]," & ") for model in ["CESM104","CNRMCM61","GISSE2R","HadCM3L","IPSLCM5A","MPIESM12"]]," & ") * " \\\\"
averaged_special_error_line(time_i,grid) = join([join([rpad(round(to_special_score_2_per_model(model,"abrupt4x","N",averaging,grid,method,false,time_i),digits=2),4,"0") for (method,averaging) in [(:full,"seasonal"),(:local_global,"seasonal"),(:local_local,"seasonal")]]," & ") for model in ["CESM104","CNRMCM61","GISSE2R","HadCM3L","IPSLCM5A","MPIESM12"]]," & ") * " \\\\"
averaged_special_error_change_line(grid) = join([join([rpad(round(to_special_score_2_per_model_change(model,"abrupt4x","N",averaging,grid,method,false),digits=2),4,"0") for (method,averaging) in [(:full,"seasonal"),(:local_global,"seasonal"),(:local_local,"seasonal")]]," & ") for model in ["CESM104","CNRMCM61","GISSE2R","HadCM3L","IPSLCM5A","MPIESM12"]]," & ") * " \\\\"


flux_marker_dict = Dict(
  "N" => "o",
  "LW_clear" => "D",
  "SW_clear" => "s",
  "LW_cloud" => "^",
  "SW_cloud" => "v"
)

flux_color_dict = Dict(
  "N" => "black",
  "LW_clear" => "purple",
  "SW_clear" => "blue",
  "LW_cloud" => "green",
  "SW_cloud" => "red"
)


flux_label_dict = Dict(
  "N" => "net",
  "LW_clear" => "LW clear",
  "LW_cloud" => "LW cloud",
  "SW_clear" => "SW clear",
  "SW_cloud" => "SW cloud"
)

color_dict = Dict(
  :true => "black",
  :global => "blue",
  :local => "orange",
  :full => "green"
)

short_color_dict = Dict(
  :true => "k",
  :global => "b",
  :local => "r",
  :full => "g"
)

period_dict = Dict(
  :early => "2 to 20",
  :late => "21 to end",
  :change => "change"
)

period_label_dict = Dict(
  :early => "years 2-20",
  :late => "years 21-end",
  :change => "change"
)

method_dict = Dict(
  :true => "true",
  :global => "global",
  :local => "local",
  :full => "MR",
  :full_tropics => "MR (tropics)"
)

averaging_dict = Dict(
  :monthly => ":",
  :seasonal => "--",
  :annual => "-."
)

# averaging_method = [("annual",:global),("annual",:local),("monthly",:full)]
# averaging_method = [("annual",:global),("annual",:local),("annual",:full)]
# averaging_method = [("monthly",:global),("monthly",:local),("monthly",:full)]
# averaging_method = [("month_agnostic",:global),("month_agnostic",:local),("month_agnostic",:full)]
averaging_method = [("seasonal",:full),("seasonal",:global),("seasonal",:local)]
default_grid = "12x24_equaldist" # "16x32_equaldist"
averaging_types = ["annual","monthly","month_agnostic","seasonal"]
# make Gregory plots

global_color = "blue"

# control gregories
if false
  figure("control gregories",figsize=(20,12))
  for (gregory_subplot_i,model_name) in enumerate(model_names_abrupt4x)
    subplot(2,3,gregory_subplot_i)
    plot([-.75,.75],[0,0],"k-")

    true_dset = to_dset(to_anomalized_fn(model_name,"control","12x24_equaldist"))
    greg_dict = to_gregory_dict(true_dset,"N")
    years = length(greg_dict["T_surf"])

    plot(greg_dict["T_surf"][end:-1:2],greg_dict["flux"][end:-1:2],"ko",markersize=6)

    if gregory_subplot_i == 6
      # plot(greg_dict["T_surf"][1],greg_dict["flux"][1],color="#FF3DFF","o",markersize=12)
      plot(greg_dict["T_surf"][1],greg_dict["flux"][1],color="#CCCCCC","o",markersize=12)
    end

    plot(greg_dict["T_surf"][1],greg_dict["flux"][1],"ko",markersize=6)

    # F,l = polyfit_exps(greg_dict["T_surf"],greg_dict["flux"],1)
    # plot(extrema(greg_dict["T_surf"]),F+l*collect(extrema(greg_dict["T_surf"])),color=global_color,"-",linewidth=2)

    title(model_name * " ($(years) years)",fontsize=28)
    xlim(-.75,.75)
    ylim(-1.59375,1.59375)

    if gregory_subplot_i > 3
      xlabel("\$\\bar{\\mathrm{T}}\$ ' (K)",fontsize=22)
    else
      xticks([],[])
    end

    if gregory_subplot_i in [1,4]
      ylabel("\$\\bar{\\mathrm{N}}\$ (Wm\$^{-2}\$)",fontsize=22)
    else
      yticks([],[])
    end
  end
  # suptitle("control",fontsize=22)
  subplots_adjust(left=0.06,right=0.98,top=0.91,bottom=0.07,hspace=0.2)
  savefig("control_gregories.eps")
end

# abrupt4x gregories & feedback scatters
if false
  areas = [:all,:special]
  var_names = ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
  function nested_dict()
    dict = Dict()
    for area in areas
      dict[area] = Dict()
      for flux in var_names
        dict[area][flux] = []
      end
    end
    dict
  end

  function averaged_dict()
    dict = Dict()
    for method in [:global,:local,:full]
      dict[method] = Dict()
      for averaging in averaging_types
        dict[method][averaging] = nested_dict()
      end
    end
    dict
  end

  true_flux_dict = nested_dict()
  flux_dict = averaged_dict()
  flux_diff_dict = averaged_dict()

  true_forcing_dict = Dict()
  forcing_dict = Dict()
  for period in [:early,:late]
    true_forcing_dict[period] = nested_dict()
    forcing_dict[period] = averaged_dict()
  end

  true_feedback_dict = Dict()
  feedback_dict = Dict()
  for period in [:early,:late,:change]
    true_feedback_dict[period] = nested_dict()
    feedback_dict[period] = averaged_dict()
  end

  for area in areas
    figure("$(area) feedbacks",figsize=(20,12))
    for i in 331:339
      subplot(i)
      annotate(string("abcdefghi"[i-330]), xy=(0.02, 0.91), xycoords="axes fraction")
      plot([-5,5],[-5,5],"k--")
    end

    for flux in var_names
      gregory_subplot_i = 0
      for model_name in model_names_abrupt4x
        grid = default_grid # model_name == "GISSE2R" ? "24x48_equaldist" : default_grid
        markersize = flux in ["N","LW_cloud","SW_cloud"] ? 6 : 4

        figure("$(area) abrupt4x gregories $(flux)",figsize=(20,12))
        # grid = default_grid
        gregory_subplot_i += 1
        subplot(2,3,gregory_subplot_i)

        true_dset = to_dset(to_anomalized_fn(model_name,"abrupt4x","12x24_equaldist"))
        greg_dict = to_gregory_dict(true_dset,"N")
        years = Int(floor(length(greg_dict["T_surf"])/100)*100)
        true_list_T = greg_dict["T_surf"][1:years]
        greg_dict["T_surf"] = fancy_average(true_list_T)

        true_flux_dict[area][flux] = fancy_average(to_annual_average(to_array(to_dset(to_anomalized_fn(model_name,"abrupt4x",grid)),flux,"$(area == :special ? "special_" : "")global"))[1:years])
        F_true_early, l_true_early = polyfit_exps(greg_dict["T_surf"][2:20],true_flux_dict[area][flux][2:20],1)
        F_true_late, l_true_late = polyfit_exps(greg_dict["T_surf"][21:end],true_flux_dict[area][flux][21:end],1)
        push!(true_forcing_dict[:early][area][flux],F_true_early)
        push!(true_forcing_dict[:late][area][flux],F_true_late)
        push!(true_feedback_dict[:early][area][flux],l_true_early)
        push!(true_feedback_dict[:late][area][flux],l_true_late)
        push!(true_feedback_dict[:change][area][flux],l_true_late-l_true_early)

        for averaging in averaging_types, method in [:full,:global,:local] #[:global,:local,:full]
          # println(to_estimated_fn(model_name,"abrupt4x",(method != :full ? Symbol("local_$(method)") : method),grid))
          flux_dict[method][averaging][area][flux] = fancy_average(to_array(to_dset(to_estimated_fn(model_name,"abrupt4x",(method != :full ? Symbol("local_$(method)") : method),grid)),"$(flux)_$(averaging)","$(area == :special ? "special_" : "")global")[1:years])
          push!(flux_diff_dict[method][averaging][area][flux],true_flux_dict[area][flux][2]-flux_dict[method][averaging][area][flux][2])
          F_early, l_early = polyfit_exps(greg_dict["T_surf"][2:20],flux_diff_dict[method][averaging][area][flux][end].+flux_dict[method][averaging][area][flux][2:20],1)
          F_late, l_late = polyfit_exps(greg_dict["T_surf"][21:end],flux_diff_dict[method][averaging][area][flux][end].+flux_dict[method][averaging][area][flux][21:end],1)
          push!(forcing_dict[:early][method][averaging][area][flux],F_early)
          push!(forcing_dict[:late][method][averaging][area][flux],F_late)
          push!(feedback_dict[:early][method][averaging][area][flux],l_early)
          push!(feedback_dict[:late][method][averaging][area][flux],l_late)
          push!(feedback_dict[:change][method][averaging][area][flux],l_late-l_early)
        end

        plot([0,8],[0,0],"k-")
        plot(greg_dict["T_surf"],true_flux_dict[area][flux],markersize = markersize+2,"o",markeredgecolor="none",color="black",label="true")
        plot(extrema(greg_dict["T_surf"][2:20]),true_forcing_dict[:early][area][flux][end].+true_feedback_dict[:early][area][flux][end].*collect(extrema(greg_dict["T_surf"][2:20])),linewidth=2,"k-")
        plot(extrema(greg_dict["T_surf"][21:end]),true_forcing_dict[:late][area][flux][end].+true_feedback_dict[:late][area][flux][end].*collect(extrema(greg_dict["T_surf"][21:end])),linewidth=2,"k-")

        for (averaging,method) in averaging_method
          plot(greg_dict["T_surf"],flux_diff_dict[method][averaging][area][flux][end].+flux_dict[method][averaging][area][flux],markersize = markersize,"o",markeredgecolor="none",color=color_dict[method],label="$(method_dict[method]) est.")
          plot(extrema(greg_dict["T_surf"][2:20]),forcing_dict[:early][method][averaging][area][flux][end].+feedback_dict[:early][method][averaging][area][flux][end].*collect(extrema(greg_dict["T_surf"][2:20])),color=color_dict[method],linewidth=2,"-")
          plot(extrema(greg_dict["T_surf"][21:end]),forcing_dict[:late][method][averaging][area][flux][end].+feedback_dict[:late][method][averaging][area][flux][end].*collect(extrema(greg_dict["T_surf"][21:end])),color=color_dict[method],linewidth=2,"-")
        end

        # feedback_dict_true = to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales))
        # feedback_dict_global = to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid,estimate_type=:local_global))
        # feedback_dict_local = to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid,estimate_type=:local_local))
        # feedback_dict_full = to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,grid,estimate_type=:global_local_multiple))
        xlim(0,8)

        if flux in ["N","SW_clear"]
          ylim(-1,16)
        elseif flux == "LW_clear"
          ylim(-5,12)
        elseif flux == "LW_cloud"
          ylim(-3,14)
        elseif flux == "SW_cloud"
          ylim(-7,10)
        end


        if gregory_subplot_i in [1,4]
          ylabel("\$\\bar{\\mathrm{N}}\$ (Wm\$^{-2}\$)",fontsize=18)
        end
        if gregory_subplot_i > 3
          xlabel("\$\\bar{\\mathrm{T}}\$ ' (K)",fontsize=18)
        end
        if gregory_subplot_i == 1
          legend(loc="upper left",numpoints=1)
        end
        title(model_name * " ($(years) years)")

        figure("$(area) feedbacks")

        feedback_subplot_i = 330
        for period in [:early,:late,:change], (averaging,method) in averaging_method
          feedback_subplot_i += 1
          subplot(feedback_subplot_i)
          if period == :early
            title(method_dict[method])
          end

          true_unct_fn = to_uncertainty_fn(model_name,period_dict[period],flux,1000,:true,grid,area)
          unct_fn = to_uncertainty_fn(model_name,period_dict[period],flux,1000,method,grid,area)
          # plot([load(true_unct_fn)[dir] for dir in ["low","high"]],feedback_dict[period][method][averaging][area][flux][end]*[1,1],color=flux_color_dict[flux],linestyle="-",zorder=(flux == "N" ? 5 : 4))
          # if isfile(unct_fn)
          #   plot(true_feedback_dict[period][area][flux][end]*[1,1],[load(unct_fn)[dir] for dir in ["low","high"]],color=flux_color_dict[flux],linestyle="-",zorder=(flux == "N" ? 5 : 4))
          # end

          if (model_name == "CESM104") && (feedback_subplot_i == 336)
            plot(true_feedback_dict[period][area][flux][end],feedback_dict[period][method][averaging][area][flux][end],color=flux_color_dict[flux],marker=flux_marker_dict[flux],markersize=markersize+2,markeredgecolor="none",label=flux_label_dict[flux],zorder=(flux == "N" ? 5 : 4))
          else
            plot(true_feedback_dict[period][area][flux][end],feedback_dict[period][method][averaging][area][flux][end],color=flux_color_dict[flux],marker=flux_marker_dict[flux],markersize=markersize+2,markeredgecolor="none",zorder=(flux == "N" ? 5 : 4))
          end

          if method == :full
            if period == :early
              ylabel("\$\\lambda_{4x,early}\$ (years 2-20)\nest. feedback (Wm\$^{-2}\$K\$^{-1}\$)")
            elseif period == :late
              ylabel("\$\\lambda_{4x,late}\$ (years 21-end)\nest. feedback (Wm\$^{-2}\$K\$^{-1}\$)")
            else
              ylabel("\$\\Delta \\lambda_{4x}\$ (change)\nest. feedback (Wm\$^{-2}\$K\$^{-1}\$)")
            end
          end
          if period == :change
            xlabel("true feedback (Wm\$^{-2}\$K\$^{-1}\$)")
          end
        end
      end

      figure("$(area) abrupt4x gregories $(flux)")
      subplots_adjust(left=0.05,right=0.98)
      if area == :all
        # suptitle("abrupt4x",fontsize=22)
        savefig("../figs/gregory/abrupt4x_revised_$(flux)_$(default_grid).eps")
      else
        # suptitle("abrupt4x (30ºS to 90ºN)",fontsize=22)
        savefig("../figs/gregory/special_abrupt4x_revised_$(flux)_$(default_grid).eps")
      end
    end

    figure("$(area) feedbacks")
    for i in 331:336
      subplot(i)
      # if area == :all
        xlim(-2.6,1.5)
        ylim(-2.6,3)
      # else
        # xlim(-3,1.5)
        # ylim(-3,3)

        # xlim(-2,0)
        # ylim(-2,1)

      # end
    end

    for i in 337:339
      subplot(i)
      xlim(-.5,1.5)
      ylim(-.5,1.)
    end

    subplot(336)
    legend(loc="upper right",fontsize=16,bbox_to_anchor=(1.4,0.85),frameon=false,numpoints=1)
    subplots_adjust(left=0.07,bottom=0.06)

    if area == :all
      # suptitle("abrupt4x",fontsize=22)
      savefig("../figs/scatter_plots/feedbacks/revised_$(default_grid).eps")
    else
      # suptitle("abrupt4x (excluding 30ºS and south)",fontsize=22)
      savefig("../figs/scatter_plots/feedbacks/revised_special_$(default_grid).eps")
    end
  end

  println("& " * join(["\\multicolumn{3}{c|}{\\emph{$(flux)}}" for flux in ["net","LW clear","SW clear","LW cloud","SW cloud"]]," & ") * " \\\\")
  println("             & global & local & MR & global & local & MR & global & local & MR & global & local & MR & global & local & MR \\\\")
  println("global error:")
  # for area in [:all,:special], period in [:early,:late,:change], averaging_type in ["annual","month_agnostic","monthly"]
  #   println("~~$(period) ($(averaging_type[[1,end]])): & " * rms_line(feedback_dict,true_feedback_dict,area,period,averaging_type))
  # end
  for area in [:all,:special], period in [:early,:late,:change]
    println("~~$(period) & " * rms_line(feedback_dict,true_feedback_dict,area,period))
  end
end

# make spatial error figures
if false
  println("& " * join(["\\multicolumn{3}{c|}{\\emph{$(flux)}}" for flux in ["net","LW clear","SW clear","LW cloud","SW cloud"]]," & ") * " \\\\")
  println("             & global & local & MR & global & local & MR & global & local & MR & global & local & MR & global & local & MR \\\\")
  println("spatial error:")
  # for averaging_type in ["annual","month_agnostic","monthly"]
  #   println("~~early ($(averaging_type[[1,end]])): & " * error_line(1,averaging_type,default_grid))
  # end
  # for averaging_type in ["annual","month_agnostic","monthly"]
  #   println("~~late ($(averaging_type[[1,end]])): & " * error_line(2,averaging_type,default_grid))
  # end
  # for averaging_type in ["annual","month_agnostic","monthly"]
  #   println("~~change ($(averaging_type[[1,end]])): & " * error_line(averaging_type,default_grid))
  # end
  println("~~early: & " * averaged_error_line(1,default_grid))
  println("~~late: & " * averaged_error_line(2,default_grid))
  println("~~change: & " * averaged_error_line(default_grid))

  println("special spatial error:")
  # for averaging_type in ["annual","month_agnostic","monthly"]
  #   println("~~early ($(averaging_type[[1,end]])): & " * special_error_line(1,averaging_type,default_grid))
  # end
  # for averaging_type in ["annual","month_agnostic","monthly"]
  #   println("~~late ($(averaging_type[[1,end]])): & " * special_error_line(2,averaging_type,default_grid))
  # end
  # for averaging_type in ["annual","month_agnostic","monthly"]
  #   println("~~change ($(averaging_type[[1,end]])): & " * special_error_line(averaging_type,default_grid))
  # end
  println("~~early: & " * averaged_special_error_line(1,default_grid))
  println("~~late: & " * averaged_special_error_line(2,default_grid))
  println("~~change: & " * averaged_special_error_line(default_grid))

  function to_special_spatial_average(vals,lats)
    weights = sin.(vcat([-30],(lats[2:end]+lats[1:end-1])/2,[90])*π/180)
    weights = (weights[2:end] - weights[1:end-1])/1.5
    if length(size(vals)) == 2
      squeeze(mean(vals,length(size(vals))),length(size(vals)))' * weights
    else
      squeeze(mean(vals,length(size(vals))),length(size(vals))) * weights
    end
  end

  function spatial_var(fn,time_i,flux,averaging="")
    dset = to_dset(fn)
    if flux == "net"
      sum([dset[tmp_flux * averaging * "_$(time_i)"][:values][:,:] for tmp_flux in ["LW_clear","SW_clear","LW_cloud","SW_cloud"]])
    else
      dset[flux * averaging * "_$(time_i)"][:values][:,:]
    end
  end
  function to_R2(xs,ys)
    F,l = polyfit_exps(xs,ys,1)
    fs = F+l*xs
    round(1 - sum((ys-fs).^2)/(sum((ys - mean(ys)).^2)),2)
  end

  figure("spatial correlation",figsize=(20,12))

  xlim_dict = Dict(
    "net" => [-4,3],
    "LW_clear" => [-4,1.5],
    "SW_clear" => [-0.2,1.0],
    "LW_cloud" => [-6,5],
    "SW_cloud" => [-6,6]
  )
  ylim_dict = Dict(
    "net" => [-12,12],
    "LW_clear" => [-12,6],
    "SW_clear" => [-2,2],
    "LW_cloud" => [-20,30],
    "SW_cloud" => [-25,20]
  )

  flux_i = 0
  for flux in ["net","LW_clear","SW_clear","LW_cloud","SW_cloud"]
    flux_i += 1

    subplot(3,5,flux_i)
    title(replace(flux,"_"," "))
    true_list_1 = reshape(spatial_var(to_finite_difference_fn("mmm","abrupt4x",default_grid),1,flux),12*24);
    true_ds = to_dset(to_finite_difference_fn("mmm","abrupt4x",default_grid))

    global_list_1 = reshape(to_dset("netcdfs/12x24_equaldist/feedbacks/mmm_control_local_global.nc")["$(flux == "net" ? "N" : flux)_month_agnostic_feedback"][:values],12*24);
    local_list_1 = reshape(spatial_var(to_finite_difference_fn("mmm","abrupt4x",default_grid,:local_local),1,flux,"_annual")/true_ds["T_surf_1_global"][:values][1],12*24);
    full_list_1 = reshape(spatial_var(to_finite_difference_fn("mmm","abrupt4x",default_grid,:full),1,flux,"_monthly")/true_ds["T_surf_1_global"][:values][1],12*24);
    plot(true_list_1,global_list_1,color="blue","o",label=string(to_R2(true_list_1,global_list_1)))
    plot(true_list_1,local_list_1,color="red","o",label=string(to_R2(true_list_1,local_list_1)))
    plot(true_list_1,full_list_1,color="green","o",label=string(to_R2(true_list_1,full_list_1)))
    xls = xlim()
    plot(xls,xls,"k--")
    xlim(xls)
    # xlim(xlim_dict[flux])
    # ylim(ylim_dict[flux])
    # if flux == "LW_clear"
    #   xticks(-3:1.5:1.5)
    # end
    legend(fontsize=12,numpoints=1,loc="lower right")
    if flux == "net"
      ylabel("early\nest. normalized\nflux change (Wm\$^{-2}\$K\$^{-1}\$)")
    end

    subplot(3,5,5+flux_i)
    true_list_2 = reshape(spatial_var(to_finite_difference_fn("mmm","abrupt4x",default_grid),2,flux),12*24);
    global_list_2 = reshape(to_dset("netcdfs/12x24_equaldist/feedbacks/mmm_control_local_global.nc")["$(flux == "net" ? "N" : flux)_month_agnostic_feedback"][:values],12*24);
    local_list_2 = reshape(spatial_var(to_finite_difference_fn("mmm","abrupt4x",default_grid,:local_local),2,flux,"_annual")/true_ds["T_surf_2_global"][:values][1],12*24);
    full_list_2 = reshape(spatial_var(to_finite_difference_fn("mmm","abrupt4x",default_grid,:full),2,flux,"_monthly")/true_ds["T_surf_2_global"][:values][1],12*24);
    plot(true_list_2,global_list_2,color="blue","o",label=string(to_R2(true_list_2,global_list_2)))
    plot(true_list_2,local_list_2,color="red","o",label=string(to_R2(true_list_2,local_list_2)))
    plot(true_list_2,full_list_2,color="green","o",label=string(to_R2(true_list_2,full_list_2)))
    xls = xlim()
    plot(xls,xls,"k--")
    xlim(xls)
    # xlim(xlim_dict[flux])
    # ylim(ylim_dict[flux])
    # if flux == "LW_clear"
    #   xticks(-3:1.5:1.5)
    # end
    legend(fontsize=12,numpoints=1,loc="lower right")
    if flux == "net"
      ylabel("late\nest. normalized\nflux change (Wm\$^{-2}\$K\$^{-1}\$)")
    end

    subplot(3,5,10+flux_i)
    plot(true_list_2-true_list_1,global_list_2-global_list_1,color="blue","o",label=string(to_R2(true_list_2-true_list_1,global_list_2-global_list_1)))
    plot(true_list_2-true_list_1,local_list_2-local_list_1,color="red","o",label=string(to_R2(true_list_2-true_list_1,local_list_2-local_list_1)))
    plot(true_list_2-true_list_1,full_list_2-full_list_1,color="green","o",label=string(to_R2(true_list_2-true_list_1,full_list_2-full_list_1)))
    xls = xlim()
    plot(xls,xls,"k--")
    xlim(xls)
    # xlim(xlim_dict[flux])
    # ylim(ylim_dict[flux])
    # if flux == "LW_clear"
    #   xticks(-3:1.5:1.5)
    # end
    legend(fontsize=12,numpoints=1,loc="lower right")
    xlabel("true normalized\nflux change (Wm\$^{-2}\$K\$^{-1}\$)")
    if flux == "net"
      ylabel("change\nest. normalized\nflux change (Wm\$^{-2}\$K\$^{-1}\$)")
    end
  end
  suptitle("abrupt4x normalized flux change",fontsize=22)
  subplots_adjust(left=0.1,right=0.97)
  savefig("../figs/scatter_plots/spatial_$(default_grid).eps")
end

# years used
if false
  figure("errors",(20,12))

  var_name_dict = Dict(
    "N" => "net",
    "SW_clear" => "SW clear",
    "LW_clear" => "LW clear",
    "SW_cloud" => "SW cloud",
    "LW_cloud" => "LW cloud"
  )

  flux_dict = Dict(
    "N" => "ko",
    "LW_clear" => "rv",
    "SW_clear" => "bv",
    "LW_cloud" => "rD",
    "SW_cloud" => "bD"
  )

  period_dict = Dict(
    :early => "early (years 2-20)",
    :late => "late (years 21-end)",
    :change => "change"
  )

  yus = [350,500,750,1000] #,:all]
  subplot_i = 0

  for score_i in 1:1 # 1:2
    println("score ",score_i)
    run_i = 0
    for period in [:early] # [:early,:late,:change]
    println("period ",period)
      subplot_i += 1
      subplot(230 + subplot_i)
      if subplot_i < 4
        title(period_dict[period])
      end
      run_i += 1

      # suptitle("score $(score_i), years $(time_scales_dict["abrupt4x"][run_i])")
      i = 0
      for flux in ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
        i += 1
        # title(run_i == 1 ? "early period" : "late period")
        if score_i == 1
          plot(yus,[to_special_score_1(flux,period,years_used) for years_used in yus],"$(flux_dict[flux])-",label=var_name_dict[flux])
        else
          if period in [:early,:late]
            # (model_name,run_type,var_name,averaging_type,grid,feedback_type,regridded_estimate,i,years_used=:all,adjacence=:all,tropics_width=:all)
            plot(yus,[to_score_2_per_model("mmm","abrupt4x",flux,"monthly","12x24_equaldist","full",false,period == :early ? 1 : 2,years_used) for years_used in yus],"$(flux_dict[flux])-",label=var_name_dict[flux])
          else
            plot(yus,[to_score_2_per_model_change("mmm","abrupt4x",flux,"monthly","12x24_equaldist","full",false,years_used) for years_used in yus],"$(flux_dict[flux])-",label=var_name_dict[flux])
          end
        end
      end
      # if score_i == 1
      #   ylim(0,0.6)
      # elseif score_i == 2
      #   ylim(0,7)
      # end
      xscale("log")
      if subplot_i > 3
        xlabel("years used")
      end
      if subplot_i in [1,4]
        ylabel("$(score_i == 1 ? "feedback error (30ºS to 90ºN)" : "spatial error") (\$Wm^{-2}K^{-1}\$)") #,fontsize=12)
      end
      xlim(extrema(yus)...)
      xticks(yus,yus)
      if score_i == 1
        ylim(0,0.7)
      else
        ylim(0,2.5)
      end
      if subplot_i == 3
        legend(loc="center right", bbox_to_anchor=(1.5,0.0),numpoints=1,frameon=false,fontsize=16)
      end
    end
  end

  subplots_adjust(left=0.05,right=0.87,bottom=0.07,wspace=0.25,top=0.93)
  savefig("../figs/years_used/all.eps")

end

# state dependence (just late)
if false
  multi_data = OrderedDict()
  model_name = "MPIESM12"
  runs = ["abrupt2x","abrupt4x","abrupt8x","abrupt16x"] # ,"abrupt32x"]
  for run_type in runs
    multi_data[run_type] = OrderedDict(
      "true" => to_feedback_dict(to_feedback_fn(model_name,run_type,:global_global_time_scales,"12x24_equaldist")),
      "full" => to_feedback_dict(to_feedback_fn(model_name,run_type,:global_global_time_scales,"12x24_equaldist",estimate_type=:full))
    )
  end

  flux_dict = Dict(
    "N" => "ko",
    "LW_clear" => "rv",
    "SW_clear" => "bv",
    "LW_cloud" => "rD",
    "SW_cloud" => "bD"
  )

  figure(figsize=(9,6))
  markersize = 7
  fontsize = 10
  offset = 0.02
  title("late period (years 21-end)")
  for flux in ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]

    # fill_between([-2.5,1.5],-.3+[-2.5,1.5],.3+[-2.5,1.5],color="#DDDDDD")
    true_late = [multi_data[run_type]["true"][flux]["21 to end"]["special_annual"]["feedback"] for run_type in runs]
    est_late = [multi_data[run_type]["full"][flux]["21 to end"]["special_monthly"]["feedback"] for run_type in runs]
    plot(true_late,est_late,flux_dict[flux],markersize=markersize,label=replace(replace(flux,"_"," "),"N","net"))
    for run_type in runs
      text(offset+multi_data[run_type]["true"][flux]["21 to end"]["special_annual"]["feedback"],offset+multi_data[run_type]["full"][flux]["21 to end"]["special_monthly"]["feedback"],replace(replace(run_type,"x",""),"abrupt",""),fontsize=fontsize)
    end

    if flux == "N"
      plot([-2.5,1.5],[-2.5,1.5],"k--")
    end

    # for model_name in model_names_abrupt4x
    #   unct_dict = load(to_uncertainty_fn(model_name,"51 to end",flux,10000))
    #   x = data[model_name]["true"][flux]["51 to end"]["annual"]["feedback"]
    #   x_unct = data[model_name]["true"][flux]["51 to end"]["annual"]["feedback_unct"]
    #   y = data[model_name]["full"][flux]["51 to end"]["monthly"]["feedback"]
    #   plot([x,x],[unct_dict["low"],unct_dict["high"]],color=flux_dict[flux][1:1])
    #   plot([x-x_unct[1],x+x_unct[1]],[y,y],color=flux_dict[flux][1:1])
    # end
    xlim(-2.5,1.5)
    ylim(-2.5,1.5)
    ylabel("estimated feedback (\$Wm^{-2}K^{-1}\$)")
    xlabel("true feedback (\$Wm^{-2}K^{-1}\$)")
  end

  legend(loc="center right", bbox_to_anchor=(1.28,0.51),numpoints=1,frameon=false,fontsize=14)
  subplots_adjust(left = 0.12,right=0.8,bottom=0.12)
  savefig("../figs/scatter_plots/multi_revised.eps")
end

# state dependence (with early and late)
if false
  multi_data = OrderedDict()
  model_name = "MPIESM12"
  runs = ["abrupt2x","abrupt4x","abrupt8x","abrupt16x"] # ,"abrupt32x"]
  for run_type in runs
    multi_data[run_type] = OrderedDict(
      "true" => to_feedback_dict(to_feedback_fn(model_name,run_type,:global_global_time_scales,"12x24_equaldist")),
      "full" => to_feedback_dict(to_feedback_fn(model_name,run_type,:global_global_time_scales,"12x24_equaldist",estimate_type=:full))
    )
  end

  flux_dict = Dict(
    "N" => "ko",
    "LW_clear" => "rv",
    "SW_clear" => "bv",
    "LW_cloud" => "rD",
    "SW_cloud" => "bD"
  )

  figure(figsize=(20,6))
  markersize = 7
  fontsize = 10
  offset = 0.02
  for flux in ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
    subplot(131)
    title("years 2-20")
    # fill_between([-2.5,1.5],-.3+[-2.5,1.5],.3+[-2.5,1.5],color="#DDDDDD")
    true_early = [multi_data[run_type]["true"][flux]["2 to 20"]["special_annual"]["feedback"] for run_type in runs]
    est_early = [multi_data[run_type]["full"][flux]["2 to 20"]["special_monthly"]["feedback"] for run_type in runs]
    plot(true_early,est_early,flux_dict[flux],markersize=markersize)
    for run_type in runs
      text(offset+multi_data[run_type]["true"][flux]["2 to 20"]["special_annual"]["feedback"],offset+multi_data[run_type]["full"][flux]["2 to 20"]["special_monthly"]["feedback"],replace(replace(run_type,"x",""),"abrupt",""),fontsize=fontsize)
    end


    # for model_name in model_names_abrupt4x
    #   unct_dict = load(to_uncertainty_fn(model_name,"6 to 50",flux,10000))
    #   x = data[model_name]["true"][flux]["6 to 50"]["annual"]["feedback"]
    #   x_unct = data[model_name]["true"][flux]["6 to 50"]["annual"]["feedback_unct"]
    #   y = data[model_name]["full"][flux]["6 to 50"]["monthly"]["feedback"]
    #   plot([x,x],[unct_dict["low"],unct_dict["high"]],color=flux_dict[flux][1:1])
    #   plot([x-x_unct[1],x+x_unct[1]],[y,y],color=flux_dict[flux][1:1])
    #   # text(x,unct_dict["low"],model_name)
    # end
    if flux == "N"
      plot([-2.5,1.5],[-2.5,1.5],"k--")
    end
    xlim(-2.5,1.5)
    ylim(-2.5,1.5)
    ylabel("estimated feedback (\$Wm^{-2}K^{-1}\$)")
    xlabel("true feedback (\$Wm^{-2}K^{-1}\$)")
    subplot(132)
    title("years 21-end")
    # fill_between([-2.5,1.5],-.3+[-2.5,1.5],.3+[-2.5,1.5],color="#DDDDDD")
    true_late = [multi_data[run_type]["true"][flux]["21 to end"]["special_annual"]["feedback"] for run_type in runs]
    est_late = [multi_data[run_type]["full"][flux]["21 to end"]["special_monthly"]["feedback"] for run_type in runs]
    plot(true_late,est_late,flux_dict[flux],markersize=markersize,label=replace(replace(flux,"_"," "),"N","net"))
    for run_type in runs
      text(offset+multi_data[run_type]["true"][flux]["21 to end"]["special_annual"]["feedback"],offset+multi_data[run_type]["full"][flux]["21 to end"]["special_monthly"]["feedback"],replace(replace(run_type,"x",""),"abrupt",""),fontsize=fontsize)
    end

    if flux == "N"
      plot([-2.5,1.5],[-2.5,1.5],"k--")
    end

    # for model_name in model_names_abrupt4x
    #   unct_dict = load(to_uncertainty_fn(model_name,"51 to end",flux,10000))
    #   x = data[model_name]["true"][flux]["51 to end"]["annual"]["feedback"]
    #   x_unct = data[model_name]["true"][flux]["51 to end"]["annual"]["feedback_unct"]
    #   y = data[model_name]["full"][flux]["51 to end"]["monthly"]["feedback"]
    #   plot([x,x],[unct_dict["low"],unct_dict["high"]],color=flux_dict[flux][1:1])
    #   plot([x-x_unct[1],x+x_unct[1]],[y,y],color=flux_dict[flux][1:1])
    # end
    xlim(-2.5,1.5)
    ylim(-2.5,1.5)
    xlabel("true feedback (\$Wm^{-2}K^{-1}\$)")

    subplot(133)
    title("change")
    # fill_between([-2.5,1.5],-.3+[-2.5,1.5],.3+[-2.5,1.5],color="#DDDDDD")
    plot(true_late-true_early,est_late-est_early,flux_dict[flux],markersize=markersize,label=replace(replace(flux,"_"," "),"N","net"))
    for run_type in runs
      text(offset+multi_data[run_type]["true"][flux]["21 to end"]["special_annual"]["feedback"]-multi_data[run_type]["true"][flux]["2 to 20"]["special_annual"]["feedback"],offset+multi_data[run_type]["full"][flux]["21 to end"]["special_monthly"]["feedback"]-multi_data[run_type]["full"][flux]["2 to 20"]["special_monthly"]["feedback"],replace(replace(run_type,"x",""),"abrupt",""),fontsize=fontsize)
    end

    if flux == "N"
      plot([-2.5,1.5],[-2.5,1.5],"k--")
    end

    # for model_name in model_names_abrupt4x
    #   unct_dict = load(to_uncertainty_fn(model_name,"51 to end",flux,10000))
    #   x = data[model_name]["true"][flux]["51 to end"]["annual"]["feedback"]
    #   x_unct = data[model_name]["true"][flux]["51 to end"]["annual"]["feedback_unct"]
    #   y = data[model_name]["full"][flux]["51 to end"]["monthly"]["feedback"]
    #   plot([x,x],[unct_dict["low"],unct_dict["high"]],color=flux_dict[flux][1:1])
    #   plot([x-x_unct[1],x+x_unct[1]],[y,y],color=flux_dict[flux][1:1])
    # end
    xlim(-2.5,1.5)
    ylim(-2.5,1.5)
    xlabel("true feedback (\$Wm^{-2}K^{-1}\$)")
  end


  legend(loc="center right", bbox_to_anchor=(1.4,0.51),numpoints=1,frameon=false,fontsize=14)
  subplots_adjust(left = 0.07,right=0.9,bottom=0.1)
  # savefig("../figs/scatter_plots/ox_multi.eps")
end

# nonfull table
if false
  for area in [:all] #,:special]
    println("area" * (area == :all ? "" : " excluding southern ocean"))
    grid_tmp = "12x24_equaldist" # (model_name) -> (model_name == "GISSE2R") ? "24x48_equaldist" :
    fluxes = ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
    println(join(vcat(["model"],[flux_label_dict[flux] for flux in fluxes],[flux_label_dict[flux] for flux in fluxes],[flux_label_dict[flux] for flux in fluxes])," & ") * " \\\\\n")

    # for period in [:early,:late]
    scales = ["early","late","change"]
    for scale_i in 1:3
      println(scales[scale_i] * " &"^8 * " \\\\")

      for model_name in all_model_names
        print(model_name)
        for method in [:true,:global,:local], flux in ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
          print(" & $(round(feedback_scale(model_name,method,flux,scale_i,"seasonal",area),digits=2))")
        end
        print(" \\\\\n")
      end
      println("\\hline")
      print("mean")
      for method in [:true,:global,:local], flux in ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
        print(" & $(round(mean([feedback_scale(model_name,method,flux,scale_i,"seasonal",area) for model_name in all_model_names]),digits=2))")
      end
      print(" \\\\\n")
      println("\\hline")
      print("mean (no GISSE2R)")
      for method in [:true,:global,:local], flux in ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
        print(" & $(round(mean([feedback_scale(model_name,method,flux,scale_i,"seasonal",area) for model_name in all_model_names[[1,2,4,5,6]]]),digits=2))")
      end
      print(" \\\\\n")
      println("\\hline")
    end
  end
end

# local & nonlocal table
if false
  for area in [:all] #,:special]
    println("area" * (area == :all ? "" : " excluding southern ocean"))
    grid_tmp = "12x24_equaldist" # (model_name) -> (model_name == "GISSE2R") ? "24x48_equaldist" :
    fluxes = ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
    println(join(vcat(["model"],[flux_label_dict[flux] for flux in fluxes],[flux_label_dict[flux] for flux in fluxes],[flux_label_dict[flux] for flux in fluxes])," & ") * " \\\\\n")

    # for period in [:early,:late]
    scales = ["early","late","change"]
    for scale_i in 1:3
      println(scales[scale_i] * " &"^8 * " \\\\")

      for model_name in all_model_names
        print(model_name)
        for method in [:full,:full_local,:full_nonlocal], flux in ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
          print(" & $(round(feedback_scale(model_name,method,flux,scale_i,"seasonal",area),digits=2))")
        end
        print(" \\\\\n")
      end
      println("\\hline")
      print("mean")
      for method in [:full,:full_local,:full_nonlocal], flux in ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
        print(" & $(round(mean([feedback_scale(model_name,method,flux,scale_i,"seasonal",area) for model_name in all_model_names]),digits=2))")
      end
      print(" \\\\\n")
      println("\\hline")
      print("mean (no GISSE2R)")
      for method in [:full,:full_local,:full_nonlocal], flux in ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
        print(" & $(round(mean([feedback_scale(model_name,method,flux,scale_i,"seasonal",area) for model_name in all_model_names[[1,2,4,5,6]]]),digits=2))")
      end
      print(" \\\\\n")
      println("\\hline")
    end
  end
end

# who makes more change
if false
  for model_name in all_model_names
    print(model_name)
    print(" & ")
    a = feedback_scale(model_name,:true,"N",3,"seasonal",:all)
    s = feedback_scale(model_name,:true,"N",3,"seasonal",:special)
    print(round(a,digits=2))
    print(" & ")
    print(round(a-s,digits=2))
    print(" & ")
    print(round(s,digits=2))
    for flux in ["LW_clear","SW_clear","LW_cloud","SW_cloud"]
      print(" & ")
      print(round(feedback_scale(model_name,:true,flux,3,"seasonal",:special),digits=2))
    end
    print(" \\\\\n")
  end

  for model_name in all_model_names
    print(model_name)
    print(" & ")
    a = feedback_scale(model_name,:full,"N",3,"seasonal",:all)
    s = feedback_scale(model_name,:full,"N",3,"seasonal",:special)
    print(round(a,digits=2))
    print(" & ")
    print(round(a-s,digits=2))
    print(" & ")
    print(round(s,digits=2))
    for flux in ["LW_clear","SW_clear","LW_cloud","SW_cloud"]
      print(" & ")
      print(round(feedback_scale(model_name,:full,flux,3,"seasonal",:special),digits=2))
    end
    print(" \\\\\n")
  end
end

# uncertainties
if false
  for method in ["full","full_local","full_nonlocal"]
    fdicts = [to_feedback_dict(to_feedback_fn(model_name,"abrupt4x",:global_global_time_scales,"12x24_equaldist",estimate_type=Symbol(method)))["N"] for model_name in setdiff(all_model_names,["GISSE2R"])]
    for time_period in ["2 to 20","21 to end","change"]

      if time_period == "change"

        println(round(mean([fdict["21 to end"]["special_monthly"]["feedback"][1] - fdict["2 to 20"]["special_monthly"]["feedback"][1] for fdict in fdicts]),2))
      else
        println(round(mean([round(fdict[time_period]["special_monthly"]["feedback"][1],2) for fdict in fdicts]),2))
      end

      # fn = "jlds/uncertainties/mmm_no_GISSE2R_$(time_period)_N_1000_$(method)_12x24_equaldist_special.jld"
      # println((round(load(fn)["low"],2),round(load(fn)["high"],2)))
    end
  end
end

# feedback time series
if false
  using GeoStats
  using LocallyWeightedRegression
  using Interpolations
  var_names = ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
  for area in [:all] #,:special]
    for (subplot_i,model_name) in enumerate(["CESM104","CNRMCM61","GISSE2R","HadCM3L","IPSLCM5A","MPIESM12"]) # model_names_abrupt4x "CCSM3",
      println(model_name)
      dset = to_dset(to_anomalized_fn(model_name,"abrupt4x","12x24_equaldist"))

      monthly_Ts = to_array(dset,"T_surf","$(area == :special ? "special_" : "")global")
      Ts = to_gregory_dict(dset,"N")["T_surf"]
      fa = fancy_average_maker(1000,50)
      fancy_T = fa(Ts)[2:end] # fancy_average(Ts)
      fancy_T = convert(Array{Float64},fancy_T)
      extra = 0.1
      min_T = min(fancy_T...) - extra
      max_T = max(fancy_T...) + extra
      points = 400
      ns = 30
      varss = 20
      T_grid = range(min_T, stop=max_T, length=points)

      ts = 1:length(fancy_T)

      get_T(t) = fancy_T[t]
      get_T_i(T) = 1+((T-min_T)/(max_T-min_T))*(points-1)

      for flux in var_names
        true_flux = 0
        for (method_i,method) in enumerate([:true,:full,:global,:local]) # ,:full_tropics
          figure("$(area) feedback $(flux)",figsize=(20,12))
          subplot(2,3,subplot_i)

          if method == :true
            monthly_flux = to_array(to_dset(to_anomalized_fn(model_name,"abrupt4x","12x24_equaldist")),flux,"$(area == :special ? "special_" : "")global")
            fluxes = to_annual_average(monthly_flux)
            plot([20,20],[-3.0,3.0],"k--",zorder=0.5)
          else
            if method != :full_tropics
              efn = to_estimated_fn(model_name,"abrupt4x",(method != :full ? Symbol("local_$(method)") : method),"12x24_equaldist")
            else
              efn = to_estimated_fn(model_name,"abrupt4x",:full,"12x24_equaldist",tropics=true)
            end
            # fluxes = to_array(to_dset(efn),"$(flux)_$(method != :global ? "monthly" : "annual")","$(area == :special ? "special_" : "")global")
            fluxes = to_array(to_dset(efn),"$(flux)_seasonal","$(area == :special ? "special_" : "")global")
          end
          println(method)

          fancy_flux = fa(fluxes)[2:end]
          fancy_flux = convert(Array{Float64},fancy_flux)
          if method == :true
            true_flux = fancy_flux
          end


          color_line = method == :full_tropics ? color_dict[:full] : color_dict[method]
          average_toggle = false
          if !average_toggle
            # smooth 'it
            solution = solve(EstimationProblem(PointSetData(Dict(:y => fancy_flux), transpose(fancy_T)), RegularGrid((min_T,), (max_T,), dims=(points,)), :y), LocalWeightRegress(:y => (neighbors=ns,variogram=ExponentialVariogram(range=varss/10),))) #()))
            Nhat, Nvar = solution[:y]

            # get feedback
            Nhat_interp = interpolate(Nhat,BSpline(Linear()))
            # interp_linear = LinearInterpolation((x[2:end].+x[1:end-1])./2, (vs[2:end].-vs[1:end-1])./(x[2:end].-x[1:end-1]), extrapolation_bc = Line())
            feedbacks = [(Nhat_interp(get_T_i(get_T(t)+.01))-Nhat_interp(get_T_i(get_T(t))))/.01 for t in ts]

            if subplot_i == 6
              plot(fa(1:length(Ts))[2:end],feedbacks,color=color_line,linestyle=(method == :full_tropics ? "--" : "-"),linewidth=(method in [:true,:full] ? 2 : 1),zorder=(2-.0001*method_i)) #averaging_dict[averaging])
              if subplot_i == 6
                plot([1000,1000],[1000,1001],linestyle="-",color=color_dict[method],marker="o",label=method_dict[method],zorder=(2-.0001*method_i),markersize=6)
              end
            else
              plot(fa(1:length(Ts))[2:end],feedbacks,color=color_line,linewidth=(method in [:true,:full] ? 2 : 1),linestyle=(method == :full_tropics ? "--" : "-"),zorder=(2-.0001*method_i)) #averaging_dict[averaging])
            end
            scale_1 = exp(mean(log.([2,20])))
            scale_2 = exp(mean(log.([21,1000])))
            plot(scale_1,feedback_scale(model_name,method,flux,1),color=color_line,linestyle="none",marker="o",zorder=(2-.0001*method_i),markersize=6)
            plot(scale_2,feedback_scale(model_name,method,flux,2),color=color_line,linestyle="none",marker="o",zorder=(2-.0001*method_i),markersize=6)
          end
          ylim(-3,3)
          xscale("log")
          xlim(2,1000)

          if method != :full_tropics

            figure("$(area) gregories $(flux)",figsize=(20,12))
            subplot(2,3,subplot_i)
            plot(Ts[1:1000],fluxes[1:1000] .+ (true_flux[1] - fancy_flux[1]),linestyle="none",color="gray",marker="o",markersize=1,zorder=0.5)
            method_markersize = Dict(
              # :true => 8,
              :true => 6,
              :global => 4,
              :local => 4,
              :full => 6
            )
            plot(fancy_T,fancy_flux .+ (true_flux[1] - fancy_flux[1]),linestyle="none",color=color_dict[method],marker="o",markersize=method_markersize[method],zorder=(2-.0001*method_i))
            if !average_toggle
              plot(T_grid,Nhat.+ (true_flux[1] - fancy_flux[1]),color=color_line,linestyle="-",zorder=(2-.0001*method_i))
            end
            if subplot_i == 6
              plot([1000,1000],[1000,1001],linestyle="-",color=color_dict[method],marker="o",label=method_dict[method],zorder=(2-.0001*method_i))
            end
            ylim(-1,10)
            xlim(0,9)
            plot([0,9],[0,0],"k-",zorder=0.5)
            # ylim(-1,6)
            # xlim(2,8)
            # plot([2,8],[0,0],"k-",zorder=0.5)
            # F,l_20 = polyfit_exps(Ts[2:20],fluxes[2:20],1)
            # plot([Ts[2],Ts[20]],F.+l_20*[Ts[2],Ts[20]].+ (true_flux[1] - fancy_flux[1]),color="pink",linestyle="--")
            # F,l_end = polyfit_exps(Ts[21:end],fluxes[21:end],1)
            # plot([Ts[21],Ts[end]],F.+l_end*[Ts[21],Ts[end]].+ (true_flux[1] - fancy_flux[1]),color="pink",linestyle="--")
            #
            # figure("$(flux_label_dict[flux]) feedback",figsize=(20,12))
            # subplot(2,3,subplot_i)
            # plot([2,20],[l_20,l_20],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="--") #averaging_dict[averaging])
            # plot([21,length(Ts)],[l_end,l_end],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="--") #averaging_dict[averaging])

            if average_toggle
              figure("$(area) gregories $(flux)",figsize=(20,12))
              F_20,l_20 = polyfit_exps(fancy_T[1:14],fancy_flux[1:14],1)
              plot([fancy_T[1],fancy_T[14]],F_20 .+l_20*[fancy_T[1],fancy_T[14]].+ (true_flux[1] - fancy_flux[1]),color=color_line,linestyle="-",zorder=4)
              F_end,l_end = polyfit_exps(fancy_T[15:end],fancy_flux[15:end],1)
              if subplot_i == 6
                plot([fancy_T[15],fancy_T[end]],F_end .+l_end*[fancy_T[15],fancy_T[end]].+ (true_flux[1] - fancy_flux[1]),color=color_line,linestyle="-",label=method_dict[method],zorder=4)
              else
                plot([fancy_T[15],fancy_T[end]],F_end .+l_end*[fancy_T[15],fancy_T[end]].+ (true_flux[1] - fancy_flux[1]),color=color_line,linestyle="-",zorder=4)
              end

              figure("$(area) feedback $(flux)",figsize=(20,12))
              subplot(2,3,subplot_i)
              plot([2,20],[l_20,l_20],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="-") #averaging_dict[averaging])
              if subplot_i == 6
                plot([21,length(Ts)],[l_end,l_end],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="-",label=method_dict[method]) #averaging_dict[averaging])
              else
                plot([21,length(Ts)],[l_end,l_end],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="-") #averaging_dict[averaging])
              end
            end
            if method == :true
              if average_toggle
                figure("$(area) residual $(flux)",figsize=(20,12))
                subplot(2,3,subplot_i)

                interp_linear_1 = LinearInterpolation([fancy_T[1],fancy_T[14]],F_20.+l_20*[fancy_T[1],fancy_T[14]],extrapolation_bc = Line())
                interp_linear_2 = LinearInterpolation([fancy_flux[15],Ts[end]],F_end .+l_end*[fancy_flux[15],Ts[end]],extrapolation_bc = Line())
                title("piecewise")
                plot([-3,43],[0,0],"k-")
                plot(vcat(fancy_flux[1:14] .- interp_linear_1.(fancy_T[1:14]),fancy_flux[15:end] .- interp_linear_2.(fancy_T[15:end])),"ko")
                ylim(-.7,.7)
              else
                figure("$(area) residual $(flux)",figsize=(20,12))
                subplot(2,3,subplot_i)
                title("loess")
                plot([-3,43],[0,0],"k-")
                plot(fancy_flux .- Nhat_interp.(get_T_i.(fancy_T)),"ko")
                ylim(-.7,.7)
              end
            end
          end
        end

        # figure("$(area) abrupt4x feedback time series $(flux)",figsize=(20,12))
        # plot([3,150],[0,0],"k--")
        # # markersize = flux in ["N","LW_cloud","SW_cloud"] ? 6 : 4
        #
        # subplot(2,3,subplot_i)
        #
        #
        #
        # plot(1:length(true_flux),(true_flux .- true_flux[1])./(Ts .- Ts[1]),"k-")
        #
        # for method in [:global,:local,:full]
          # est_flux = to_array(to_dset(to_estimated_fn(model_name,"abrupt4x",(method != :full ? Symbol("local_$(method)") : method),"12x24_equaldist")),"$(flux)_$(method != :global ? "monthly" : "annual")","$(area == :special ? "special_" : "")global")
        #   plot(1:length(est_flux),((est_flux .- est_flux[1])./(Ts .- Ts[1])),color=color_dict[method],linestyle="-")
        # end

        figure("$(area) feedback $(flux)")
        title(model_name)
        # ylim(-4,4)
        # xlim(2,150)
        if subplot_i in [2,3,5,6]
          yticks([])
        end
        if subplot_i < 4
          xticks([])
        end
        if subplot_i > 3
          xlabel("time (yr)")
        end
        if subplot_i in [1,4]
          ylabel("feedback (\$Wm^{-2}K^{-1}\$)")
        end
        figure("$(area) gregories $(flux)")
        title(model_name)
        # ylim(-4,4)
        # xlim(2,150)
        if subplot_i in [2,3,5,6]
          yticks([])
        end
        if subplot_i < 4
          xticks([])
        end
        if subplot_i > 3
          xlabel("\$\\bar{T}'\$ (K)")
        end
        if subplot_i in [1,4]
          ylabel("\$\\bar{N}\$ (\$Wm^{-2}\$)")
        end
        figure("$(area) residual $(flux)")
      end
    end

    for flux in var_names
      for fig in ["$(area) feedback $(flux)","$(area) gregories $(flux)","$(area) residual $(flux)"]
        figure(fig)
        # suptitle(fig)
        legend(loc="upper right",fontsize=20,bbox_to_anchor=(1.55,1.4),frameon=false,numpoints=1)
        subplots_adjust(left=0.06,hspace=0.15,right=0.85,wspace=0.1)
        savefig("../figs/$(replace(fig," "=>"_")).eps")
      end
    end
  end
end

# cfmip feedback time series
if false
  using GeoStats
  using LocallyWeightedRegression
  using Interpolations
  var_names = ["N","LW_clear","SW_clear","LW_cloud","SW_cloud"]
  for area in [:all] #,:special]
    for (subplot_i,model_name) in enumerate(["CESM104","CNRMCM61","GISSE2R","HadCM3L","IPSLCM5A","MPIESM12"]) # model_names_abrupt4x "CCSM3",
      println(model_name)
      dset = to_dset(to_anomalized_fn(model_name,"abrupt4x","12x24_equaldist"))

      monthly_Ts = to_array(dset,"T_surf","$(area == :special ? "special_" : "")global")
      Ts = to_gregory_dict(dset,"N")["T_surf"]
      fa = fancy_average_maker(1000,50)
      fancy_T = fa(Ts)[2:end] # fancy_average(Ts)
      fancy_T = convert(Array{Float64},fancy_T)
      extra = 0.1
      min_T = min(fancy_T...) - extra
      max_T = max(fancy_T...) + extra
      points = 400
      ns = 30
      varss = 20
      T_grid = range(min_T, stop=max_T, length=points)

      ts = 1:length(fancy_T)

      get_T(t) = fancy_T[t]
      get_T_i(T) = 1+((T-min_T)/(max_T-min_T))*(points-1)

      for flux in var_names
        true_flux = 0
        for (method_i,method) in enumerate([:true,:full]) #,:global,:local]) # ,:full_tropics
          # if flux == "N"
            figure("$(area) feedback $(flux)",figsize=(20,12))
          # else
            # figure("$(area) feedback others",figsize=(20,12))
          # end
          subplot(2,3,subplot_i)

          if method == :true
            monthly_flux = to_array(to_dset(to_anomalized_fn(model_name,"abrupt4x","12x24_equaldist")),flux,"$(area == :special ? "special_" : "")global")
            fluxes = to_annual_average(monthly_flux)
            # plot([20,20],[-3.0,3.0],"k--",zorder=0.5)
          else
            if method != :full_tropics
              efn = to_estimated_fn(model_name,"abrupt4x",(method != :full ? Symbol("local_$(method)") : method),"12x24_equaldist")
            else
              efn = to_estimated_fn(model_name,"abrupt4x",:full,"12x24_equaldist",tropics=true)
            end
            # fluxes = to_array(to_dset(efn),"$(flux)_$(method != :global ? "monthly" : "annual")","$(area == :special ? "special_" : "")global")
            fluxes = to_array(to_dset(efn),"$(flux)_seasonal","$(area == :special ? "special_" : "")global")
          end
          println(method)

          fancy_flux = fa(fluxes)[2:end]
          fancy_flux = convert(Array{Float64},fancy_flux)
          if method == :true
            true_flux = fancy_flux
          end


          color_line = method == :full_tropics ? color_dict[:full] : color_dict[method]
          average_toggle = false
          linestyle_dict = Dict(
            "N" => "-",
            "LW_clear" => "-",
            "SW_clear" => "-",
            "LW_cloud" => "-",
            "SW_cloud" => "-"
          )
          if !average_toggle
            # smooth 'it
            solution = solve(EstimationProblem(PointSetData(Dict(:y => fancy_flux), transpose(fancy_T)), RegularGrid((min_T,), (max_T,), dims=(points,)), :y), LocalWeightRegress(:y => (neighbors=ns,variogram=ExponentialVariogram(range=varss/10),))) #()))
            Nhat, Nvar = solution[:y]

            # get feedback
            Nhat_interp = interpolate(Nhat,BSpline(Linear()))
            # interp_linear = LinearInterpolation((x[2:end].+x[1:end-1])./2, (vs[2:end].-vs[1:end-1])./(x[2:end].-x[1:end-1]), extrapolation_bc = Line())
            feedbacks = [(Nhat_interp(get_T_i(get_T(t)+.01))-Nhat_interp(get_T_i(get_T(t))))/.01 for t in ts]

            if subplot_i == 6
              plot(fa(1:length(Ts))[2:end],feedbacks,color=color_line,linestyle=linestyle_dict[flux],linewidth=(method in [:true] ? 3 : 2),zorder=(2+.0001*method_i)) #averaging_dict[averaging])
              if subplot_i == 6
                plot([1000,1000],[1000,1001],linestyle=linestyle_dict[flux],color=color_dict[method],label=method_dict[method] * " " * flux,zorder=(2+.0001*method_i),markersize=6)
              end
            else
              plot(fa(1:length(Ts))[2:end],feedbacks,color=color_line,linewidth=(method in [:true] ? 3 : 2),linestyle=linestyle_dict[flux],zorder=(2+.0001*method_i)) #averaging_dict[averaging])
            end
            # scale_1 = exp(mean(log.([2,20])))
            # scale_2 = exp(mean(log.([21,1000])))
            # plot(scale_1,feedback_scale(model_name,method,flux,1),color=color_line,linestyle="none",marker="o",zorder=(2-.0001*method_i),markersize=6)
            # plot(scale_2,feedback_scale(model_name,method,flux,2),color=color_line,linestyle="none",marker="o",zorder=(2-.0001*method_i),markersize=6)
          end
          # if flux == "SW_cloud"
          #   ylim(-1.5,1)
          # elseif flux == "LW_coud"
          #   ylim(-.5,2)
          # elseif flux == "SW_clear"
          #   ylim(0,2.5)
          # elseif flux == "LW_clear"
          #   ylim()
          # else
          if flux == "N"
            ylim(-2.5,0)
          else
            ylim(-1.5,1)
          end
          xscale("log")
          xlim(2,1000)

          if method != :full_tropics

            figure("$(area) gregories $(flux)",figsize=(20,12))
            subplot(2,3,subplot_i)
            plot(Ts[1:1000],fluxes[1:1000] .+ (true_flux[1] - fancy_flux[1]),linestyle="none",color="gray",marker="o",markersize=1,zorder=0.5)
            method_markersize = Dict(
              # :true => 8,
              :true => 6,
              :global => 4,
              :local => 4,
              :full => 6
            )
            plot(fancy_T,fancy_flux .+ (true_flux[1] - fancy_flux[1]),linestyle="none",color=color_dict[method],marker="o",markersize=method_markersize[method],zorder=(2-.0001*method_i))
            if !average_toggle
              plot(T_grid,Nhat.+ (true_flux[1] - fancy_flux[1]),color=color_line,linestyle="-",zorder=(2-.0001*method_i))
            end
            if subplot_i == 6
              plot([1000,1000],[1000,1001],linestyle="-",color=color_dict[method],marker="o",label=method_dict[method],zorder=(2-.0001*method_i))
            end
            ylim(-1,7)
            xlim(0,9)
            plot([0,9],[0,0],"k-",zorder=0.5)
            # ylim(-1,6)
            # xlim(2,8)
            # plot([2,8],[0,0],"k-",zorder=0.5)
            # F,l_20 = polyfit_exps(Ts[2:20],fluxes[2:20],1)
            # plot([Ts[2],Ts[20]],F.+l_20*[Ts[2],Ts[20]].+ (true_flux[1] - fancy_flux[1]),color="pink",linestyle="--")
            # F,l_end = polyfit_exps(Ts[21:end],fluxes[21:end],1)
            # plot([Ts[21],Ts[end]],F.+l_end*[Ts[21],Ts[end]].+ (true_flux[1] - fancy_flux[1]),color="pink",linestyle="--")
            #
            # figure("$(flux_label_dict[flux]) feedback",figsize=(20,12))
            # subplot(2,3,subplot_i)
            # plot([2,20],[l_20,l_20],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="--") #averaging_dict[averaging])
            # plot([21,length(Ts)],[l_end,l_end],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="--") #averaging_dict[averaging])

            if average_toggle
              figure("$(area) gregories $(flux)",figsize=(20,12))
              F_20,l_20 = polyfit_exps(fancy_T[1:14],fancy_flux[1:14],1)
              plot([fancy_T[1],fancy_T[14]],F_20 .+l_20*[fancy_T[1],fancy_T[14]].+ (true_flux[1] - fancy_flux[1]),color=color_line,linestyle="-",zorder=4)
              F_end,l_end = polyfit_exps(fancy_T[15:end],fancy_flux[15:end],1)
              if subplot_i == 6
                plot([fancy_T[15],fancy_T[end]],F_end .+l_end*[fancy_T[15],fancy_T[end]].+ (true_flux[1] - fancy_flux[1]),color=color_line,linestyle="-",label=method_dict[method],zorder=4)
              else
                plot([fancy_T[15],fancy_T[end]],F_end .+l_end*[fancy_T[15],fancy_T[end]].+ (true_flux[1] - fancy_flux[1]),color=color_line,linestyle="-",zorder=4)
              end



              figure("$(area) feedback $(flux)",figsize=(20,12))
              subplot(2,3,subplot_i)
              plot([2,20],[l_20,l_20],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="-") #averaging_dict[averaging])
              if subplot_i == 6
                plot([21,length(Ts)],[l_end,l_end],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle=linestyle_dict[flux],label=method_dict[method]) #averaging_dict[averaging])
              else
                plot([21,length(Ts)],[l_end,l_end],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle=linestyle_dict[flux]) #averaging_dict[averaging])
              end
            end
            if method == :true
              if average_toggle
                figure("$(area) residual $(flux)",figsize=(20,12))
                subplot(2,3,subplot_i)

                interp_linear_1 = LinearInterpolation([fancy_T[1],fancy_T[14]],F_20.+l_20*[fancy_T[1],fancy_T[14]],extrapolation_bc = Line())
                interp_linear_2 = LinearInterpolation([fancy_flux[15],Ts[end]],F_end .+l_end*[fancy_flux[15],Ts[end]],extrapolation_bc = Line())
                title("piecewise")
                plot([-3,43],[0,0],"k-")
                plot(vcat(fancy_flux[1:14] .- interp_linear_1.(fancy_T[1:14]),fancy_flux[15:end] .- interp_linear_2.(fancy_T[15:end])),"ko")
                ylim(-.7,.7)
              else
                figure("$(area) residual $(flux)",figsize=(20,12))
                subplot(2,3,subplot_i)
                title("loess")
                plot([-3,43],[0,0],"k-")
                plot(fancy_flux .- Nhat_interp.(get_T_i.(fancy_T)),"ko")
                ylim(-.7,.7)
              end
            end
          end
        end

        # figure("$(area) abrupt4x feedback time series $(flux)",figsize=(20,12))
        # plot([3,150],[0,0],"k--")
        # # markersize = flux in ["N","LW_cloud","SW_cloud"] ? 6 : 4
        #
        # subplot(2,3,subplot_i)
        #
        #
        #
        # plot(1:length(true_flux),(true_flux .- true_flux[1])./(Ts .- Ts[1]),"k-")
        #
        # for method in [:global,:local,:full]
          # est_flux = to_array(to_dset(to_estimated_fn(model_name,"abrupt4x",(method != :full ? Symbol("local_$(method)") : method),"12x24_equaldist")),"$(flux)_$(method != :global ? "monthly" : "annual")","$(area == :special ? "special_" : "")global")
        #   plot(1:length(est_flux),((est_flux .- est_flux[1])./(Ts .- Ts[1])),color=color_dict[method],linestyle="-")
        # end

        if flux == "N"
          figure("$(area) feedback $(flux)",figsize=(20,12))
        else
          figure("$(area) feedback others",figsize=(20,12))
        end
        title(model_name)
        # ylim(-4,4)
        # xlim(2,150)
        if subplot_i in [2,3,5,6]
          yticks([])
        end
        if subplot_i < 4
          xticks([])
        end
        if subplot_i > 3
          xlabel("time (yr)")
        end
        if subplot_i in [1,4]
          ylabel("feedback (\$Wm^{-2}K^{-1}\$)")
        end
        figure("$(area) gregories $(flux)")
        title(model_name)
        # ylim(-4,4)
        # xlim(2,150)
        if subplot_i in [2,3,5,6]
          yticks([])
        end
        if subplot_i < 4
          xticks([])
        end
        if subplot_i > 3
          xlabel("\$\\bar{T}'\$ (K)")
        end
        if subplot_i in [1,4]
          ylabel("\$\\bar{N}\$ (\$Wm^{-2}\$)")
        end
        figure("$(area) residual $(flux)")
      end
    end

    for flux in ["N","SW_cloud"]
      for fig in ["$(area) feedback $(flux)","$(area) gregories $(flux)","$(area) residual $(flux)"]
        figure(fig)
        # suptitle(fig)
        legend(loc="upper right",fontsize=20,bbox_to_anchor=(1.55,1.4),frameon=false,numpoints=1)
        subplots_adjust(left=0.06,hspace=0.15,right=0.85,wspace=0.1)
        savefig("../figs/$(replace(fig," "=>"_"))_cfmip.eps")
      end
    end
  end
end

# fractional feedback time series
if false
  using GeoStats
  using LocallyWeightedRegression
  using Interpolations
  var_names = ["N"] #,"LW_clear","SW_clear","LW_cloud","SW_cloud"]
  for area in [:all] #,:special]
    for (subplot_i,model_name) in enumerate(["CCSM3","CESM104","GISSE2R","HadCM3L","IPSLCM5A","MPIESM12"]) # model_names_abrupt4x
      println(model_name)
      dset = to_dset(to_anomalized_fn(model_name,"abrupt4x","12x24_equaldist"))
      monthly_Ts = to_array(dset,"T_surf","$(area == :special ? "special_" : "")global")
      Ts = to_gregory_dict(dset,"N")["T_surf"]
      fa = fancy_average_maker(1000,50)
      fancy_T = fa(Ts)[2:end] # fancy_average(Ts)
      extra = 0.1
      min_T = min(fancy_T...) - extra
      max_T = max(fancy_T...) + extra
      points = 400
      ns = 30
      varss = 20
      T_grid = range(min_T, stop=max_T, length=points)

      ts = 1:length(fancy_T)

      get_T(t) = fancy_T[t]
      get_T_i(T) = 1+((T-min_T)/(max_T-min_T))*(points-1)

      for flux in var_names
        true_flux = 0
        true_feedbacks = 0
        for (method_i,method) in enumerate([:true,:full,:global,:local]) # ,:full_tropics
          figure("$(area) feedback $(flux)",figsize=(20,12))
          subplot(2,3,subplot_i)

          if method == :true
            monthly_flux = to_array(to_dset(to_anomalized_fn(model_name,"abrupt4x","12x24_equaldist")),flux,"$(area == :special ? "special_" : "")global")
            fluxes = to_annual_average(monthly_flux)
            plot([20,20],[-3.0,3.0],"k--",zorder=0.5)
          else
            if method != :full_tropics
              efn = to_estimated_fn(model_name,"abrupt4x",(method != :full ? Symbol("local_$(method)") : method),"12x24_equaldist")
            else
              efn = to_estimated_fn(model_name,"abrupt4x",:full,"12x24_equaldist",tropics=true)
            end
            # fluxes = to_array(to_dset(efn),"$(flux)_$(method != :global ? "monthly" : "annual")","$(area == :special ? "special_" : "")global")
            fluxes = to_array(to_dset(efn),"$(flux)_seasonal","$(area == :special ? "special_" : "")global")
          end
          println(method)

          fancy_flux = fa(fluxes)[2:end]
          if method == :true
            true_flux = fancy_flux
          end


          color_line = method == :full_tropics ? color_dict[:full] : color_dict[method]
          average_toggle = false
          if !average_toggle
            # smooth 'it
            solution = solve(EstimationProblem(PointSetData(Dict(:y => fancy_flux), transpose(fancy_T)), RegularGrid((min_T,), (max_T,), dims=(points,)), :y), LocalWeightRegress(:y => (neighbors=ns,variogram=ExponentialVariogram(range=varss/10),))) #()))
            Nhat, Nvar = solution[:y]

            # get feedback
            Nhat_interp = interpolate(Nhat,BSpline(Linear()))
            # interp_linear = LinearInterpolation((x[2:end].+x[1:end-1])./2, (vs[2:end].-vs[1:end-1])./(x[2:end].-x[1:end-1]), extrapolation_bc = Line())
            feedbacks = [(Nhat_interp(get_T_i(get_T(t)+.01))-Nhat_interp(get_T_i(get_T(t))))/.01 for t in ts]

            if method == :true
              true_feedbacks = feedbacks
            else
              if subplot_i == 6
                plot(fa(1:length(Ts))[2:end],abs.(feedbacks.-true_feedbacks)./abs.(true_feedbacks),color=color_line,linestyle=(method == :full_tropics ? "--" : "-"),linewidth=(method in [:true,:full] ? 2 : 1),zorder=(2-.0001*method_i)) #averaging_dict[averaging])
                if subplot_i == 6
                  plot([1000,1000],[1000,1001],linestyle="-",color=color_dict[method],marker="o",label=method_dict[method],zorder=(2-.0001*method_i),markersize=6)
                end
              else
                plot(fa(1:length(Ts))[2:end],abs.(feedbacks.-true_feedbacks)./abs.(true_feedbacks),color=color_line,linewidth=(method in [:true,:full] ? 2 : 1),linestyle=(method == :full_tropics ? "--" : "-"),zorder=(2-.0001*method_i)) #averaging_dict[averaging])
              end
              # scale_1 = exp(mean(log.([2,20])))
              # scale_2 = exp(mean(log.([21,1000])))
              # plot(scale_1,feedback_scale(model_name,method,flux,1),color=color_line,linestyle="none",marker="o",zorder=(2-.0001*method_i),markersize=6)
              # plot(scale_2,feedback_scale(model_name,method,flux,2),color=color_line,linestyle="none",marker="o",zorder=(2-.0001*method_i),markersize=6)
            end
          end
          ylim(-3,3)
          xscale("log")
          xlim(2,1000)

          if method != :full_tropics

            figure("$(area) gregories $(flux)",figsize=(20,12))
            subplot(2,3,subplot_i)
            plot(Ts[1:1000],fluxes[1:1000] .+ (true_flux[1] - fancy_flux[1]),linestyle="none",color="gray",marker="o",markersize=1,zorder=0.5)
            method_markersize = Dict(
              # :true => 8,
              :true => 6,
              :global => 4,
              :local => 4,
              :full => 6
            )
            plot(fancy_T,fancy_flux .+ (true_flux[1] - fancy_flux[1]),linestyle="none",color=color_dict[method],marker="o",markersize=method_markersize[method],zorder=(2-.0001*method_i))
            if !average_toggle
              plot(T_grid,Nhat.+ (true_flux[1] - fancy_flux[1]),color=color_line,linestyle="-",zorder=(2-.0001*method_i))
            end
            if subplot_i == 6
              plot([1000,1000],[1000,1001],linestyle="-",color=color_dict[method],marker="o",label=method_dict[method],zorder=(2-.0001*method_i))
            end
            ylim(-1,10)
            xlim(0,8)
            plot([0,8],[0,0],"k-",zorder=0.5)
            # ylim(-1,6)
            # xlim(2,8)
            # plot([2,8],[0,0],"k-",zorder=0.5)
            # F,l_20 = polyfit_exps(Ts[2:20],fluxes[2:20],1)
            # plot([Ts[2],Ts[20]],F.+l_20*[Ts[2],Ts[20]].+ (true_flux[1] - fancy_flux[1]),color="pink",linestyle="--")
            # F,l_end = polyfit_exps(Ts[21:end],fluxes[21:end],1)
            # plot([Ts[21],Ts[end]],F.+l_end*[Ts[21],Ts[end]].+ (true_flux[1] - fancy_flux[1]),color="pink",linestyle="--")
            #
            # figure("$(flux_label_dict[flux]) feedback",figsize=(20,12))
            # subplot(2,3,subplot_i)
            # plot([2,20],[l_20,l_20],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="--") #averaging_dict[averaging])
            # plot([21,length(Ts)],[l_end,l_end],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="--") #averaging_dict[averaging])

            if average_toggle
              figure("$(area) gregories $(flux)",figsize=(20,12))
              F_20,l_20 = polyfit_exps(fancy_T[1:14],fancy_flux[1:14],1)
              plot([fancy_T[1],fancy_T[14]],F_20 .+l_20*[fancy_T[1],fancy_T[14]].+ (true_flux[1] - fancy_flux[1]),color=color_line,linestyle="-",zorder=4)
              F_end,l_end = polyfit_exps(fancy_T[15:end],fancy_flux[15:end],1)
              if subplot_i == 6
                plot([fancy_T[15],fancy_T[end]],F_end .+l_end*[fancy_T[15],fancy_T[end]].+ (true_flux[1] - fancy_flux[1]),color=color_line,linestyle="-",label=method_dict[method],zorder=4)
              else
                plot([fancy_T[15],fancy_T[end]],F_end .+l_end*[fancy_T[15],fancy_T[end]].+ (true_flux[1] - fancy_flux[1]),color=color_line,linestyle="-",zorder=4)
              end

              figure("$(area) feedback $(flux)",figsize=(20,12))
              subplot(2,3,subplot_i)
              plot([2,20],[l_20,l_20],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="-") #averaging_dict[averaging])
              if subplot_i == 6
                plot([21,length(Ts)],[l_end,l_end],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="-",label=method_dict[method]) #averaging_dict[averaging])
              else
                plot([21,length(Ts)],[l_end,l_end],color=color_dict[method],linewidth=(method in [:true,:full] ? 2 : 1),linestyle="-") #averaging_dict[averaging])
              end
            end
            if method == :true
              if average_toggle
                figure("$(area) residual $(flux)",figsize=(20,12))
                subplot(2,3,subplot_i)

                interp_linear_1 = LinearInterpolation([fancy_T[1],fancy_T[14]],F_20.+l_20*[fancy_T[1],fancy_T[14]],extrapolation_bc = Line())
                interp_linear_2 = LinearInterpolation([fancy_flux[15],Ts[end]],F_end .+l_end*[fancy_flux[15],Ts[end]],extrapolation_bc = Line())
                title("piecewise")
                plot([-3,43],[0,0],"k-")
                plot(vcat(fancy_flux[1:14] .- interp_linear_1.(fancy_T[1:14]),fancy_flux[15:end] .- interp_linear_2.(fancy_T[15:end])),"ko")
                ylim(-.7,.7)
              else
                figure("$(area) residual $(flux)",figsize=(20,12))
                subplot(2,3,subplot_i)
                title("loess")
                plot([-3,43],[0,0],"k-")
                plot(fancy_flux .- Nhat_interp.(get_T_i.(fancy_T)),"ko")
                ylim(-.7,.7)
              end
            end
          end
        end

        # figure("$(area) abrupt4x feedback time series $(flux)",figsize=(20,12))
        # plot([3,150],[0,0],"k--")
        # # markersize = flux in ["N","LW_cloud","SW_cloud"] ? 6 : 4
        #
        # subplot(2,3,subplot_i)
        #
        #
        #
        # plot(1:length(true_flux),(true_flux .- true_flux[1])./(Ts .- Ts[1]),"k-")
        #
        # for method in [:global,:local,:full]
          # est_flux = to_array(to_dset(to_estimated_fn(model_name,"abrupt4x",(method != :full ? Symbol("local_$(method)") : method),"12x24_equaldist")),"$(flux)_$(method != :global ? "monthly" : "annual")","$(area == :special ? "special_" : "")global")
        #   plot(1:length(est_flux),((est_flux .- est_flux[1])./(Ts .- Ts[1])),color=color_dict[method],linestyle="-")
        # end

        figure("$(area) feedback $(flux)")
        title(model_name)
        # ylim(-4,4)
        # xlim(2,150)
        if subplot_i in [2,3,5,6]
          yticks([])
        end
        if subplot_i < 4
          xticks([])
        end
        if subplot_i > 3
          xlabel("time (yr)")
        end
        if subplot_i in [1,4]
          ylabel("feedback (\$Wm^{-2}K^{-1}\$)")
        end
        figure("$(area) gregories $(flux)")
        title(model_name)
        # ylim(-4,4)
        # xlim(2,150)
        if subplot_i in [2,3,5,6]
          yticks([])
        end
        if subplot_i < 4
          xticks([])
        end
        if subplot_i > 3
          xlabel("\$\\bar{T}'\$ (K)")
        end
        if subplot_i in [1,4]
          ylabel("\$\\bar{N}\$ (\$Wm^{-2}\$)")
        end
        figure("$(area) residual $(flux)")
      end
    end

    for flux in var_names
      for fig in ["$(area) feedback $(flux)","$(area) gregories $(flux)","$(area) residual $(flux)"]
        figure(fig)
        # suptitle(fig)
        legend(loc="upper right",fontsize=20,bbox_to_anchor=(1.55,1.4),frameon=false,numpoints=1)
        subplots_adjust(left=0.06,hspace=0.15,right=0.85,wspace=0.1)
        savefig("../figs/$(replace(fig," "=>"_")).eps")
      end
    end
  end
end

# seasonal zonal plot
if false
  rc("font", size=12)
  full_dset =  to_dset("/project2/abbot/jsbj/spatial/netcdfs/12x24_equaldist/feedbacks/mmm_control_global_local_multiple_regridded.nc")
  full_dset_local_nonlocal = to_dset("/project2/abbot/jsbj/spatial/netcdfs/12x24_equaldist/feedbacks/mmm_control_full_local_nonlocal.nc")
  plot(full_dset.lat.values,mean((full_dset.N_seasonal_02_feedback - full_dset.N_seasonal_04_feedback).values,dims=2),"k-",label="MR")
  plot(full_dset.lat.values,mean((full_dset_local_nonlocal.N_seasonal_02_feedback_local - full_dset_local_nonlocal.N_seasonal_04_feedback_local).values,dims=2),"r-",label="MR (local)")
  plot(full_dset.lat.values,mean((full_dset_local_nonlocal.N_seasonal_02_feedback_nonlocal - full_dset_local_nonlocal.N_seasonal_04_feedback_nonlocal).values,dims=2),"b-",label="MR (nonlocal)")
  title("JJA - DJF")
  legend(loc="upper right")
  xlabel("latitude")
  ylabel("feedback (\$Wm^{-2}K^{-1}\$)")
  subplots_adjust(left=0.17)
  savefig("../figs/zonal_mean_seasons.eps")
end

# seasonal zonal plot
if false
    plot_lat(values,format) = plot(sin.(true_dset.lat.values*π/180),mean(values,dims=2),format,label="true",linewidth=2)
  contribution_dset =  to_dset("/project2/abbot/jsbj/spatial/netcdfs/12x24_equaldist/contributions/mmm_abrupt4x.nc")
  plot([-1,1],[0,0],"k-")
  plot(contribution_dset.lat.values,mean(contribution_dset.N_seasonal.values,dims=2),"k-",label="MR")
  # plot(full_dset.lat.values,mean((full_dset_local_nonlocal.N_seasonal_02_feedback_local - full_dset_local_nonlocal.N_seasonal_04_feedback_local).values,dims=2),"r-",label="MR (local)")
  # plot(full_dset.lat.values,mean((full_dset_local_nonlocal.N_seasonal_02_feedback_nonlocal - full_dset_local_nonlocal.N_seasonal_04_feedback_nonlocal).values,dims=2),"b-",label="MR (nonlocal)")
  title("contributions")
  legend(loc="lower right")
  xlabel("latitude")
  ylabel("feedback (\$Wm^{-2}K^{-1}\$)")
  subplots_adjust(left=0.1)
  savefig("../figs/zonal_mean_seasons.eps")
end

if false
  true_dset = to_dset("/project2/abbot/jsbj/spatial/netcdfs/12x24_equaldist/finite_difference/mmm_abrupt4x.nc")
  MR_dset = to_dset("/project2/abbot/jsbj/spatial/netcdfs/12x24_equaldist/finite_difference/mmm_abrupt4x_full.nc")
  # MR_tropics_dset = to_dset("/project2/abbot/jsbj/spatial/netcdfs/12x24_equaldist/finite_difference/mmm_abrupt4x_full_tropics.nc")
  true_net(i) = sum([@eval true_dset.$(Symbol(flux,"_",i)).values for flux in [:LW_clear,:SW_clear,:LW_cloud,:SW_cloud]])
  MR_net(i,extra) = sum([@eval MR_dset.$(Symbol(flux,extra*"_",i)).values for flux in [:LW_clear,:SW_clear,:LW_cloud,:SW_cloud]])


  plot_lat(values,format) = plot(sin.(true_dset.lat.values*π/180),mean(values,dims=2),format,label="true",linewidth=2)
  # MR_tropics_net(i,extra) = sum([@eval MR_tropics_dset.$(Symbol(flux,extra*"_",i)).values for flux in [:LW_clear,:SW_clear,:LW_cloud,:SW_cloud]])

  plot_lat(true_net(2)./true_dset.T_surf_2_global.values-true_net(1)./true_dset.T_surf_1_global.values,"k-")
  plot_lat(MR_net(2,"_seasonal")./true_dset.T_surf_2_global.values-net(1,"_seasonal")./true_dset.T_surf_1_global.values,"g-")


  plot([-1,1],[0,0],"k-")
  plot_lat(MR_dset.LW_clear_seasonal_2.values./true_dset.T_surf_2_global.values-MR_dset.LW_clear_seasonal_1.values./true_dset.T_surf_1_global.values,"g--")
  plot_lat(MR_dset.LW_cloud_seasonal_2.values./true_dset.T_surf_2_global.values-MR_dset.LW_cloud_seasonal_1.values./true_dset.T_surf_1_global.values,"g:")
  plot_lat(MR_dset.SW_clear_seasonal_2.values./true_dset.T_surf_2_global.values-MR_dset.SW_clear_seasonal_1.values./true_dset.T_surf_1_global.values,"g-.")
  plot_lat(MR_dset.SW_cloud_seasonal_2.values./true_dset.T_surf_2_global.values-MR_dset.SW_cloud_seasonal_1.values./true_dset.T_surf_1_global.values,"g-")


  plot_lat(true_dset.LW_clear_2.values./true_dset.T_surf_2_global.values-true_dset.LW_clear_1.values./true_dset.T_surf_1_global.values,"k--")
  plot_lat(true_dset.LW_cloud_2.values./true_dset.T_surf_2_global.values-true_dset.LW_cloud_1.values./true_dset.T_surf_1_global.values,"k:")
  plot_lat(true_dset.SW_clear_2.values./true_dset.T_surf_2_global.values-true_dset.SW_clear_1.values./true_dset.T_surf_1_global.values,"k-.")
  plot_lat(true_dset.SW_cloud_2.values./true_dset.T_surf_2_global.values-true_dset.SW_cloud_1.values./true_dset.T_surf_1_global.values,"k-")
  ylim(-5,5)
  xlim(-1,1)
  # plot(sin.(true_dset.lat.values*π/180),mean(MR_tropics_net(2,"_seasonal")./true_dset.T_surf_2_global.values-net(1,"_seasonal")./true_dset.T_surf_1_global.values,dims=2),"g--",label="MR (tropics)")
  lats = vcat(-90,(true_dset.lat.values[2:end]+true_dset.lat.values[1:(end-1)])/2,90)
  xticks(sin.(lats*π/180),lats)
  # plot(full_dset.lat.values,mean((full_dset.N_seasonal_02_feedback - full_dset.N_seasonal_04_feedback).values,dims=2),"k-",label="MR")
  # plot(full_dset.lat.values,mean((full_dset_local_nonlocal.N_seasonal_02_feedback_local - full_dset_local_nonlocal.N_seasonal_04_feedback_local).values,dims=2),"r-",label="MR (local)")
  # plot(full_dset.lat.values,mean((full_dset_local_nonlocal.N_seasonal_02_feedback_nonlocal - full_dset_local_nonlocal.N_seasonal_04_feedback_nonlocal).values,dims=2),"b-",label="MR (nonlocal)")
  # title("JJA - DJF")
  # legend(loc="lower right")
  # xlabel("latitude")
  # ylabel("feedback (\$Wm^{-2}K^{-1}\$)")
  # subplots_adjust(left=0.1)
  # savefig("../figs/zonal_mean_seasons.eps")
end

# tables
if false
  function feedback_error(method,flux,scale_i,averaging="seasonal",area=:all)
    sqrt(mean([(feedback_scale(model_name,method,flux,scale_i,averaging,area) - feedback_scale(model_name,:true,flux,scale_i,averaging,area)).^2 for model_name in all_model_names]))
  end

  averaged_error_line(time_i,grid) = join([join([rpad(round(to_score_2_per_model("mmm","abrupt4x",flux,averaging,grid,method,false,time_i),digits=2),4,"0") for (method,averaging) in [(:local_global,"month_agnostic"),(:local_local,"annual"),(:full,"monthly")]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"
  averaged_error_line(grid) = join([join([rpad(round(to_score_2_per_model_change("mmm","abrupt4x",flux,averaging,grid,method,false),digits=2),4,"0") for (method,averaging) in [(:local_global,"month_agnostic"),(:local_local,"annual"),(:full,"monthly")]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"
  averaged_special_error_line(time_i,grid) = join([join([rpad(round(to_special_score_2_per_model("mmm","abrupt4x",flux,averaging,grid,method,false,time_i),digits=2),4,"0") for (method,averaging) in [(:local_global,"month_agnostic"),(:local_local,"annual"),(:full,"monthly")]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"
  averaged_special_error_line(grid) = join([join([rpad(round(to_special_score_2_per_model_change("mmm","abrupt4x",flux,averaging,grid,method,false),digits=2),4,"0") for (method,averaging) in [(:local_global,"month_agnostic"),(:local_local,"annual"),(:full,"monthly")]]," & ") for flux in ["N","LW_clear", "SW_clear", "LW_cloud", "SW_cloud"]]," & ") * " \\\\"

  scales = ["early","late","change"]

  # for area in [:special] #[:all] #,:special]
  #   println(area)
  #   for scale_i in 1:3
  #     println(scales[scale_i]," & ",join([join([round(feedback_error(method,flux,scale_i,:seasonal,area),digits=2) for method in [:full,:global,:local]]," & ") for flux in ["N","LW_clear","LW_cloud","SW_clear","SW_cloud"]]," & "))
  #     # for averaging in [:annual,:seasonal,:monthly,:month_agnostic]
  #     #   println(averaging,": ",join([round(feedback_error(method,flux,scale_i,averaging,area),digits=2) for method in [:global,:local,:full]],", "))
  #     # end
  #   end
  # end

  # # which averaging?
  # println(" & ",join(["\\multicolumn{3}{c|}{\\emph{$averaging}}" for averaging in ["seasonal","annual","all months","monthly clim."]]," & ")," \\\\")
  # println(" & ",join(["MR & global & local" for averaging in ["seasonal","annual","all months","monthly clim."]]," & ")," \\\\")
  # println("score 1")
  # for scale_i in 1:3
  #   println(scales[scale_i]," & ",join([join([round(feedback_error(method,"N",scale_i,averaging),digits=2) for method in [:full,:global,:local]]," & ") for averaging in [:seasonal,:annual,:monthly,:month_agnostic]]," & "))
  # end
  # println("score 2")
  # for scale_i in 1:2
  #   println(scales[scale_i]," & ",join([join([round(to_score_2_per_model("mmm","abrupt4x","N",averaging,"12x24_equaldist",method,false,scale_i),digits=2) for method in [:full,:local_global,:local_local]]," & ") for averaging in [:seasonal,:annual,:monthly,:month_agnostic]]," & "))
  # end
  # println(scales[3]," & ",join([join([round(to_score_2_per_model_change("mmm","abrupt4x","N",averaging,"12x24_equaldist",method,false),digits=2) for method in [:full,:local_global,:local_local]]," & ") for averaging in [:seasonal,:annual,:monthly,:month_agnostic]]," & "))

  # println("all")
  # for scale_i in 1:2
  #   println(scales[scale_i]," & ",join([join([round(to_score_2_per_model("mmm","abrupt4x",flux,:seasonal,"12x24_equaldist",method,false,scale_i),digits=2) for method in [:full,:local_global,:local_local]]," & ") for flux in ["N","LW_clear","LW_cloud","SW_clear","SW_cloud"]]," & "))
  #   # for averaging in [:annual,:seasonal,:monthly,:month_agnostic]
  #   #   println(averaging,": ",join([round(feedback_error(method,flux,scale_i,averaging,area),digits=2) for method in [:global,:local,:full]],", "))
  #   # end
  # end
  # println(scales[3]," & ",join([join([round(to_score_2_per_model_change("mmm","abrupt4x",flux,:seasonal,"12x24_equaldist",method,false),digits=2) for method in [:full,:local_global,:local_local]]," & ") for flux in ["N","LW_clear","LW_cloud","SW_clear","SW_cloud"]]," & "))
  #
  # println("special")
  # for scale_i in 1:2
  #   println(scales[scale_i]," & ",join([join([round(to_special_score_2_per_model("mmm","abrupt4x",flux,:seasonal,"12x24_equaldist",method,false,scale_i),digits=2) for method in [:full,:local_global,:local_local]]," & ") for flux in ["N","LW_clear","LW_cloud","SW_clear","SW_cloud"]]," & "))
  #   # for averaging in [:annual,:seasonal,:monthly,:month_agnostic]
  #   #   println(averaging,": ",join([round(feedback_error(method,flux,scale_i,averaging,area),digits=2) for method in [:global,:local,:full]],", "))
  #   # end
  # end
  # println(scales[3]," & ",join([join([round(to_special_score_2_per_model_change("mmm","abrupt4x",flux,:seasonal,"12x24_equaldist",method,false),digits=2) for method in [:full,:local_global,:local_local]]," & ") for flux in ["N","LW_clear","LW_cloud","SW_clear","SW_cloud"]]," & "))
end
