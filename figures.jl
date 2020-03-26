ion()
rc("font", weight="light")
rc("font", size=14)
rc("mathtext", fontset="stixsans")


# @pyimport cartopy.crs as ccrs
# @pyimport cartopy.io.shapereader as shpreader

subplot_dims(n) = (Int(ceil((sqrt(1+4n)-1)/2)), Int(ceil((sqrt(1+4n)-1)/2))+1)
color_dict = Dict(
  "true" => "black",
  # "global_local" => "blue",
  # "global_global" => "green",
  # "global" => "blue",
  # "full" => "orange",
  "global annual" => "blue",
  "global monthly" => "green",
  "full monthly" => "orange",
  "local annual" => "red"
)

lighter_color_dict = Dict(
  "true" => "#888888",
  "global annual" => "#8888FF",
  "global monthly" => "#88D388",
  "full monthly" => "#FFD388",
  "local annual" => "#FF8888"
)

function time_series_plots(data::OrderedDict,anomaly,fn="")
  var_names = anomaly ? ["T","LW","LWCS","LWCRF","SW","SWCS","SWCRF","N"] : ["tas","rlut","rsut","rlutcs","rsutcs"]
  for var_name in var_names
    figure(var_name,figsize=(16,9))
    ndata = length(data)
    rows,cols = subplot_dims(ndata)
  
    run_i = 0
    for (run_name,run_data) in data
      run_i += 1
      subplot(rows,cols,run_i)
      title(run_name)
      for (est_name,est_data) in run_data
        color = color_dict[est_name]
        plot(1:length(est_data[var_name]),est_data[var_name],"k-")
      end
    
      if run_i % cols == 1
        ylabel("$(var_name) ($(var_name in ["T","tas"] ? "K" : "W/m2"))")
      end
      if ndata - run_i < cols
        xlabel("time (years)")
      end
    end
  
    if !isempty(fn)
      savefig("figs/time_series/" * fn * "_" * var_name * ".eps")
    end
  end
end

# gregory_data = OrderedDict{String,Dict}()
# for i in 1:6
#   gregory_data["run$(i)"] = Dict()
#   for est in ["true","full"]
#     gregory_data["run$(i)"][est] = Dict(
#       "temp" => rand(10),
#       "flux" => rand(10)
#     )
#   end
# end
function gregory_plots(data::OrderedDict,fn="")
  println("making Gregory plots of ",fn)
  figure(fn,figsize=(16,7))
  
  ax = nothing
  
  cols = 4
  rows = 3
  for (run_types,model_dict) in data
    run_i = 0
    for (model_name,run_dict) in model_dict
      if fn[1:9] == "long_runs"
        if run_types == "abrupts"
          subplot_i = 1 + (div(run_i,(cols-2)) * cols) + run_i%(cols-2)
        elseif run_types == "ramps"
          subplot_i = 3 + (div(run_i,(cols-3)) * cols) + run_i%(cols-3)
        elseif run_types == "rcps"
          subplot_i = 4
        end
      else
        subplot_i = 1 + run_i
      end
      run_i += 1
      ax = subplot(rows,cols,subplot_i)
      for (run_type,greg_dict) in run_dict
        if run_type == "control"
          plot([-2,15],[0,0],"k-")
          title(model_name)
        end
        
        if fn[end-7:end] == "all_ests"
          for (est_name,est_data) in greg_dict
            color = color_dict[est_name]
            # plot(est_data["T_surf"],est_data["flux"],markersize=1,marker="o",color=color,linestyle="none")
            plot(est_data["T_surf"],est_data["flux"],markersize=1,marker="o",markeredgecolor="none",color="0.5",linestyle="none")
            if run_type == "control"
              flux_0, flux_slope = [ones(est_data["T_surf"]) est_data["T_surf"]]\est_data["flux"]
              # println(run_name," ",flux_slope)
              exts = extrema(est_data["T_surf"])
              plot(exts,flux_0 + collect(exts)*flux_slope,color=color_dict["global"],linewidth=2)
            else
              for time_scale in time_scales_dict[run_type]
                if to_full_time_scale_boolean(est_data["T_surf"],time_scale)
                  Ts = to_time_scale(est_data["T_surf"],time_scale)
                  flux_0, flux_slope = [ones(length(Ts)) Ts]\to_time_scale(est_data["flux"],time_scale)
                  # println(run_name," ",flux_slope)
                  exts = extrema(Ts)
                  plot(exts,flux_0 + collect(exts)*flux_slope,color=color,linewidth=2)
                end
              end
            end
          end
        else
          est_name = "true"
          est_data = greg_dict[est_name]

          color = color_dict[est_name]
          # plot(est_data["T_surf"],est_data["flux"],markersize=1,marker="o",color=color,linestyle="none")
          plot(est_data["T_surf"],est_data["flux"],markersize=1,marker="o",markeredgecolor="none",color="0.5",linestyle="none")
          if run_type == "control"
            flux_0, flux_slope = [ones(est_data["T_surf"]) est_data["T_surf"]]\est_data["flux"]
            # println(run_name," ",flux_slope)
            exts = extrema(est_data["T_surf"])
            plot(exts,flux_0 + collect(exts)*flux_slope,color=color_dict["global"],linewidth=2)
          else
            for time_scale in time_scales_dict[run_type]
              if to_full_time_scale_boolean(est_data["T_surf"],time_scale)
                Ts = to_time_scale(est_data["T_surf"],time_scale)
                flux_0, flux_slope = [ones(length(Ts)) Ts]\to_time_scale(est_data["flux"],time_scale)
                # println(run_name," ",flux_slope)
                exts = extrema(Ts)
                plot(exts,flux_0 + collect(exts)*flux_slope,color=color,linewidth=2)
              end
            end
          end
        end
      end
      if run_types == "abrupts"
        # xlim(-1,8)
        xlim(1,8)
        ylim(-2,7.5)
      elseif run_types == "ramps"
        xlim(-1,5.5)
        ylim(-2,4.8611)
      elseif run_types == "rcps"
        xlim(-1,13)
        ylim(-2,12.77775)
      end
      if fn[end-7:end] == "all_ests"
        old_ylims = ylim()
        ylim(old_ylims[1]+(-(old_ylims...)),old_ylims[2])
      end
    end
  end
  for row_i in 1:rows
    subplot(rows,cols,(row_i-1)*cols+1)
    ylabel("\$R_{net}\$ (W/m2)")
  end
  
  # for (run_type,run_dict) in data_abrupt4x
  #   run_i = 0
  #   for (run_name,run_data) in run_dict
  #     subplot_i = 1 + (div(run_i,(cols-1)) * cols) + run_i%(cols-1)
  #     run_i += 1
  #     subplot(rows,cols,subplot_i)
  #     if run_type == "control"
  #       plot([-100,100],[0,0],"k-")
  #     end
  #     # plot([0,0],"k-")
  #     title(run_name)
  #     for (est_name,est_data) in run_data
  #       color = color_dict[est_name]
  #       # plot(est_data["T_surf"],est_data["flux"],markersize=1,marker="o",color=color,linestyle="none")
  #       plot(est_data["T_surf"],est_data["flux"],markersize=1,marker="o",markeredgecolor="none",color="0.5",linestyle="none")
  #       if run_type == "control"
  #         flux_0, flux_slope = [ones(est_data["T_surf"]) est_data["T_surf"]]\est_data["flux"]
  #         # println(run_name," ",flux_slope)
  #         exts = extrema(est_data["T_surf"])
  #         plot(exts,flux_0 + collect(exts)*flux_slope,"r-",linewidth=2)
  #       else
  #         for time_scale in time_scales_dict["abrupt4x"]
  #           if to_full_time_scale_boolean(est_data["T_surf"],time_scale)
  #             Ts = to_time_scale(est_data["T_surf"],time_scale)
  #             flux_0, flux_slope = [ones(length(Ts)) Ts]\to_time_scale(est_data["flux"],time_scale)
  #             # println(run_name," ",flux_slope)
  #             exts = extrema(Ts)
  #             plot(exts,flux_0 + collect(exts)*flux_slope,"k-",linewidth=2)
  #           end
  #         end
  #       end
  #     end
  #
  #     if subplot_i % cols == 1
  #       ylabel("TOA flux (W/m2)")
  #     end
  #     if subplot_i > cols*(rows-1)
  #       xlabel("surface temp (K)")
  #     end
  #     xlim(-1,8)
  #     ylim(-2,7.5)
  #     # if split(fn," ")[2] == "control"
  #     #   xlim(-0.75,0.75)
  #     #   ylim(-1.5,1.5)
  #     # elseif split(fn," ")[2] == "abrupt4x"
  #     #   xlim(0,10)
  #     #   ylim(-10,10)
  #     # else
  #     #   if run_name == "control"
  #     #     xlim(-0.75,0.75)
  #     #     ylim(-1.5,1.5)
  #     #   else
  #     #     quadruplings = log(4,Int(float(match(r"\d+",run_name).match)))
  #     #     xlim(0,10quadruplings)
  #     #     ylim(-10quadruplings,10quadruplings)
  #     #   end
  #     # end
  #   end
  # end
  #
  # for (run_type,run_dict) in data_1pct2x
  #   run_i = 0
  #   for (run_name,run_data) in run_dict
  #     run_i += 1
  #     subplot_i = run_i*cols
  #     ax = subplot(rows,cols,subplot_i)
  #     if run_type == "control"
  #       plot([-100,100],[0,0],"k-")
  #     end
  #     # plot([0,0],"k-")
  #     title(run_name)
  #     for (est_name,est_data) in run_data
  #       color = color_dict[est_name]
  #       # plot(est_data["T_surf"],est_data["flux"],markersize=1,marker="o",color=color,linestyle="none")
  #       plot(est_data["T_surf"],est_data["flux"],markersize=1,marker="o",markeredgecolor="none",color="0.5",linestyle="none")
  #       if run_type == "control"
  #         flux_0, flux_slope = [ones(est_data["T_surf"]) est_data["T_surf"]]\est_data["flux"]
  #         # println(run_name," ",flux_slope)
  #         exts = extrema(est_data["T_surf"])
  #         plot(exts,flux_0 + collect(exts)*flux_slope,"r-",linewidth=2)
  #       else
  #         for time_scale in time_scales_dict["1pct2x"]
  #           if to_full_time_scale_boolean(est_data["T_surf"],time_scale)
  #             Ts = to_time_scale(est_data["T_surf"],time_scale)
  #             flux_0, flux_slope = [ones(length(Ts)) Ts]\to_time_scale(est_data["flux"],time_scale)
  #             # println(run_name," ",flux_slope)
  #             exts = extrema(Ts)
  #             plot(exts,flux_0 + collect(exts)*flux_slope,"k-",linewidth=2)
  #           end
  #         end
  #       end
  #     end
  #
  #
  #     # if run_i % cols == 1
  #     #   ylabel("TOA flux (W/m2)")
  #     # end
  #     if subplot_i > cols*(rows-1)
  #       xlabel("surface temp (K)")
  #     end
  #     xlim(-1,5.5)
  #     ylim(-2,4.8611)
  #     # if split(fn," ")[2] == "control"
  #     #   xlim(-0.75,0.75)
  #     #   ylim(-1.5,1.5)
  #     # elseif split(fn," ")[2] == "abrupt4x"
  #     #   xlim(0,10)
  #     #   ylim(-10,10)
  #     # else
  #     #   if run_name == "control"
  #     #     xlim(-0.75,0.75)
  #     #     ylim(-1.5,1.5)
  #     #   else
  #     #     quadruplings = log(4,Int(float(match(r"\d+",run_name).match)))
  #     #     xlim(0,10quadruplings)
  #     #     ylim(-10quadruplings,10quadruplings)
  #     #   end
  #     # end
  #   end
  # end
  #
  # ax[:add_line](matplotlib[:lines][:Line2D](-2.5*[1,1],[-42,14],color="black",linewidth=2,clip_on=false))
  # ax[:add_line](matplotlib[:lines][:Line2D](-19.4*[1,1],[-42,14],color="black",linewidth=2,clip_on=false))
  subplots_adjust(hspace=.3,left=0.05,right=0.97)
  # text(-14.7,24,"abrupt4x",fontsize=22)
  # text(1.3,24,"1pct2x",fontsize=22)
  # # ax[:add_line](matplotlib[:lines][:Line2D]([0,10],[.5,.5],clip_on=false))
  # if !isempty(fn)
    # savefig("figs/gregory/" * join(split(fn," "),"_") * ".eps")
    println("saving figs/gregory/$(fn).eps")
    savefig("figs/gregory/$(fn).eps")
  # end
  #
  ax
end

function gregory_plots_abrupts(data::OrderedDict,fn="")
  println("making Gregory plots of ",fn)
  # figure(fn,figsize=(16,9))
  figure(fn,figsize=(12,6))
  
  ax = nothing
  
  markersize = 2
  cols = 3
  rows = 2
  
  for forcing in [false] #[false,true]
    for (run_types,model_dict) in data
      run_i = 0
      for (model_name,run_dict) in model_dict
        subplot_i = 1 + run_i
        run_i += 1
        ax = subplot(rows,cols,subplot_i)
        for (run_type,greg_dict) in run_dict
        
          if fn[end-7:end] == "all_ests"
            if forcing
              F_true = 0
              F_other = 0
            end
            for (est_name,est_data) in greg_dict
              if (run_type == "abrupt4x") && (est_name == "true")
                plot([-2,15],[0,0],"k-")
                title(model_name)
              end
              
              println(est_name)
              color = color_dict[est_name]
              lighter_color = lighter_color_dict[est_name]

              if run_type == "control"
                plot(est_data["T_surf"],est_data["flux"],markersize=markersize,marker="o",markeredgecolor="none",linestyle="none",color=lighter_color_dict["global annual"])
                flux_0, flux_slope = [ones(est_data["T_surf"]) est_data["T_surf"]]\est_data["flux"]
                exts = extrema(est_data["T_surf"])
                plot(exts,flux_0 + collect(exts)*flux_slope,color=color_dict["global annual"],linewidth=2)
              else
                if forcing
                  for time_scale in time_scales_dict[run_type]
                    if to_full_time_scale_boolean(est_data["T_surf"],time_scale)
                      Ts = to_time_scale(est_data["T_surf"],time_scale)
                      if time_scale == time_scales_dict[run_type][1]
                        flux_0, flux_slope = [ones(length(Ts)) Ts]\to_time_scale(est_data["flux"],time_scale)
                        if est_name == "true"
                          F_true = flux_0 + flux_slope * est_data["T_surf"][6]
                        end
                        F_other = flux_0 + flux_slope * est_data["T_surf"][6]
                      end
                      flux_0, flux_slope = [ones(length(Ts)) Ts]\to_time_scale(est_data["flux"],time_scale)
                      # println(run_name," ",flux_slope)
                      exts = extrema(Ts)
                      plot(exts,flux_0 + collect(exts)*flux_slope+F_true-F_other,color=color,linewidth=2,zorder=((est_name == "true") ? 10 : 1)+10)
                    end
                  end
                  plot(est_data["T_surf"],est_data["flux"]+F_true-F_other,markersize=markersize,marker="o",markeredgecolor="none",linestyle="none",color=lighter_color,zorder=((est_name == "true") ? 10 : 1))
                else
                  plot(est_data["T_surf"],est_data["flux"],markersize=markersize,marker="o",markeredgecolor="none",linestyle="none",color=lighter_color,zorder=((est_name == "true") ? 10 : 1))
                  first = true
                  for time_scale in time_scales_dict[run_type]
                    if to_full_time_scale_boolean(est_data["T_surf"],time_scale)
                      Ts = to_time_scale(est_data["T_surf"],time_scale)
                      flux_0, flux_slope = [ones(length(Ts)) Ts]\to_time_scale(est_data["flux"],time_scale)
                      # println(run_name," ",flux_slope)
                      exts = extrema(Ts)
                      plot(exts,flux_0 + collect(exts)*flux_slope,color=color,linewidth=2,zorder=((est_name == "true") ? 10 : 1)+10,label=first ? est_name : nothing)
                      first = false
                    end
                  end
                end
              end
            end
          # else
          #   est_name = "true"
          #   est_data = greg_dict[est_name]
          #
          #   color = color_dict[est_name]
          #   # plot(est_data["T_surf"],est_data["flux"],markersize=1,marker="o",color=color,linestyle="none")
          #   plot(est_data["T_surf"],est_data["flux"],markersize=2,marker="o",markeredgecolor="none",color="0.5",linestyle="none")
          #   if run_type == "control"
          #     flux_0, flux_slope = [ones(est_data["T_surf"]) est_data["T_surf"]]\est_data["flux"]
          #     # println(run_name," ",flux_slope)
          #     exts = extrema(est_data["T_surf"])
          #     plot(exts,flux_0 + collect(exts)*flux_slope,color="green",linewidth=2,linestyle=":")
          #     monthly_slope = mean([to_feedback_dict(to_feedback_fn(model_name,"control",:global_global))[flux_name]["monthly"][i]["feedback"] for i in 1:12])[1]
          #     plot(exts,flux_0 + collect(exts)*monthly_slope,color=color_dict["global"],linewidth=2,linestyle="-")
          #   else
          #     for time_scale in time_scales_dict[run_type]
          #       if to_full_time_scale_boolean(est_data["T_surf"],time_scale)
          #         Ts = to_time_scale(est_data["T_surf"],time_scale)
          #         flux_0, flux_slope = [ones(length(Ts)) Ts]\to_time_scale(est_data["flux"],time_scale)
          #         # println(run_name," ",flux_slope)
          #         exts = extrema(Ts)
          #         plot(exts,flux_0 + collect(exts)*flux_slope,color=color,linewidth=2)
          #       end
          #     end
          #   end
          end
        end

        if fn[1:2] == "ox"
          xlim(0,8)
          ylim(-8,8)
        else
          xlim(-1,8)
          ylim(-10,16)
        end
        # ylim(-1,16)
        # ylim(-2,10)
        # if fn[end-7:end] == "all_ests"
        #   old_ylims = ylim()
        #   ylim(old_ylims[1]+(-(old_ylims...)),old_ylims[2])
        # end
      end
    end


    if fn[1:2] != "ox"
      legend(loc="center right", bbox_to_anchor=(1.55,1.15),numpoints=1,frameon=false,fontsize=12)

      for row_i in 1:rows
        subplot(rows,cols,(row_i-1)*cols+1)
        ylabel("\$R_{4x}\$ (\$W/m^2\$)",labelpad=-10)
      end
  
      for col_i in 1:cols
        subplot(rows,cols,cols*(rows-1)+col_i)
        xlabel("\$\\Delta T_{4x}\$ (\$K\$)")
      end
    else
      for row_i in 1:rows
        subplot(rows,cols,(row_i-1)*cols+1)
        ylabel("\$R'\$ (\$W/m^2\$)",labelpad=-5)
      end
  
      for col_i in 1:cols
        subplot(rows,cols,cols*(rows-1)+col_i)
        xlabel("\$T'\$ (\$K\$)")
      end
    end

    # subplots_adjust(hspace=.3,left=0.04,right=0.87,top=0.95,bottom=0.06)
    subplots_adjust(left=0.08,bottom=0.09,right=0.98,top=0.95,hspace=0.4,wspace=0.25)
    
    println("saving figs/gregory/$(fn).eps")
    savefig("figs/gregory/$(fn).eps")
    ax
  end
end

function gregory_plots_abrupts_w_forcing(data::OrderedDict,fn="")
  println("making Gregory plots of ",fn)
  figure(fn,figsize=(16,9))
  
  ax = nothing
  
  markersize = 2
  
  cols = 3 #4
  rows = 2
  for (run_types,model_dict) in data
    run_i = 0
    for (model_name,run_dict) in model_dict
      subplot_i = 1 + run_i
      run_i += 1
      ax = subplot(rows,cols,subplot_i)
      for (run_type,greg_dict) in run_dict
        if run_type == "control"
          plot([-2,15],[0,0],"k-")
          title(model_name)
        end
        
        if fn[end-7:end] == "all_ests"
          # F_true = 0
          # F_other = 0
          for (est_name,est_data) in greg_dict
            println(est_name)
            color = color_dict[est_name]
            lighter_color = lighter_color_dict[est_name]
            # plot(est_data["T_surf"],est_data["flux"],markersize=1,marker="o",color=color,linestyle="none")

            if run_type == "control"
              plot(est_data["T_surf"],est_data["flux"],markersize=markersize,marker="o",markeredgecolor="none",linestyle="none",color=lighter_color_dict["global annual"])
              flux_0, flux_slope = [ones(est_data["T_surf"]) est_data["T_surf"]]\est_data["flux"]
              # println(run_name," ",flux_slope)
              exts = extrema(est_data["T_surf"])
              plot(exts,flux_0 + collect(exts)*flux_slope,color=color_dict["global annual"],linewidth=2)
              # plot(exts,flux_0 + collect(exts)*flux_slope,color="green",linewidth=2,linestyle=":")
              # monthly_slope = mean([to_feedback_dict(to_feedback_fn(model_name,"control",:global_global))[flux_name]["monthly"][i]["feedback"] for i in 1:12])[1]
              # plot(exts,flux_0 + collect(exts)*monthly_slope,color=color_dict["control"],linewidth=2,linestyle="-")
            else
              for time_scale in time_scales_dict[run_type]
                if to_full_time_scale_boolean(est_data["T_surf"],time_scale)
                  Ts = to_time_scale(est_data["T_surf"],time_scale)
                  if time_scale == time_scales_dict[run_type][1]
                    flux_0, flux_slope = [ones(length(Ts)) Ts]\to_time_scale(est_data["flux"],time_scale)
                    if est_name == "true"
                      F_true = flux_0 + flux_slope * est_data["T_surf"][6]
                    end
                    F_other = flux_0 + flux_slope * est_data["T_surf"][6]
                  end
                  flux_0, flux_slope = [ones(length(Ts)) Ts]\to_time_scale(est_data["flux"],time_scale)
                  # println(run_name," ",flux_slope)
                  exts = extrema(Ts)
                  plot(exts,flux_0 + collect(exts)*flux_slope+F_true-F_other,color=color,linewidth=2,zorder=((est_name == "true") ? 10 : 1)+10)
                end
              end
              plot(est_data["T_surf"],est_data["flux"]+F_true-F_other,markersize=markersize,marker="o",markeredgecolor="none",linestyle="none",color=lighter_color,zorder=((est_name == "true") ? 10 : 1))
              # plot(est_data["T_surf"],est_data["flux"],markersize=markersize,marker="o",markeredgecolor="none",linestyle="none",color=lighter_color,zorder=((est_name == "true") ? 10 : 1))

              # for time_scale in time_scales_dict[run_type]
              #   if to_full_time_scale_boolean(est_data["T_surf"],time_scale)
              #     Ts = to_time_scale(est_data["T_surf"],time_scale)
              #     flux_0, flux_slope = [ones(length(Ts)) Ts]\to_time_scale(est_data["flux"],time_scale)
              #     # println(run_name," ",flux_slope)
              #     exts = extrema(Ts)
              #     plot(exts,flux_0 + collect(exts)*flux_slope,color=color,linewidth=2,zorder=((est_name == "true") ? 10 : 1)+10)
              #   end
              # end
              
              # flux_0, flux_slope = [ones(length(est_data["T_surf"])) est_data["T_surf"]]\est_data["flux"]
              # # println(run_name," ",flux_slope)
              # exts = extrema(est_data["T_surf"])
              # plot(exts,flux_0 + collect(exts)*flux_slope,color=color,linewidth=2,zorder=((est_name == "true") ? 10 : 1)+10)

            end
          end
        # else
        #   est_name = "true"
        #   est_data = greg_dict[est_name]
        #
        #   color = color_dict[est_name]
        #   # plot(est_data["T_surf"],est_data["flux"],markersize=1,marker="o",color=color,linestyle="none")
        #   plot(est_data["T_surf"],est_data["flux"],markersize=2,marker="o",markeredgecolor="none",color="0.5",linestyle="none")
        #   if run_type == "control"
        #     flux_0, flux_slope = [ones(est_data["T_surf"]) est_data["T_surf"]]\est_data["flux"]
        #     # println(run_name," ",flux_slope)
        #     exts = extrema(est_data["T_surf"])
        #     plot(exts,flux_0 + collect(exts)*flux_slope,color="green",linewidth=2,linestyle=":")
        #     monthly_slope = mean([to_feedback_dict(to_feedback_fn(model_name,"control",:global_global))[flux_name]["monthly"][i]["feedback"] for i in 1:12])[1]
        #     plot(exts,flux_0 + collect(exts)*monthly_slope,color=color_dict["global"],linewidth=2,linestyle="-")
        #   else
        #     for time_scale in time_scales_dict[run_type]
        #       if to_full_time_scale_boolean(est_data["T_surf"],time_scale)
        #         Ts = to_time_scale(est_data["T_surf"],time_scale)
        #         flux_0, flux_slope = [ones(length(Ts)) Ts]\to_time_scale(est_data["flux"],time_scale)
        #         # println(run_name," ",flux_slope)
        #         exts = extrema(Ts)
        #         plot(exts,flux_0 + collect(exts)*flux_slope,color=color,linewidth=2)
        #       end
        #     end
        #   end
        end
      end
      xlim(-1,8)
      ylim(-10,10)
      # ylim(-1,16)
      # ylim(-2,10)
      # if fn[end-7:end] == "all_ests"
      #   old_ylims = ylim()
      #   ylim(old_ylims[1]+(-(old_ylims...)),old_ylims[2])
      # end
    end
  end
  for row_i in 1:rows
    subplot(rows,cols,(row_i-1)*cols+1)
    ylabel("\$R'\$ (\$W/m^2\$)")
  end
  
  for col_i in 1:cols
    subplot(rows,cols,cols*(rows-1)+col_i)
    xlabel("\$T'\$ (K)")    
  end

  subplots_adjust(hspace=.3,left=0.05,right=0.97)
  println("saving figs/gregory/$(fn)_w_forcing.eps")
  savefig("figs/gregory/$(fn)_w_forcing.eps")
  ax
end


function normal_score_plots(data)
  plot(a,quantile(Normal(),ecdf(a)(a)),"bo")
end

function residual_plots(data,plot_type,fn="")
  figure(figsize=(16,9))
  ndata = length(data)
  rows,cols = subplot_dims(ndata)
  
  run_i = 0
  for (run_name,run_data) in data
    run_i += 1
    subplot(rows,cols,run_i)
    # plot([-100,100],[0,0],"k-")
    # plot([0,0],"k-")
    title(run_name)
    if plot_type == :normal_score
      plot(run_data["residuals"],quantile(Normal(),ecdf(run_data["residuals"])(run_data["residuals"])),markersize=1,"ko")
    elseif plot_type == :against_T_surfs
      plot(run_data["T_surfs"],run_data["residuals"],markersize=1,"ko")      
    end
  end
  
  if !isempty(fn)
    savefig("figs/residual_plots/$(plot_type)/" * join(split(fn," "),"_") * ".eps")
  end
end

feedback_limits = Dict(
  "LW_clear" => (-3,-1),
  "SW_clear" => (-.5,1.5),
  "N_clear" => (-2.5,-.5),
  "LW_cloud" => (-.5,1.5),
  "SW_cloud" => (-1,1),
  "N_cloud" => (-1,1),
  "LW" => (-2.5,-.5),
  "SW" => (-.5,1.5),
  "N" => (-2.25,-.25)
)

# var_markers = Dict(
# "LW_clear" => "^",
# "SW_clear" => "<",
# "LW_cloud" => ">",
# "SW_cloud" => "v",
# "N" => "o"
# )
#
#
# function scatter_plot(data::OrderedDict,full_type=:global_local,fn="")
#   long_run = fn == "long_run"
#
#   # for avg_type in ["monthly","annual"]
#   avg_type = "monthly"
#
#   figure(figsize=(8,8))
#
#   score_dict(1,avg_type,"12x24_equaldist",fn,full_type)
#   tmp_score_dict = to_score_dict(to_score_fn(1,avg_type,"12x24_equaldist",fn,full_type))
#
#   counter = 150
#   for var_name in ["LW_clear","SW_clear","LW_cloud","SW_cloud"] #var_names
#     # subplot(counter+=1)
#     for est_type in ["global","local","full"]
#       # println(collect(keys(data[all_model_names[1]]["true"][var_name])))
#       # est_part = collect(est_dict[var_name])[long_run ? end : 2][2]
#
#       plot(
#         [collect(data[model_name]["true"][var_name])[end][2]["annual"]["feedback"][1] for model_name in model_names_abrupt4x],
#         [collect(data[model_name][est_type][var_name])[end][2][avg_type]["feedback"][1] for model_name in model_names_abrupt4x],
#         color=color_dict[est_type],
#         linestyle="none",marker=var_markers[var_name],
#         label = est_type != "true" ? est_type * " $(round(tmp_score_dict[var_name][est_type],2))" : nothing
#       )
#     end
#   end
#
#
#   # for i in 1:5
#   #   subplot(1,5,i)
#     # xlim(.5,length(data)+.5)
#     # title(replace(var_names[i],"_"," "))
#     # plot([-2,3],[-2,3],"k--")
#     # xlim(-2,3)
#     # ylim(-2,3)
#     xlim_save = xlim()
#     ylim_save = ylim()
#     plot(xlim_save,xlim_save,"k--")
#     xlim(xlim_save)
#     ylim(ylim_save)
#
#
#     # if i == 1
#     #   # ylabel("feedback (W/m2/K)")
#     # else
#     #   # yticks([],[])
#     # end
#     # xticks(1:length(data),collect(keys(data)),rotation=90)
#     # legend(numpoints=1,loc=(var_names[i] == "LW_clear" ? "upper right" : "lower right"))
#     # legend(numpoints=1,loc=(var_names[i] == "SW_clear" ? "lower left" : "lower right"))
#     # ylim(feedback_limits[vars[i]])
#   # end
#
#   subplots_adjust(left=.05,top=.93,right=0.99,bottom=0.13,wspace=0.12)
#
#   if !isempty(fn)
#     savefig("figs/scatter_plots/feedbacks/" * join(split(fn," "),"_") * "_" * avg_type * "_" * string(full_type) * ".eps")
#   end
# end

function scatter_plot(data::OrderedDict,full_type=:global_local,fn="")
  long_run = fn == "long_run"
  
  # for avg_type in ["monthly","annual"]
  avg_type = "monthly"
  
  figure(figsize=(19,5.5))
  
  score_dict(1,avg_type,"12x24_equaldist",fn,full_type)
  tmp_score_dict = to_score_dict(to_score_fn(1,avg_type,"12x24_equaldist",fn,full_type))

  counter = 150
  for var_name in var_names
    subplot(counter+=1)
    for est_type in ["global","local","full"]
      # println(collect(keys(data[all_model_names[1]]["true"][var_name])))
      # est_part = collect(est_dict[var_name])[long_run ? end : 2][2]
      
      plot(
        [collect(data[model_name]["true"][var_name])[end][2]["annual"]["feedback"][1] for model_name in model_names_abrupt4x],
        [collect(data[model_name][est_type][var_name])[end][2][avg_type]["feedback"][1] for model_name in model_names_abrupt4x],
        color=color_dict[est_type],
        linestyle="none",marker="o",
        label = est_type != "true" ? est_type * " $(round(tmp_score_dict[var_name][est_type],2))" : nothing
      )
    end
  end  
  

  for i in 1:5
    subplot(1,5,i)
    # xlim(.5,length(data)+.5)
    title(replace(var_names[i],"_"," "))
    # plot([-2,3],[-2,3],"k--")
    # xlim(-2,3)
    # ylim(-2,3)
    xlim_save = xlim()
    ylim_save = ylim()
    plot(xlim_save,xlim_save,"k--")
    xlim(xlim_save)
    ylim(ylim_save)


    # if i == 1
    #   # ylabel("feedback (W/m2/K)")
    # else
    #   # yticks([],[])
    # end
    # xticks(1:length(data),collect(keys(data)),rotation=90)
    # legend(numpoints=1,loc=(var_names[i] == "LW_clear" ? "upper right" : "lower right"))
    legend(numpoints=1,loc=(var_names[i] == "SW_clear" ? "lower left" : "lower right"))
    # ylim(feedback_limits[vars[i]])
  end

  subplots_adjust(left=.05,top=.93,right=0.99,bottom=0.13,wspace=0.12)
  
  if !isempty(fn)
    savefig("figs/scatter_plots/feedbacks/" * join(split(fn," "),"_") * "_" * avg_type * "_" * string(full_type) * ".eps")
  end
end

function feedback_dot_plots(data::OrderedDict,full_type=:global_local,fn="")
  if fn[1:2] == "ox"
    ox = true
    fn = fn[4:end]
  else
    ox = false
  end
  long_run = fn == "long_run"
  
  var_name_dict = Dict(
    "N" => "net",
    "SW_clear" => "SW clear",
    "LW_clear" => "LW clear",
    "SW_cloud" => "SW cloud",
    "LW_cloud" => "LW cloud"
  )
  
  label_dict = Dict(
    "true" => "true",
    "full" => "MR",
    "global" => "global",
    "local" => "local"
  )
        
  j = 0
  for time_scale in time_scales_dict["abrupt4x"]
    first_time_scale, last_time_scale = split(time_scale," to ")
    # if last_time_scale == "end"
    #   last_time_scale = "\$n_{years,4x}\$"
    # end
    
    # est_name_with_symbol = Dict(
    #   "true" => "true \$\\lambda_{4x,$(first_time_scale),$(last_time_scale)}\$",
    #   "global" => "global \$\\hat{\\lambda}_{global}\$",
    #   "local" => "local \$\\hat{\\lambda}_{4x,local,$(first_time_scale),$(last_time_scale)}\$",
    #   "full" => "full \$\\hat{\\lambda}_{4x,full,$(first_time_scale),$(last_time_scale)}\$"
    # )
    # est_name_with_symbol = Dict(
    #   "true" => "true \$\\lambda_{4x}\$",
    #   "global_annual" => "global \$\\hat{\\lambda}_{4x,glob-ann}\$",
    #   "global_monthly" => "global \$\\hat{\\lambda}_{4x,glob-mon}\$",
    #   "local_annual" => "local \$\\hat{\\lambda}_{4x,loc-ann}\$",
    #   "full_monthly" => "full \$\\hat{\\lambda}_{4x,full-mon}\$"
    # )
    
    j += 1
    # if j== 1
    #   continue
    # end
    # for avg_type in ["monthly","annual"]
    # avg_type = "monthly"
      if ox
        figure(1,figsize=(19,8.5))
      else
        figure(figsize=(19,8.5))
      end
      model_i = 0

      for i in 1:5
        subplot(1,5,i)
        plot([.5,10.5],[0,0],"k-")
        fill_between([1.5,2.5],[-5,-5],[5,5],color="0.8")
        fill_between([3.5,4.5],[-5,-5],[5,5],color="0.8")
        fill_between([5.5,6.5],[-5,-5],[5,5],color="0.8")
        fill_between([7.5,8.5],[-5,-5],[5,5],color="0.8")
        if fn == "long_run"
          fill_between([9.5,10.5],[-5,-5],[5,5],color="0.8")          
        end
      end


      for (model_name,model_data) in data
        model_i += 1
        last = (model_i == 6)

        counter = 150
        for var_name in var_names
          subplot(counter+=1)
        
          for (est_type,est_dict) in model_data
            function tmp(avg_type)
              # if est_type != "true"
              #   score_dict(1,avg_type,j == 1 ? "abrupt4x_starts" : fn,full_type)
              #   tmp_score_dict = to_score_dict(to_score_fn(1,avg_type,j == 1 ? "abrupt4x_starts" : fn,full_type))
              # end
              
              ratio = .39
                        
              est_part = est_dict[var_name][time_scale]
              color = color_dict[est_type == "true" ? "true" : est_type * " " * avg_type]
              # label = (last ? est_name_with_symbol[est_type * "_" * avg_type] * ((est_type != "true") ? "\n($(round(tmp_score_dict[var_name][est_type],2)) \$W/m^2/K\$)" : "") : nothing)
              # label = ((last && (est_type != "true")) ? est_type * " " * avg_type * "\n($(round(tmp_score_dict[var_name][est_type],2)) \$W/m^2/K\$)" : nothing)

              # if last  && ox
              #   if est_type != "true" && j == 1
              #     label = time_scale #* "\n($(round(tmp_score_dict[var_name][est_type],2)) \$W/m^2/K\$)"
              #     plot(100,100,label=label,linestyle="none")
              #   end
              # end

              if last
                if !ox
                  label = label_dict[est_type] * ((est_type != "true") ? " " : "") # * "\n($(round(tmp_score_dict[var_name][est_type],2)) \$W/m^2/K\$)" : "")
                else
                  if j == 1
                    if est_type == "true"
                      label = "true"
                    else
                      label = "estimate" # time_scale * "\n($(round(tmp_score_dict[var_name][est_type],2)) \$W/m^2/K\$)"
                    end
                  else
                    label = nothing
                  end
                end
              else
                label = nothing
              end
              

              
              plot(model_i+(ox ? ratio*(j-1.5) : 0),est_part[avg_type]["feedback"][1],marker="o",color=color,linestyle="none",markersize=(est_type == "true" ? 9 : 6),label=label)
              if last
                # if (est_type == "true")
                #   plot(100,100,label="  ",linestyle="none")
                # elseif ((est_type == "global") && (avg_type == "monthly"))
                #   for m in 1:1 #2
                #     plot(100,100,label="  ",linestyle="none")
                #   end
                # end
                
                if ox
                  if est_type == "true" && j == 1
                    # plot(100,100,label="  ",linestyle="none")
                  else
                    if est_type != "true" && j == 2
                      label = time_scale * "" #"\n($(round(tmp_score_dict[var_name][est_type],2)) \$W/m^2/K\$)"
                      plot(100,100,label=label,linestyle="none")                      
                    end 
                  end
                end
              end
              plot([model_i+(ox ? ratio*(j-1.5) : 0),model_i+(ox ? ratio*(j-1.5) : 0)],est_part[avg_type]["feedback"][1]+[-1,1]*est_part[(est_type == "true") ? "annual" : avg_type]["feedback_unct"][1],color=color)
            end

      
            avg_type = (est_type in ["true","global","local"]) ? "annual" : "monthly"
            tmp(avg_type)
          end
        end
      end

      for i in 1:5
        subplot(1,5,i)
        xlim(.5,length(data)+.5)
        title(var_name_dict[var_names[i]])
        # ylim(-2,2)
        if ox
          ylim(-2.5,1.5)
        else
          ylim(-3,3)
        end
        # ylim(-4,4)
        if i == 1
          # ylabel("\$abrupt4_x\$ feedback, years $(first_time_scale) to $(last_time_scale) (\$W/m^2/K\$)")
          # ylabel("abrupt4\$_x\$ feedback, years $(time_scale) (\$W/m^2/K\$)")
          # ylabel("feedback, years $(time_scale) (\$W/m^2/K\$)")
          ylabel("feedback (\$Wm^{-2}K^{-1}\$)")
        else
          yticks([],[])
        end
        xticks(1:length(data),collect(keys(data)),rotation=45)
        legend(numpoints=1,loc=(var_names[i] == "LW_clear" ? "upper right" : "lower right"), fontsize=12,handlelength=1,handletextpad=0.5,labelspacing=0.1,ncol=2,columnspacing=0.5) #,title="\$\\bullet\$true")
        # ylim(feedback_limits[vars[i]])
      end

      subplots_adjust(left=.05,top=.9,right=0.99,bottom=0.13,wspace=0.12)
      # suptitle(replace(fn,"_"," ")* " " * string(full_type))
      if !ox
        suptitle("years $(time_scale)",fontsize=18)
      end
      if !isempty(fn)
        println("figs/dot_plots/feedbacks/$(ox ? "ox_" : "")" * join(split(fn," "),"_") * "_" * replace(time_scale," ","_") * ".eps")
        savefig("figs/dot_plots/feedbacks/$(ox ? "ox_" : "")" * join(split(fn," "),"_") * "_" * replace(time_scale," ","_") * ".eps")
      end
  end
end

function feedback_change_dot_plots(data::OrderedDict,full_type=:global_local,fn="")
  long_run = fn == "long_run"
        
  var_name_dict = Dict(
    "N" => "net",
    "SW_clear" => "SW clear",
    "LW_clear" => "LW clear",
    "SW_cloud" => "SW cloud",
    "LW_cloud" => "LW cloud"
  )
        
  figure(figsize=(19,8.5))
  model_i = 0

  for i in 1:5
    subplot(1,5,i)
    plot([.5,10.5],[0,0],"k-")
    fill_between([1.5,2.5],[-5,-5],[5,5],color="0.8")
    fill_between([3.5,4.5],[-5,-5],[5,5],color="0.8")
    fill_between([5.5,6.5],[-5,-5],[5,5],color="0.8")
    fill_between([7.5,8.5],[-5,-5],[5,5],color="0.8")
    if fn == "long_run"
      fill_between([9.5,10.5],[-5,-5],[5,5],color="0.8")          
    end
  end


  for (model_name,model_data) in data
    model_i += 1
    last = (model_i == 6)

    counter = 150
    for var_name in var_names
      subplot(counter+=1)
        
      for (est_type,est_dict) in model_data
        function tmp(avg_type)
          # if est_type != "true"
          #   score_dict(3,avg_type,fn,full_type)
          #   tmp_score_dict = to_score_dict(to_score_fn(3,avg_type,fn,full_type))
          # end
                      
          est_part = collect(est_dict[var_name])[2][2][avg_type]["feedback"][1] - collect(est_dict[var_name])[1][2][avg_type]["feedback"][1]
          color = color_dict[est_type == "true" ? "true" : est_type * " " * avg_type]
          # label = (last ? est_name_with_symbol[est_type * "_" * avg_type] * ((est_type != "true") ? "\n($(round(tmp_score_dict[var_name][est_type],2)) \$W/m^2/K\$)" : "") : nothing)
          # label = ((last && (est_type != "true")) ? est_type * " " * avg_type * "\n($(round(tmp_score_dict[var_name][est_type],2)) \$W/m^2/K\$)" : nothing)
          label = (last ? (est_type == "full" ? "MR" : est_type) : nothing)  # * ((est_type != "true") ? "" : "" * avg_type * "\n($(round(tmp_score_dict[var_name][est_type],2)) \$W/m^2/K\$)" : "") : nothing)
          plot(model_i,est_part,marker="o",color=color,linestyle="none",markersize=(est_type == "true" ? 9 : 6),label=label)
          if last
            if (est_type == "true")
              # plot(100,100,label="  ",linestyle="none")
            elseif ((est_type == "global") && (avg_type == "monthly"))
              for m in 1:1 #2
                plot(100,100,label="  ",linestyle="none")
              end
            end
          end
          # plot([model_i,model_i],est_part[avg_type]["feedback"][1]+[-1,1]*est_part[(est_type == "true") ? "annual" : avg_type]["feedback_unct"][1],color=color)
        end

        avg_type = (est_type in ["true","global","local"]) ? "annual" : "monthly"
        tmp(avg_type)
      end
    end
  end

  for i in 1:5
    subplot(1,5,i)
    xlim(.5,length(data)+.5)
    title(var_name_dict[var_names[i]])
    ylim(-1.,1.)
    if i == 1
      ylabel("change in feedback (\$Wm^{-2}K^{-1}\$)")
    else
      yticks([],[])
    end
    xticks(1:length(data),collect(keys(data)),rotation=45)
    legend(numpoints=1,loc="lower right", fontsize=12,handlelength=1,handletextpad=0.5,labelspacing=0.1,ncol=2,columnspacing=0.5)      # ylim(feedback_limits[vars[i]])
  end

  subplots_adjust(left=.05,top=.9,right=0.99,bottom=0.13,wspace=0.12)
  # suptitle(replace(fn,"_"," ")* " " * string(full_type))
  if !isempty(fn)
    savefig("figs/dot_plots/feedbacks/abrupt4x_change.eps")
  end
end

function error_dot_plots(data::OrderedDict,fn="")
  long_run = fn == "long_run"
  
  var_name_dict = Dict(
    "N" => "net",
    "SW_clear" => "SW clear",
    "LW_clear" => "LW clear",
    "SW_cloud" => "SW cloud",
    "LW_cloud" => "LW cloud"
  )
    
  j = 0
  for time_scale in time_scales_dict["abrupt4x"]
    j += 1

    # for avg_type in ["annual"] #["monthly"] # ,"annual"]
      # avg_type = "monthly"
      figure(figsize=(19,7))
      model_i = 0

      for i in 1:5
        subplot(1,5,i)
        plot([.5,10.5],[0,0],"k-")
        fill_between([1.5,2.5],[-5,-5],[65,65],color="0.8")
        fill_between([3.5,4.5],[-5,-5],[65,65],color="0.8")
        fill_between([5.5,6.5],[-5,-5],[65,65],color="0.8")
        fill_between([7.5,8.5],[-5,-5],[65,65],color="0.8")
        if fn == "long_run"
          fill_between([9.5,10.5],[-5,-5],[65,65],color="0.8")          
        end
      end
        
      for (model_name,model_data) in data
        model_i += 1
        last = (model_i == 6)

        counter = 150
        for var_name in var_names
          subplot(counter+=1)
          for (both_type,est_dict) in model_data
            est_type,avg_type = split(both_type," ")
            score_dict(2,avg_type,j == 1 ? "abrupt4x_starts" : fn)
            tmp_score_dict = to_score_dict(to_score_fn(2,avg_type,j == 1 ? "abrupt4x_starts" : fn))
            
            color = color_dict[est_type == "true" ? "true" : both_type]
            label = ((est_type != "true") && last) ? both_type * "\n($(round(tmp_score_dict[var_name][est_type],2)) \$W/m^2/K\$)" : nothing
            # label = (last ? est_type * ((est_type != "true") ? " " * avg_type * "\n($(round(tmp_score_dict[var_name][est_type],2)) \$W/m^2/K\$)" : "") : nothing)
            plot(model_i,est_dict[var_name][j],marker="o",color=color,linestyle="none",markersize=(est_type == "true" ? 9 : 6),label=label)
              # plot([x,x],dict[(est_type == "true") ? "annual" : avg_type]["feedback"][1]+[-1,1]*dict[(est_type == "true") ? "annual" : avg_type]["feedback_unct"][1],color=color)

            # text(x-.09,est_type == "true" ? feedback_limits[var_name][2] : feedback_limits[var_name][1]-.05,round(dict["R2"][1],1) == 1 ? "1" : string(round(dict["R2"][1],1))[2:end],color=color)
          
          end
        end
      end

      for i in 1:5
        subplot(1,5,i)
        xlim(.5,length(data)+.5)
        title(var_name_dict[var_names[i]])
        ylim(0,16) # 60) # avg_type == "monthly" ? 7 : 9)
        if i == 1
          ylabel("finite difference error (\$W/m^2/K\$)")
          # ylabel("area-weighted root mean square error difference in \$R\$ (\$W/m^2/K\$)")
        else
          yticks([],[])
        end
        xticks(1:length(data),collect(keys(data)),rotation=45)
        legend(numpoints=1,loc="upper right",fontsize=12,handlelength=1,handletextpad=0.5,labelspacing=0.1,ncol=2,columnspacing=0.5)
        # legend(numpoints=1,loc=(var_names[i] == "LW_clear" ? "upper right" : "lower right"), fontsize=12,handlelength=1,handletextpad=0.5,labelspacing=0.1,ncol=2,columnspacing=0.5) #,title="\$\\bullet\$true")
        # ylim(feedback_limits[vars[i]])
      end

      subplots_adjust(left=.05,top=.9,right=0.99,bottom=0.14,wspace=0.12)
      suptitle("years $(time_scale)",fontsize=18)
      if !isempty(fn)
        savefig("figs/dot_plots/errors/" * join(split(fn," "),"_") * "_$(replace(time_scale," ","_")).eps")
      end
    # end
  end
end

# function maps(data::OrderedDict,fn="")
#   figure(figsize=(16,9))
#   ndata = length(data)
#   rows,cols = subplot_dims(ndata)
#
#   run_i = 0
#   for (run_name,run_data) in data
#     run_i += 1
#     subplot(rows,cols,run_i)
#     title(run_name)
#
#     ax = axes(projection=ccrs.LambertCylindrical(central_longitude=180.0))
#     ax[:coastlines](linewidth=0.2)
#     ax[:set_adjustable]("datalim")
#     ax[:set_aspect]("auto")
#     pcolormesh(run_data["lons"],run_data["lats"],run_data["values"],cmap="bwr",transform=ccrs.PlateCarree()) # ,vmin=-0.03,vmax=0.03
#   end
#
#   if !isempty(fn)
#     savefig(fn)
#   end
# end

function asya_gregories(data::OrderedDict,fn="")
  
  if !isempty(fn)
    savefig(fn)
  end
end

function asya_maps(data::OrderedDict,fn="")

  if !isempty(fn)
    savefig(fn)
  end
end

function area_weighted_cdf(model_name,run_type,grid="orig",estimate_type=:none)
  dset = to_dset(to_finite_difference_fn(model_name,run_type,grid,estimate_type))
  orig_dset = to_dset(to_finite_difference_fn(model_name,run_type))
  lats = dset["lat"][:values]
  nlat = length(lats)
  nlon = length(dset["lon"][:values])
  
  areas = fill(1/nlon,(nlat,nlon))
  weights = sin.(vcat([-90],(lats[2:end]+lats[1:end-1])/2,[90])*π/180)
  areas.*=(weights[2:end] - weights[1:end-1])/2
  areas,nlat,nlon = to_vector(areas)
  
  if (grid == "orig") && (estimate_type == :none)
  
    for t in 1:3
      figure(figsize=(16,9))

      suptitle("$(model_name) $(run_type) $(annual_time_scales[t])")
    

      flux = "T_surf_$(t)"
      fluxes,nlat,nlon = to_vector(to_array(orig_dset,flux))
      subplot(3,4,1)
      plot([0,1],[0,0],"k-")
      plot([0,1],[dset[flux*"_global"][:values][1],dset[flux*"_global"][:values][1]],"b--")

      s = sortperm(fluxes)
      xs = zeros(2*(nlat*nlon))
      xs[2:2:end] = areas[s]
      ys = zeros(2*(nlat*nlon))
      ys[1:2:end] = ys[2:2:end] = fluxes[s]
      plot(cumsum(xs),ys,"b-")
      fill_between(cumsum(areas[s]),fluxes[s],zeros(nlat*nlon),color="0.8")        
      xlim(0,1)
    
      for freq_i in 1:3
        freq = freqs[freq_i]
        for sky_i in 1:3
          sky = skies[sky_i]
          flux = freq*sky*"_$(t)"
          fluxes,nlat,nlon = to_vector(to_array(dset,flux))
          subplot(3,4,(4*(sky_i-1))+1+freq_i)

          plot([0,1],[0,0],"k-")
          plot([0,1],[to_array(dset,flux,:global)[1],to_array(dset,flux,:global)[1]],"b--")

          s = sortperm(fluxes)
          xs = zeros(2*(nlat*nlon))
          xs[2:2:end] = areas[s]
          ys = zeros(2*(nlat*nlon))
          ys[1:2:end] = ys[2:2:end] = fluxes[s]
          plot(cumsum(xs),ys,"b-")
          fill_between(cumsum(areas[s]),fluxes[s],zeros(nlat*nlon),color="0.8")        
          xlim(0,1)
        end
      end
    
      savefig("figs/area_weighted_cdfs/$(model_name)_$(run_type)_$(estimate_type == :none ? "true" : estimate_type).eps")
    end
  else
    for avg_type in ["annual","monthly"]
      for t in 1:3
        figure(figsize=(16,9))

        suptitle("$(model_name) $(run_type) $(annual_time_scales[t])")
    

        flux = "T_surf_$(t)"
        fluxes,nlat,nlon = to_vector(to_array(orig_dset,flux))
        subplot(3,4,1)
        plot([0,1],[0,0],"k-")
        plot([0,1],[orig_dset[flux*"_global"][:values][1],orig_dset[flux*"_global"][:values][1]],"b--")

        s = sortperm(fluxes)
        xs = zeros(2*(nlat*nlon))
        xs[2:2:end] = areas[s]
        ys = zeros(2*(nlat*nlon))
        ys[1:2:end] = ys[2:2:end] = fluxes[s]
        plot(cumsum(xs),ys,"b-")
        fill_between(cumsum(areas[s]),fluxes[s],zeros(nlat*nlon),color="0.8")        
        xlim(0,1)
    
        for freq_i in 1:3
          freq = freqs[freq_i]
          for sky_i in 1:3
            sky = skies[sky_i]
            flux = freq*sky*"_$(avg_type)_$(t)"
            # println(dset)
            fluxes,nlat,nlon = to_vector(to_array(dset,flux))
            subplot(3,4,(4*(sky_i-1))+1+freq_i)

            plot([0,1],[0,0],"k-")
            plot([0,1],[to_array(dset,flux,:global)[1],to_array(dset,flux,:global)[1]],"b--")

            s = sortperm(fluxes)
            xs = zeros(2*(nlat*nlon))
            xs[2:2:end] = areas[s]
            ys = zeros(2*(nlat*nlon))
            ys[1:2:end] = ys[2:2:end] = fluxes[s]
            plot(cumsum(xs),ys,"b-")
            fill_between(cumsum(areas[s]),fluxes[s],zeros(nlat*nlon),color="0.8")        
            xlim(0,1)
          end
        end
    
        savefig("figs/area_weighted_cdfs/$(model_name)_$(run_type)_$(estimate_type == :none ? "true" : estimate_type)_$(avg_type).eps")
      end
    end
  end
end
