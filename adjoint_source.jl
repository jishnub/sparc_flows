using PyCall
using FITSIO
using Polynomials
using StatsBase
# using Dierckx
# @pyimport matplotlib.pyplot as plt
using PyPlot

macro forward_dir(x)
  return :(joinpath(directory,"forward_src$(dec(srcno,2))_ls$(dec(jobno,2))",$x))
end
macro adjoint_dir(x)
  return :(joinpath(directory,"adjoint_src$(dec(srcno,2))",$x))
end
macro data_dir(x)
  return :(joinpath(directory,"data",$x))
end

function main()

  unshift!(PyVector(pyimport("sys")["path"]), "")
  @pyimport read_params
  directory = read_params.get_directory()

  procno::AbstractString = get(ENV,"PBS_VNODENUM","0")
  srcno::Int = parse(Int,procno)%8+1
  jobno::Int = div(parse(Int,procno),8)
  linesearch::Bool =jobno>0

  iwls::Bool =false
  prev_iter_exist::Bool = false
  prevtau::Float64 = 1
  iwls_pow::Float64 = 2
  iwls_eps::Float64 = 0.25/60.

  Lx::Float64 = read_params.get_xlength()
  nx::Int = read_params.get_nx()
  outputcad::Float64 = read_params.get_dt()

  x = linspace(0,1,nx+1)[1:end-1]

  dx::Float64 = x[2]-x[1]
  eyekh::Array{Complex{Float64},1}= fftfreq(nx)*2π*1im*nx

  master_pixels::Array{Float64,1} = squeezeto1D(readdlm(joinpath(directory,"master.pixels")))
  x00 = master_pixels[srcno]

  distances::Array{Float64,1} = abs((x-0.5)*Lx - x00)
  signed_distances::Array{Float64,1} = (x-0.5)*Lx - x00

  acc::Array{Float64,3} = fitsread(@forward_dir "vz_cc.fits")

  dat::Array{Float64,3} = fitsread(@data_dir "$(dec(srcno,2)).fits")

  nt = size(dat)[end]


  t = readdlm(@forward_dir "timeline")[:,2]
  t = -t[end:-1:1]/60
  dt_min = t[2]-t[1]
  dt_sec = dt_min*60
  total_time_sec = dt_sec*nt

  vzcc_fft::Array{Complex{Float64},3} = fft(acc)
  data_fft::Array{Complex{Float64},3} = fft(dat)

  spectral_ridges=Dict{Int64,Dict{ASCIIString,Any}}(0=>Dict("name"=>"fmode"))
  for i=1:7
    spectral_ridges[i]=Dict("name"=>"p$(i)mode")
  end

  for modeno in keys(spectral_ridges)

    if !isfile(joinpath(directory,"params.$modeno"))
      continue
    end

    filter_name = "$(spectral_ridges[modeno]["name"])_filter"
    if isfile("$(filter_name).fits")
      spectral_ridges[modeno]["filter"] = fitsread("$(filter_name).fits")
    else
      spectral_ridges[modeno]["filter"] = eval(parse("$filter_name"))(nt,outputcad,nx,Lx)::Array{Float64,3}
      fitswrite("$(filter_name).fits",spectral_ridges[modeno]["filter"])
    end

    if isfile(joinpath(directory,"filter.params.$(spectral_ridges[modeno]["name"])"))
      leftcorner,rightcorner = squeezeto1D(readdlm(
                  joinpath(directory,"filter.params.$(spectral_ridges[modeno]["name"])")))
      filt = freq_filter(leftcorner,rightcorner,nt)
      for i=1:nt
        spectral_ridges[modeno]["filter"][i] =  spectral_ridges[modeno]["filter"][i]*filt[i]
      end
    end
  end

  for (i,vel) in enumerate([0.5,0.75,0.95,1.15,1.2,1.4,1.7,1.9])
    spectral_ridges[i-1]["velocity"] = vel
  end

  adj = zeros(nx,1,nt)
  misfit::Float64 = 0.0
  nmeasurements::Int = 0

  for modeno in keys(spectral_ridges)

    # Compute travel time difference
    if !isfile(joinpath(directory,"params.$modeno"))
      continue
    end

    mindist,maxdist,window = readdlm(joinpath(directory,"params.$modeno"))

    halftime = round(Int,window/(2.*dt_min))

    vzcc_filtered::Array{Float64,3} = real(bfft(vzcc_fft .* spectral_ridges[modeno]["filter"]))
    data_filtered::Array{Float64,3} = real(bfft(data_fft .* spectral_ridges[modeno]["filter"]))
    ccdot::Array{Float64,3} = fft_time_deriv(vzcc_filtered,total_time_sec)

    if !linesearch
      windowfile = open(joinpath(directory,"forward_src$(dec(srcno,2))_ls00","windows.$modeno"),"w")
    else
      windowfile = open(joinpath(directory,"forward_src$(dec(srcno,2))_ls00","windows.$(modeno)"),"r")
    end

    ttdiff_file = open((@forward_dir "ttdiff.$modeno"),"w")


    xindex_list = Int[]

    tau_list = Float64[]
    # tau_orig_list = Float64[]
    tau_err_low_list = Float64[]
    tau_err_high_list = Float64[]

    left_time_cutoff = Int[]
    right_time_cutoff = Int[]

    for (xindex,dist_from_src) in enumerate(distances)

      if !(mindist<dist_from_src<maxdist) continue
      end

      if !linesearch
        if spectral_ridges[modeno]["name"]=="highpmode"
          timest_index = round(Int,floor((-0.000417147774671*dist_from_src^2+
          0.313350998096*dist_from_src+17.9609631186)/dt_min))
          timefin_index = round(Int,floor((-0.00028361249034*dist_from_src^2+
          0.29270337114*dist_from_src+36.4398734578)/dt_min))
        else
          timest_index = round(Int,floor(dist_from_src/(spectral_ridges[modeno]["velocity"]*dt_min)))
          timefin_index = timest_index + 40
        end

        # vzcc_in_range = slice(vzcc_filtered,xindex,1,timest_index:timefin_index)
        # loc_orig = maxloc(abs(real(vzcc_in_range)))+timest_index-1
        # lef_orig = loc_orig - halftime
        # rig_orig = loc_orig + halftime

      end

      tau_different_window = Float64[]
      for loc=(timest_index+halftime):(timefin_index-halftime)

        lef = loc - halftime
        rig = loc + halftime

        # tau = compute_tt_quadratic_fit(
        # slice(vzcc_filtered,xindex,1,:),slice(data_filtered,xindex,1,:),lef,rig,dt_min)
        tau = compute_tt_gizonbirch(slice(vzcc_filtered,xindex,1,:),
        slice(data_filtered,xindex,1,:),lef,rig,dt_min)

        push!(tau_different_window,tau)
      end

      # tau_orig = compute_tt_gizonbirch(slice(vzcc_filtered,xindex,1,:),
      # slice(data_filtered,xindex,1,:),lef_orig,rig_orig,dt_min)

      tau_lower_quantile = quantile(tau_different_window,0.16)
      tau_upper_quantile = quantile(tau_different_window,0.84)
      tau_median = median(tau_different_window)
      tau_err_low = tau_median-tau_lower_quantile
      tau_err_high = tau_upper_quantile - tau_median

      push!(xindex_list,xindex)
      push!(tau_list,tau_median)
      # push!(tau_orig_list,tau_orig)
      push!(tau_err_low_list,tau_err_low)
      push!(tau_err_high_list,tau_err_high)
      push!(left_time_cutoff,timest_index)
      push!(right_time_cutoff,timefin_index)

      ttstr = @sprintf("%3d %14.8f %14.8f %14.8f %4d %4d\n",xindex,tau_median*60,
      tau_err_low*60,tau_err_high*60,timest_index,timefin_index)
      write(ttdiff_file,ttstr)

    end # loop over distances

    errorbar(signed_distances[xindex_list],tau_list*60,
    yerr=(tau_err_low_list*60,tau_err_high_list*60),color="red",marker=".",ls="None",label="median")

    # plot(signed_distances[xindex_list],tau_orig_list*60,color="b",marker=".",ls="None",label="original")

    plot_filename = @forward_dir "traveltime_$(spectral_ridges[modeno]["name"]).png"
    xlabel("Distance from source (Mm)")
    ylabel(L"$\Delta\tau$ (sec)")
    savefig(plot_filename)
    # legend(loc="best")
    clf()
    # plt.xlabel("Distance from source (Mm)")
    # plt.ylabel(L"$\Delta\tau$ (sec)")
    # plt.savefig(plot_filename)
    # plt.clf()

    close(ttdiff_file)
    close(windowfile)

    nmeasurements+=length(tau_list)

    tau_err_avg = (tau_err_low_list+tau_err_high_list)/2

    misfit += sum((tau_list./tau_err_avg).^2)/length(tau_list)


    # Adjoint source computation
    for (i,xindex) in enumerate(xindex_list)
      window_bandpass = zeros(nt)
      time_range = left_time_cutoff[i]:right_time_cutoff[i]
      window_bandpass[time_range] = 1.0
      filtex = (slice(ccdot,xindex,1,:).*window_bandpass)[end:-1:1]*nx*nt

      con = -(tau_list[i]/tau_err_avg[i]^2)/sum(slice(ccdot,xindex,1,time_range).^2*dt_min)
      adj[xindex,1,:] += reshape(filtex*con,1,1,nt)

    end
  end # loop over spectral ridges

  adj*=nt / 60

  if !linesearch
    fitswrite((@adjoint_dir "source.fits"),adj)

    f=open((@adjoint_dir "information"),"w")
    write(f,"$nt\n")
    write(f,"$(t[1]*60)\n")
    write(f,"$(t[nt]*60)\n")
    write(f,"2D")
    close(f)
  end

  f=open(joinpath(directory,"kernel","misfit_$(dec(srcno,2))_$(dec(jobno,2))"),"w")
  println(@sprintf "%25.16E" misfit)
  write(f,@sprintf "%12d%12d%26.16E" srcno nmeasurements misfit)
  close(f)
end

function fitswrite(filename,array)
  f=FITS(filename,"w")
  write(f,array)
  close(f)
  println("Written file $(filename)")
end

function fitsread(filename)
  f=FITS(filename,"r")
  arr=read(f[1])
  close(f)
  println("Reading $(filename)")
  return arr
end

function squeezeto1D(arr)
  singleton = find(x->(x==1),[size(arr)...])
  return squeeze(arr,tuple(singleton...))
end

function simps(f,dt::Float64)

  int_f = 0.0
  for k=1:length(f)

    if (k == 1) || (k == length(f))
        int_f += f[k]
        continue
    end

    if (mod(k,2) == 0)
        int_f += 4*f[k]
    else
        int_f += 2*f[k]
    end
  end

  int_f*=dt/3
  return int_f
end

function compute_tt_gizonbirch(u0,u,t_low_index,t_high_index,dt_sec)

  nt = size(u)[1]
  window = zeros(nt)
  window[t_low_index:t_high_index] = 1

  total_time_sec = dt_sec*nt
  u0dot = fft_time_deriv(u0,total_time_sec)

  return simps(window.*u0dot.*(u0-u),dt_sec)/simps(window.*u0dot.^2,dt_sec)
end

function compute_tt_quadratic_fit(u0, u,t_low_index,t_high_index,dt_min)
  u0 = u0[t_low_index:t_high_index];
  u  =  u[t_low_index:t_high_index];
  nt = t_high_index-t_low_index+1
  cc = crosscov(u0,u,collect(-(nt-1):(nt-1)),demean=false);

  t = float(collect(-nt+1:nt-1))*dt_min

  cc_max_index::Int = indmax(cc);

  t_wavepacket = slice(t,cc_max_index-1:cc_max_index+1);
  cc_wavepacket = slice(cc,cc_max_index-1:cc_max_index+1);

  p = polyfit(t_wavepacket,cc_wavepacket,2);
  return -p[1]/(2*p[2])
end


function fft_time_deriv(arr::Union{Array{Float64,1},
                        SubArray{Float64,1,Array{Float64,3},Tuple{Int64,Int64,Colon},3}},
                        total_time_sec::Float64)
  N = length(arr)
  darr = zeros(Complex{Float64},N)
  arr_fft = fft(arr)
  twopij = 2π*1im/total_time_sec
  for k=0:div(N,2)-1
    arr_fft[k+1]*=twopij*k
  end

  arr_fft[div(N,2)+1]=0

  for k=div(N,2)+1:N-1
    arr_fft[k+1]*=twopij*(k-N)
  end

  return real(ifft(arr_fft))
end


function fft_time_deriv(arr::Array{Float64,3},total_time_sec::Float64)

  two_pi_j_by_T = 2π*1im/total_time_sec
  N=size(arr)[3]
  freq_deriv = zeros(Complex{Float64},(1,1,div(N,2)+1))

  for k=0:div(N,2)-1
    freq_deriv[1,1,k+1]=k
  end

  freq_deriv[1,1,div(N,2)+1]=0

  return real(irfft(rfft(arr,3).*freq_deriv*two_pi_j_by_T,N,3))
end

macro modefilter()

  return quote
    dx = Lx/nx
    k = abs(fftfreq(nx,dx)) * 2π
    f = abs(fftfreq(nt,dt))*1e3

    disp_rel_0=sum([Polylow[i]*k.^(i-1) for i=1:length(Polylow)])
    disp_rel_1=sum([Poly[i]*k.^(i-1) for i=1:length(Poly)])

    filter_array = zeros(nx,1,nt)

    d = 0.0
    delta = 0.0
    for j=1:nt,i=1:nx
      delta = disp_rel_1[i] - disp_rel_0[i]
      d = f[j] - disp_rel_0[i]
      if disp_rel_0[i]<f[j]<disp_rel_1[i]
        filter_array[i,1,j] = 0.5(1+cos(π*(d-delta/2)/(delta/2)))
      end
      if f_low<=f[j]<f_low+df
        filter_array[i,1,j] *= 0.5(1.+cos(π*(f[j]-(f_low+df))/df) )
      end
      if f[j] < f_low
        filter_array[i,1,j] = 0
      end
    end

    return filter_array
  end
end

function fmode_filter(nt::Int64,dt::Float64,nx::Int64,Lx::Float64)
  Poly::Array{Float64,1}=[1.1,1.9,-0.2]
  Polylow::Array{Float64,1} = [0.7,1.7,-0.2]

  f_low::Float64 = 1.1
  df::Float64 = 0.5

  @modefilter()
end


function p1mode_filter(nt::Int64,dt::Float64,nx::Int64,Lx::Float64)
  Poly::Array{Float64,1}=[1.4,3.0,-0.5]
  Polylow::Array{Float64,1} = [1.1,2.4,-0.3]

  f_low::Float64 = 1.6
  df::Float64 = 0.5

  @modefilter()
end

function p2mode_filter(nt::Int64,dt::Float64,nx::Int64,Lx::Float64)
  Poly::Array{Float64,1}=[1.6,3.8,-0.65]
  Polylow::Array{Float64,1} = [1.4,3.3,-0.62]

  f_low::Float64 = 1.6
  df::Float64 = 0.5

  @modefilter()
end

function p3mode_filter(nt::Int64,dt::Float64,nx::Int64,Lx::Float64)
  Poly::Array{Float64,1}=[2,4.1,-0.8]
  Polylow::Array{Float64,1}= [2,3.55,-0.7]

  f_low::Float64 = 1.6
  df::Float64 = 0.5

  @modefilter()
end

function p4mode_filter(nt::Int64,dt::Float64,nx::Int64,Lx::Float64)
  Poly::Array{Float64,1}=[1.41275585491,9.54809436998,-10.5814674886,7.99827826844,-2.42768573272]
  Polylow::Array{Float64,1} = [1.25437276419,8.13839040116,-7.73561854055,4.96643235694,-1.25914661289]

  f_low::Float64 = 1.6
  df::Float64 = 0.5

  @modefilter()
end

function p5mode_filter(nt::Int64,dt::Float64,nx::Int64,Lx::Float64)
  Poly::Array{Float64,1}=[2.35,5.6,-1.1]
  Polylow::Array{Float64,1} = [2.2,4.7,-1.0]

  f_low::Float64 = 1.6
  df::Float64 = 0.5

  @modefilter()
end

function p6mode_filter(nt::Int64,dt::Float64,nx::Int64,Lx::Float64)
  Poly::Array{Float64,1}=[1.985671159,8.0910898684,-3.2031633181]
  Polylow::Array{Float64,1} = [1.800354172,7.4293910566,-2.8459576439]

  f_low::Float64 = 1.6
  df::Float64 = 0.5

  @modefilter()
end

function p7mode_filter(nt::Int64,dt::Float64,nx::Int64,Lx::Float64)
  Poly::Array{Float64,1}=[1.6216558409,13.2697412028,-14.8197349148,9.44009291408,-2.17165011845]
  Polylow::Array{Float64,1} = [1.69587030179,10.4058096492,-9.1004586514,5.14732934967,-1.19380905259]

  f_low::Float64 = 1.6
  df::Float64 = 0.5

  @modefilter()
end

function fftfreq(nx::Int64,dx::Float64)
  freq_arr::Array{Float64,1} = zeros(nx)
  for i=1:div(nx,2)-1
    freq_arr[i+1]=i/(nx*dx)
    freq_arr[nx-i+1]=-i/(nx*dx)
  end
  freq_arr[div(nx,2)+1] = -1/(2*dx)
  return freq_arr
end

function fftfreq(nx::Int64)
  freq_arr::Array{Float64,1} = zeros(nx)
  for i=1:div(nx,2)-1
    freq_arr[i+1]=i/nx
    freq_arr[nx-i+1]=-i/nx
  end
  freq_arr[div(nx,2)+1] = -1/2
  return freq_arr
end

main()
