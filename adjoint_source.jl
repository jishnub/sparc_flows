using PyCall
using FITSIO
using Polynomials
using NPZ
using StatsBase

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
  srcno::Int =(parse(Int,procno)+1)%8
  jobno::Int = parse(Int,procno)/8
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


  # acc::Array{Float64,3} = read(f[1])
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

  for (i,vel) in enumerate([0.44,0.6,0.75,0.9,1.2,1.4,1.7,1.9])
    spectral_ridges[i-1]["velocity"] = vel
  end


  # fmtstr = FormatExpr("{:3d} {:14.8f}$(" {:4d}"^5)\n")

  adj = zeros(nx,1,nt)
  misfit::Float64 = 0.0
  nmeasurements::Int = 0

  for modeno in 0:0
  # for modeno in keys(spectral_ridges)

    # Compute travel time difference
    if !isfile(joinpath(directory,"params.$modeno"))
      continue
    end

    mindist,maxdist,window = readdlm(joinpath(directory,"params.$modeno"))

    halftime = round(Int,window/(2.*dt_min))
    leng = 2*halftime+1

    vzcc_filtered::Array{Float64,3} = real(bfft(vzcc_fft .* spectral_ridges[modeno]["filter"]))
    data_filtered::Array{Float64,3} = real(bfft(data_fft .* spectral_ridges[modeno]["filter"]))
    ccdot::Array{Float64,3} = fft_time_deriv(vzcc_filtered,total_time_sec)
    # fitswrite((@adjoint_dir "ccdot_julia.fits"),ccdot)
    # return

    if !linesearch
      windowfile = open(joinpath(directory,"forward_src$(dec(srcno,2))_ls00","windows.$modeno"),"w")
    else
      windowfile = open(joinpath(directory,"forward_src$(dec(srcno,2))_ls00","windows.$(modeno)"),"r")
    end

    ttdiff_file = open((@forward_dir "ttdiff_julia.$modeno"),"w")

    for (xindex,dist_from_src) in enumerate(distances)
      if mindist<dist_from_src<maxdist

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


          vzcc_in_range = slice(vzcc_filtered,xindex,1,timest_index:timefin_index)

          loc = indmax(abs(vzcc_in_range))+timest_index-1
          lef = loc - halftime
          rig = loc + halftime

          write(windowfile,transpose([lef,rig]))
        else
          lef,rig = readdlm(IOBuffer(readline(windowfile)))
        end

        tau = compute_tt_gizonbirch(slice(vzcc_filtered,xindex,1,:),slice(data_filtered,xindex,1,:),lef,rig,dt_min)
        # tau = compute_tt_quadratic_fit(slice(vzcc_filtered,xindex,1,:),slice(data_filtered,xindex,1,:),lef,rig,dt_min)
        ttstr = @sprintf "%3d %14.8f %4d %4d %4d %4d %4d\n" xindex tau*60. lef rig loc timest_index timefin_index
        write(ttdiff_file,ttstr)

      end # if distance is in range

    end # loop over distances

    close(ttdiff_file)
    close(windowfile)

    # compute adjoint source
    tt_list = readdlm(@forward_dir "ttdiff_julia.$modeno")
    xindex_misfit = round(Int,tt_list[:,1])
    tau = tt_list[:,2]/60
    lef = round(Int,tt_list[:,3])
    rig = round(Int,tt_list[:,4])
    nmeasurements+=length(tau)
    iwls_misfit_factor = ones(size(tau))
    if iwls
      ttdiff_prev_path = joinpath(directory,"forward_src$(dec(srcno,2))_ls$(dec(jobno,2))","ttdiff_prev_iter.$modeno")
      if isfile(ttdiff_prev_path)
        ttdiff_prev_iter = readdlm(ttdiff_prev_path)[:,1]/60
        iwls_misfit_factor=(ttdiff_prev_iter.^2 +iwls_eps).^(iwls_pow/2-1)
      end
    end
    misfit += sum(tau.^2 .* iwls_misfit_factor)

    for (i,xindex) in enumerate(xindex_misfit)
      window_bandpass = zeros(nt)
      window_bandpass[lef[i]:rig[i]] = 1.0
      filtex = (slice(ccdot,xindex,1,:).*window_bandpass)[end:-1:1]*nx*nt

      con = -tau[i]/sum(slice(ccdot,xindex,1,lef[i]:rig[i]).^2*dt_min)
      adj[xindex,1,:] += reshape(filtex*con,1,1,nt)
      
    end

  end # loop over spectral ridges

  adj*=nt / 60

  if !linesearch
    fitswrite((@adjoint_dir "source_julia.fits"),adj)

    f=open((@adjoint_dir "information"),"w")
    write(f,"$nt\n")
    write(f,"$(t[1]*60)\n")
    write(f,"$(t[nt]*60)\n")
    write(f,"2D")
    close(f)
  end

  f=open(joinpath(directory,"kernel","misfit_$(dec(srcno,2))_$(dec(jobno,2))_julia"),"w")
  misfit *= 0.5
  println(@sprintf "%25.16E" misfit)
  write(f,@sprintf "%12d%12d%26.16E" srcno nmeasurements misfit)
  close(f)
end

function fitswrite(filename,array)
  f=FITS(filename,"w")
  write(f,array)
  close(f)
  println("Written file $(basename(filename))")
end

function fitsread(filename)
  f=FITS(filename,"r")
  arr=read(f[1])
  close(f)
  println("Reading $(basename(filename))")
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


# tic()
# # Useless function calls to compile code
# @time squeezeto1D(zeros(1,2));
# @time fftfreq(4);
# @time fftfreq(4,0.1);
# @time compute_tt_quadratic_fit(float(collect(1:4)),float(collect(1:4)),1,3,0.5)
# @time compute_tt_gizonbirch(float(collect(1:4)),float(collect(1:4)),1,3,0.5)
# @time fft_time_deriv(zeros(2),10.);
# @time fft_time_deriv(zeros(2,2,2),10.);
# @time fmode_filter(4,0.5,4,100.);
# toc()

@time main()

# @profile main()
# # Profile.print()
# Profile.print()
