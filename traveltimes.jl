module traveltimes

using FFTW,Optim

filt(a,f) = real(ifft(fft(a)*f))

function fft_derivative(a;dt=1)
	nt = size(a,1)
	aω = fft(a,1)
	T = nt*dt
	ω = 2π/T * [0:div(nt,2)-1 ; -div(nt,2):-1]
	ω[div(nt,2)] = 0
	return real(ifft(im .* ω .* aω))
end

fft_integrate(a;dt=1) = sum(a,1) * dt

trapz(a;dt=1) = dt/2 * ( sum(a[1:end-1]) + sum(a[2:end]) )

function GB02(ξ0::AbstractArray{<:Real,1},ξ::AbstractArray{<:Real,1};dt=1,T_window_pix = 1:size(ξ0,1))
	# ξ(t)
	dt_ξ0 = fft_derivative(ξ0,dt=dt)
	f_t = zeros(size(ξ0,1))
	f_t[T_window_pix] .= 1

	num = trapz(f_t.*dt_ξ0.*(ξ0 - ξ),dt=dt)
	den = trapz(f_t.*dt_ξ0.^2,dt=dt)

	return num./den
end

GB02(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2},x_rec_pix::Integer;dt=1) = traveltimes_GB02(ξ0[:,x_rec_pix],ξ[:,x_rec_pix],dt=dt)

function GB02(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2},(x_min_pix,x_max_pix);dt=1)
	
	δτ_pix = zeros(x_max_pix - x_min_pix + 1)

	for (ind,x_pix) in x_min_pix:x_max_pix
		δτ_pix[ind] = GB02(ξ0,ξ,x_pix,dt=dt)
	end
	return δτ_pix 
end

function GB02(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2},x_ranges...;dt=1,dx=1)
	
	nx = size(ξ0,2)
	Lx = dx*nx

	x = range(-Lx/2,stop=Lx/2,length=nx+1)[1:nx]

	x_pix_ranges = Array{UnitRange{Int64},1}(undef,length(x_ranges))
	
	for (ind,x_range) in enumerate(x_ranges)
		x_low = x_range[1]
		x_low_pix = argmin(abs.(x .- x_low))
		x_high = x_range[end]
		x_high_pix = argmin(abs.(x .- x_high))
		x_pix_ranges[ind] = x_low_pix:x_high_pix
	end

	δτ_pix = zeros(sum(length.(x_pix_ranges)))

	ind = 0
	for x_pix_range in x_pix_ranges,x_pix in x_pix_range
	
		ind += 1
		δτ_pix[ind] = GB02(ξ0,ξ,x_pix,dt=dt)
	
	end
	return δτ_pix 
end

function GB02(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2},filter_fn,x_ranges...;dt=1,dx=1)
	ξ_filt = filt(ξ,filter_fn)
	ξ0_filt = filt(ξ0,filter_fn)
	GB02(ξ0_filt,ξ_filt,x_ranges,dt=dt,dx=dx)
end

function GB02(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2},filter_fn;dt=1)
	ξ_filt = filt(ξ,filter_fn)
	ξ0_filt = filt(ξ0,filter_fn)

	GB02(ξ0_filt,ξ_filt,dt=dt)
end

function GB02(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2},filter_fn,x_rec_pix;dt=1)
	ξ_filt = filt(ξ,filter_fn)
	ξ0_filt = filt(ξ0,filter_fn)

	GB02(ξ0_filt,ξ_filt,x_rec_pix,dt=dt)
end

function shift_to_later(y::AbstractArray{<:Real,1},Δ,dt=1)
	nt = size(y,1)
	T = nt*dt
	ω = 2π/T * (0:div(nt,2))
	irfft(rfft(y).*exp.(-im.*ω.*Δ),nt)
end

function GB04(ξ0::AbstractArray{<:Real,1},ξ::AbstractArray{<:Real,1};dt=1,T_window_pix = 1:size(ξ0,1))
	
	nt = size(ξ0,1)
	T = nt*dt

	τ = range(-T*(nt-1)/nt,stop=T*(nt-1)/nt,length=10000*nt)
	X = zeros(length(τ))
	
	f_t = zeros(size(ξ0,1))
	f_t[T_window_pix] .= 1

	X(Δτ) = trapz(f_t.*(ξ .- shift_to_later(ξ0,Δτ,dt)).^2,dt=dt)

	lower = [-T/2]
	upper = [T/2]
	initial_guess = [0.]
	inner_optimizer = GradientDescent()

	res= optimize(x->X(first(x)),lower,upper,initial_guess,Fminbox(inner_optimizer))
	Δτ = first(Optim.minimizer(res))
end

GB04(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2},x_rec_pix::Integer;dt=dt) = GB04(ξ0[:,x_rec_pix],ξ[:,x_rec_pix],dt=dt)

function GB04(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2};x::Real=0,dt=1,dx=1)
	nx = size(ξ0,2)
	Lx = dx*nx
	x_arr = LinRange(-Lx/2,Lx/2,nx+1)[1:nx]
	x_rec_pix = argmin(abs.(x_arr .- x ))
	GB04(ξ0[:,x_rec_pix],ξ[:,x_rec_pix],dt=dt)
end

function GB04(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2},(x_min_pix,x_max_pix);dt=1)
	
	δτ_pix = zeros(x_max_pix - x_min_pix + 1)

	for (ind,x_pix) in x_min_pix:x_max_pix
		δτ_pix[ind] = GB04(ξ0,ξ,x_pix,dt=dt)
	end
	return δτ_pix 
end


function GB04(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2},x_ranges...;dt=1,dx=1)
	
	nx = size(ξ0,2)
	Lx = dx*nx

	x = range(-Lx/2,stop=Lx/2,length=nx+1)[1:nx]

	x_pix_ranges = Array{UnitRange{Int64},1}(undef,length(x_ranges))
	
	for (ind,x_range) in enumerate(x_ranges)
		x_low = x_range[1]
		x_low_pix = argmin(abs.(x .- x_low))
		x_high = x_range[end]
		x_high_pix = argmin(abs.(x .- x_high))
		x_pix_ranges[ind] = x_low_pix:x_high_pix
	end

	δτ_pix = zeros(sum(length.(x_pix_ranges)))

	ind = 0
	for x_pix_range in x_pix_ranges,x_pix in x_pix_range
	
		ind += 1
		δτ_pix[ind] = GB04(ξ0,ξ,x_pix,dt=dt)
	
	end
	return δτ_pix 
end

function GB04(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2},filter_fn,x_ranges...;dt=1,dx=1)
	ξ_filt = filt(ξ,filter_fn)
	ξ0_filt = filt(ξ0,filter_fn)
	GB04(ξ0_filt,ξ_filt,x_ranges,dt=dt,dx=dx)
end

function GB04(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2},filter_fn;dt=1)
	ξ_filt = filt(ξ,filter_fn)
	ξ0_filt = filt(ξ0,filter_fn)

	GB04(ξ0_filt,ξ_filt,dt=dt)
end

function GB04(ξ0::AbstractArray{<:Real,2},ξ::AbstractArray{<:Real,2},filter_fn,x_rec_pix;dt=1)
	ξ_filt = filt(ξ,filter_fn)
	ξ0_filt = filt(ξ0,filter_fn)

	GB04(ξ0_filt,ξ_filt,x_rec_pix,dt=dt)
end

end