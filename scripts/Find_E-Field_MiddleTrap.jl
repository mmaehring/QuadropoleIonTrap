using DrWatson
@quickactivate "QuadropoleIonTrap"
using Plots
using LinearAlgebra
using ProgressBars

function charged_circle(x_mid :: Int64, y_mid :: Int64, radius :: Int64, val :: Float64, potential)
    xrange = collect(x_mid-radius:x_mid+radius)[2:end-1]
    # calculates y corresponding to an x if the circle is centered at (w, h) w. radius radius
    δy(x) = sqrt( radius^2 - (x - x_mid)^2 )
    for i in xrange
        Δy = convert.(Int64, round.( δy(i)) )
        yrange = y_mid-Δy:y_mid+Δy
        potential[yrange, i] .= val
    end
    return potential
end

function generate_potential(w_sample_bc = false)
    potential = zeros((1000, 1000))
    if w_sample_bc
        potential[1, :] = lower_y
        potential[end, :] = upper_y
        potential[:, 1] = lower_x
        potential[:, end] = upper_x
    end
    return potential
end


function compute_potential_fixed(potential, fixed_indices, n_iter)
    len = size(potential)[1]
    for n in ProgressBar(1:n_iter)
        for i in 2:len-1
            for j in 2:len-1
                if (fixed_indices[j, i])
                    tmp = 1/4 * ( potential[j+1, i] + potential[j-1, i]  + potential[j, i+1] + potential[j, i-1] )
                    potential[j, i] = tmp
                end
            end
        end
    end
    return potential
end

V₀ = 1.0
potential = generate_potential(false)
potential = charged_circle(500, 800, 100, -V₀, potential)
potential = charged_circle(500, 200, 100, V₀, potential)
potential = charged_circle(200, 500, 100, V₀, potential)
potential = charged_circle(800, 500, 100, -V₀, potential)

fixed = iszero.(potential)
potential = charged_circle(200, 500, 100, 0.0, potential)
potential = charged_circle(800, 500, 100, 0.0, potential)

plot(potential, seriestype =:contourf, contours = 50, colorbar=true,
    xlabel = "X axis",
    ylabel = "Y axis",
    xticks = 0:100:1000,
    yticks = 0:100:1000,
    title = "Charged Circle",
    c = :bluesreds
)

# @view fixed[700:900, 400:600]

potential_fixed = compute_potential_fixed(potential, fixed, 1e4)

plot(potential_fixed, seriestype =:contourf, contours = 20, colorbar=true,
    xlabel = "X axis",
    ylabel = "Y axis",
    xticks = 0:100:1000,
    yticks = 0:100:1000,
    title = "Propagated Circle",
    size = (1440,1080),
    dpi = 200,
    c = :bluesreds,
    colorbar_title = "Potential"
)

##

E_mid_val0 = potential_fixed[500, 500] - potential_fixed[501, 500] # N
E_mid_val1 = potential_fixed[500, 500] - potential_fixed[501, 501] # NE
E_mid_val2 = potential_fixed[500, 500] - potential_fixed[500, 501] # E
E_mid_val3 = potential_fixed[500, 500] - potential_fixed[499, 501] # SE
E_mid_val4 = potential_fixed[500, 500] - potential_fixed[499, 500] # S
E_mid_val5 = potential_fixed[500, 500] - potential_fixed[499, 499] # SW
E_mid_val6 = potential_fixed[500, 500] - potential_fixed[500, 499] # W
E_mid_val7 = potential_fixed[500, 500] - potential_fixed[501, 499] # NW

# we assume from 0 to 1e3 -> 14 mm => 1 step = 14 / 1000 mm => 1 step = 14*10^-6 m

E_vals = [E_mid_val0, E_mid_val1, E_mid_val2, E_mid_val3, E_mid_val4, E_mid_val5,
          E_mid_val6, E_mid_val7] ./ 14e-6

E_vals_z = [E_mid_val0, E_mid_val4] ./ 14e-6

# E = α * V / d
V = 1u"V"
E_val_final_all = mean(abs.(E_vals))u"V/m"
E_val_final_ud = mean(abs.(E_vals_z))u"V/m"
α =  E_val_final_ud * (d|>u"m") / V
