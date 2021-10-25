using DrWatson
@quickactivate "QuadropoleIonTrap"
using Plots
using LinearAlgebra
using ProgressBars
# using DifferentialEquations

## Solve Laplace : ∇² V = 0


##
# function charged_square(x :: Int64, y :: Int64, radius :: Int64, val :: Int64, potential)
#     radius = convert(Int64, radius / 2)
#     xrange = x-radius:x+radius
#     yrange = y-radius:y+radius
#     # @assert (1 < x radius < 1000) && (1 < y+radius < 1000)
#     potential[xrange, yrange] .= val
#     return potential
# end
#
# function charged_rectangle(x :: Int64, y :: Int64, radius :: Int64, val :: Int64, potential)
#     xrange = collect(x-radius:x+radius)[2:end-1]
#     array = collect(1:radius)
#     app = (radius .- collect(1:radius))[1:end-1]
#     append!(array, app)
#     for (xpos,ypos) in zip(xrange, array)
#         potential[xpos, y-ypos:1:y+ypos] .= 10
#     end
#     return potential
# end


## Boundary conditions
edge = collect(range(-1, 1, length = 1000))
upper_y = cos.(π .* edge ./ 2)
lower_y = edge.^4
upper_x = 1 ./ (exp(-1) .- exp(1)) * (exp.(edge) .- exp(1))
lower_x = 0.5 .* (edge.^2 - edge)

plot(edge, lower_x)

## Meshgrid
xv = edge' .* ones(1000)
yv = ones(1000)' .* edge

## Potential Solve
function compute_potential(potential, n_iter)
    len = size(potential)[1]
    for n in ProgressBar(1:n_iter)
        for i in 2:len-1
            for j in 2:len-1
                tmp = 1/4 * ( potential[j+1, i] + potential[j-1, i]  + potential[j, i+1] + potential[j, i-1] )
                potential[j, i] = tmp
            end
        end
    end
    return potential
end

## Solve for potential
function generate_potential(w_sample_bc = true)
    potential = zeros((1000, 1000))
    if w_sample_bc
        potential[1, :] = lower_y
        potential[end, :] = upper_y
        potential[:, 1] = lower_x
        potential[:, end] = upper_x
    end
    return potential
end

# Iterate through array using discretized laplacian
potential = generate_potential(true)
potential = compute_potential(potential, 1e4)

## Contour Plot
# Before
plot(generate_potential(true), seriestype = :contourf, contours = 30, colorbar=true)
# After
plot(potential , seriestype = :contourf, contours = 30, colorbar=true)

## Block of fixed potential (not analytically possible)
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


## Compute potential with fixed zone of charge

base_potential = generate_potential(false)
plot(base_potential, seriestype =:contourf, contours = 50, colorbar=true,
    xlabel = "X axis",
    ylabel = "Y axis",
    title = "Raw Potential"
)

potential = charged_circle(500, 800, 150, 10.0, base_potential)
plot(potential, seriestype =:contourf, contours = 50, colorbar=true,
    xlabel = "X axis",
    ylabel = "Y axis",
    xticks = 0:100:1000,
    yticks = 0:100:1000,
    title = "Charged Circle"
)

fixed = iszero.(potential)

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

# potential_no_fixed = compute_potential(potential, 1e4)
potential_fixed = compute_potential_fixed(potential, fixed, 1e4)

plot(potential_fixed, seriestype =:contourf, contours = 50, colorbar=true,
    xlabel = "X axis",
    ylabel = "Y axis",
    xticks = 0:100:1000,
    yticks = 0:100:1000,
    title = "Propagated Circle"
)

# using DelimitedFiles
# writedlm("Laplace_1_Circ_Pot_100IT.csv", potential_fixed, ", ")


## Geometry of our setup
potential = generate_potential(false)
potential = charged_circle(500, 800, 100, -1.0, potential)
potential = charged_circle(500, 200, 100, 1.0, potential)
potential = charged_circle(200, 500, 100, 1.0, potential)
potential = charged_circle(800, 500, 100, -1.0, potential)

fixed = iszero.(potential)
potential = charged_circle(200, 500, 100, 0.0, potential)
potential = charged_circle(800, 500, 100, 0.0, potential)

plot(potential, seriestype =:contourf, contours = 50, colorbar=true,
    xlabel = "X axis",
    ylabel = "Y axis",
    xticks = 0:100:1000,
    yticks = 0:100:1000,
    title = "Charged Circle",
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
    c = :bluesreds
)

using Interpolations
using ForwardDiff
