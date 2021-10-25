using Plots

V₀ = 1.0
R = 300
d = 400

function charged_circle(x_mid :: Int64, y_mid, radius :: Int64, val :: Float64, potential)
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

x = collect(-10:0.01:10)
y = collect(-10:0.01:10)
potential = zeros(length(x), length(y))
potential = charged_circle(length(x)÷2, length(y)÷2 - d, R, 1.0, potential)
fixed = iszero.(potential)
# potential = charged_circle(length(x)÷2, length(y)÷2 - d, R-25, 0.0, potential)

plot(potential, seriestype =:contourf, contours = 50, colorbar=true,
    xlabel = "X axis",
    ylabel = "Y axis",
    # xticks = 0:100:1000,
    # yticks = 0:100:1000,
    title = "Propagated Circle",
    size = (1440,1080),
    dpi = 200,
    c = :bluesreds
)


function VT(x,y, potential)
    grid = potential
    for (xind, xval) in enumerate(x)
        for (yind, yval) in enumerate(y)
            if fixed[xval, yval]
                r = sqrt(xval^2 + yval^2) - d
                grid[xind, yind] = ( log.(r / d) )
            end
        end
    end
    
    factor = V₀ / (log(R / d))
    grid = factor * grid
end

gridd = VT(x, y, potential)

plot(gridd, seriestype =:contourf, contours = 50, colorbar=true,
    xlabel = "X axis",
    ylabel = "Y axis",
    # xticks = 0:100:1000,
    # yticks = 0:100:1000,
    title = "Propagated Circle",
    size = (1440,1080),
    dpi = 200,
    c = :bluesreds
)
