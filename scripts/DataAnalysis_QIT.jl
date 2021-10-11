using Measurements
using CSV
using DataFrames
using StatsPlots
using Distributions
using Unitful

## Read in data
data_path_qm_distr = "data\\exp_raw\\qM_Values_220V_AC.txt"
mass_charge_distr = CSV.File(data_path_qm_distr,
                        delim=' ',
                        ignorerepeated=false,
                        missingstring="NA") |> DataFrame
mass_charge_distr = mass_charge_distr |> Array

## Just a scatter plot
plot(mass_charge_distr, mass_charge_distr,
    xlabel="xdata",
    ylabel="ydata",
    seriestype=:scatter,
    markersize=3.25,
    dpi=300,
    size=((600, 400)),
    title=(nothing)
)

## Bar plot + fitted normal curve
histogram(mass_charge_distr, bins=range(minimum(mass_charge_distr), stop=maximum(mass_charge_distr), length=6),
    weights = ones(length(mass_charge_distr)) / length(mass_charge_distr),
    label = "Measured Data"
)
fitted_distribution = Distributions.fit(Distributions.Normal, Measurements.value.(mass_charge_distr))
μ_fitted = round(fitted_distribution.μ; digits=2)
σ_fitted = round(fitted_distribution.σ; digits=2)

plot!(fitted_distribution, xlabel = "Measured q/m values", ylabel = "Percentage of data set", label = "Gaussian Fit, \n μ=$(μ_fitted) & σ=$(σ_fitted)")

# savefig("qm_Distr.png")

## Calculate the average radius : diffraction experiment
λ = (655 ± 5)u"nm" # Wave-length of laser
R = (50 ± 0.4)u"cm" # Distance to screen
D = (3.8 ± 0.2)u"cm" # Distance between two minima 
Rₐ = D / 2 # Distance from center to one minimum

# Small angle approximation on equation : sin(θ) = 1.22 ⋅ λ / 2R
θ = Rₐ / R # small angle approximation to determine the sine of θ ; i.e. distance between minima / distance to screen
# Now solve equation for R̄ : R̄ = 1.22 ⋅ λ / 2θ
R̄ = (1.22 * λ / 2θ) |> u"μm"

rel_error = Measurements.uncertainty(R̄) / Measurements.value(R̄) * 100

## Calculate the average radius : phone screen


## Terminal velocity calculations : Setup 2 (Data set 2 & 4)
data_path_pixel_trails = "data\\exp_raw\\trail_distances.txt"
pixel_values = CSV.File(data_path_trails,
                delim=',',
                ignorerepeated=false,
                missingstring="NA"
) |> DataFrame |> Array # pixels
                    
pixels_per_2_mm = 100 # pixels / mm
Δs = pixel_values ./ pixels_per_2_mm

FPS = 30 # FPS
Δt = 1 // FPS # time under which we photograph the falling particle

v = Δs /Δt # calculate the velocity (almost terminal velocity, the particles have been falling for a relatively long time)

# Stokes law - > calculating m/R
ρ = 1.2u"kg^3/m" # (at 1 ATM)
μ = 1.8e-5u"kg/m/s"
g_con = (9.80600 ± 0.000005)u"m/s^2" # https://www.metas.ch/metas/de/home/dok/gravitationszonen.html

mR(velocity) = 6*π*velocity*μ / g_con

prop_values = mR.(v)

# Reynolds number 
Re(v, R) = 2 * ρ * v * R / μ
