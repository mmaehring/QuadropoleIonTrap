using DrWatson
# @quickactivate("")
begin
    using Measurements
    using CSV
    using DataFrames
    using StatsPlots
    using Distributions
    using Unitful
    using UnitfulUS
    using DataFrames
    using StatsBase
    using Statistics
end

# @quickactivate("QuadropoleIonTrap")
cd("C:\\Users\\marcu\\OneDrive\\Desktop\\PraktikumIII\\QuadropoleIonTrap")
## Read in data
begin
    data_path_qm_distr = "data\\exp_raw\\qM_Values_220V_AC.txt"
    mass_charge_distr_Voltages_DC = CSV.File(data_path_qm_distr,
    delim=' ',
    ignorerepeated=false,
    missingstring="NA") |> DataFrame
    mass_charge_distr_Voltages_DC = mass_charge_distr_Voltages_DC |> Array

    # Old data -> for later use
    data_path_qm_distr_OLD = "data\\exp_raw\\qM_Values_220V_AC_OLD.txt"
    mass_charge_distr_Voltages_DC_OLD = CSV.File(data_path_qm_distr_OLD,
    delim=' ',
    ignorerepeated=false,
    missingstring="NA") |> DataFrame
    mass_charge_distr_Voltages_DC_OLD = mass_charge_distr_Voltages_DC_OLD |> Array
end

## Just a scatter plot
# plot(mass_charge_distr_Voltages_DC, mass_charge_distr_Voltages_DC,
#     xlabel="xdata",
#     ylabel="ydata",
#     seriestype=:scatter,
#     markersize=3.25,
#     dpi=300,
#     size=((600, 400)),
#     title=(nothing)
# )

## Bar plot + fitted normal curve (new data)
begin
    histogram(mass_charge_distr_Voltages_DC,
    bins= 7,
    weights = ones(length(mass_charge_distr_Voltages_DC)) / length(mass_charge_distr_Voltages_DC),
    label = "Measured Data \n",
    xlims = (7, 38),
    legend = :topleft,
    dpi = 500
    )
    fitted_distribution = Distributions.fit(Distributions.Normal, Measurements.value.(mass_charge_distr_Voltages_DC))
    μ_fitted = round(fitted_distribution.μ; digits=2)
    σ_fitted = round(fitted_distribution.σ; digits=2)

    plot!(fitted_distribution, xlabel = "DC Values [V]",
    ylabel = "Percentage of data set",
    label = "Gaussian Fit: μ=$(μ_fitted) & σ=$(σ_fitted)")

    # savefig("plots/qm_distr_new_data.pdf")
end

## Old data and new data
# begin
#     old_and_new_charge_distr = vcat(mass_charge_distr_Voltages_DC_OLD |> Array, mass_charge_distr_Voltages_DC |> Array)
#
#     histogram(old_and_new_charge_distr,
#     bins= 7,
#     weights = ones(length(old_and_new_charge_distr)) / length(old_and_new_charge_distr),
#     label = "Measured Data \n",
#     title = "Old and new",
#     xlims = (10, 35),
#     legend=:topleft,
#     yticks = ([0, 0.1, 0.2, 0.3, 0.4, 0.5], string.([0, 10, 20, 30, 40, 50])),
#     dpi = 250
#     )
#     fitted_distribution = Distributions.fit(Distributions.Normal, Measurements.value.(old_and_new_charge_distr))
#     μ_fitted = round(fitted_distribution.μ; digits=2)
#     σ_fitted = round(fitted_distribution.σ; digits=2)
#
#     plot!(fitted_distribution, xlabel = "DC Values [V]", ylabel = "Percentage of data set", label = "\n Gaussian Fit: μ=$(μ_fitted) & σ=$(σ_fitted)")
# end


## NEW EXPERIMENT PART
## Calculate the average radius : diffraction experiment
begin
    const λ = (655 ± 5)u"nm" # Wave-length of laser
    const R = (50 ± 0.4)u"cm" # Distance to screen
    const D = (3.8 ± 0.2)u"cm" # Distance between two minima
    const Rₐ = D / 2 # Distance from center to one minimum
end

# Small angle approximation on equation : sin(θ) = 1.22 ⋅ λ / 2R
# small angle approximation to determine the sine of θ ; i.e. distance between minima / distance to screen
θ = Rₐ / R
# Now solve equation for R̄ : R̄ = 1.22 ⋅ λ / 2θ
R̄ = (1.22 * λ / (2θ)) |> u"μm"

rel_error = Measurements.uncertainty(R̄) / Measurements.value(R̄) * 100

## Calculate the average radius : phone screen
function analyze_data_phone_screen(PATH; PPI = 458, pixels_calibration = 1)
    df = CSV.read(PATH, DataFrame)

    # Get the first to values from the array, the calibration values for the OPhone screen pixels
    calib = df[1:2, 1:end] |> Array

    # This gives the amount of pixels on the image / 1 pixel on the iPhone's screen
    calibration_values = [hypot(i,j)±k for (i,j,k) in zip(calib[:, 1], calib[:, 2], calib[:, 3])] ./ pixels_calibration
    average_calib = mean(calibration_values)

    # Conversion factor from pixels on Picamera image to inches and then to μm
    conversion_factor = 1/(PPI * average_calib)  * 1u"sinch_us" |> u"μm"

    # get the rest of the raw data, note that these are the DIAMETER values
    data = df[3:end, begin:end] |> Array
    values = [hypot(i,j)±k for (i,j,k) in zip(data[:, 1], data[:, 2], data[:, 3])] ./ 2 # get radii

     measurements = values .* conversion_factor
     relative_errors = Measurements.uncertainty.(measurements) ./ Measurements.value.(measurements) * 100
     return (measurements, relative_errors)
end;

PATH1 = "data//exp_raw//PhoneScreen1_first2_calib_1_pixel_onPhone.txt";
data_set1, errors1 = analyze_data_phone_screen(PATH1);
println("Mean & standard deviation of the first data set: $(mean(data_set1)) ± ($(std(data_set1))")
mean_Radius_Phone1 = mean(data_set1)
println("Relative error: $(std(data_set1) / mean(data_set1) * 100) %")
println("List of relative errors of the pure data: $(round.(errors1; digits=3))")

PATH2 = "data//exp_raw//PhoneScreen2_first2_calib_1_pixel_onPhone.txt";
data_set2, error2 = analyze_data_phone_screen(PATH2)
println("Mean & standard deviation of the second data set: $(mean(data_set2)) ± ($(std(data_set2)))")
mean_Radius_Phone2 = mean(data_set2)

function analyze_data_terminal_velocity(calibration_val, T, PATH)
    #=
    calibration_val is the average diameter of the pendulum in pixels: Pixel -> 2 mm
    T is the time interval for one flash / photo
    path is the path to the data
    =#
    df = CSV.read(PATH, DataFrame)

    # Get only the good data
    df = subset(df, :OkQuality => x -> x .== "G")

    # Data where the TerminalVelocity is ascertained
    df_T = subset(df, :SureTerminal => x -> x .== 1)
    # Data where the TerminalVelocity is NOT ascertained
    df_NT = subset(df, :SureTerminal => x -> x .== 0)

    println(df_T)

    # get the raw data, note that these are lengths of one bright stripe
    data_T = mean(df_T[!, [:y1, :y2]] |> Array, dims=2)
    data_NT = mean(df_NT[!, [:y1, :y2]] |> Array, dims=2)
    data_ALL = mean(df[!, [:y1, :y2]] |> Array, dims=2)

    # Conversion factor from pixels to mm    P ⋅ (P / 2mm)^-1 = 2mm * P / P
    data_T /= calibration_val / 2u"mm"
    data_NT /= calibration_val / 2u"mm"
    data_ALL /= (calibration_val / 2u"mm")

    # Convert to velocities
    data_T /= T
    data_NT /= T
    data_ALL /= T

    # New data array
    # return DataFrame("T" => data_T, "NT" => data_NT, "ALL" => data_ALL);
    return [data_T .|> u"m/s", data_NT .|> u"m/s", data_ALL .|> u"m/s"];
end;

path_M4 = "data/exp_raw/TerminalVelocity/RawVids/Particles/M4/M4.csv"
conversion_setup_2 = mean([95 ± 2, 90 ± 2, 89 ± 1.5]) # pixel / 2mm
f = 50u"Hz"

data_M4 = analyze_data_terminal_velocity(conversion_setup_2, 1/f |> u"s", path_M4)
M4_average_velocity_T, M4_std_velocity_T = [mean(data_M4[1]), std(data_M4[1])]
M4_average_velocity_T_mm, M4_std_velocity_T_mm = [mean(data_M4[1]), std(data_M4[1])] .|> u"mm/s"
## Some quick gaussian boy action
M4_T_Cleaned = ustrip.(Measurements.value.(data_M4[1]))
min_val, max_val = extrema(M4_T_Cleaned);
# bins_M4 = range(min_val, max_val; length=4);
histogram_M4_T = histogram(ustrip.(Measurements.value.(M4_T_Cleaned)), bins = 3)

## Stokes law and things
# Stokes law - > calculating m/R
const ρ = (1.2±0.05)u"kg*m^-3" # (at 1 ATM)
const μ = 10^-5 * (1.8±0.05)u"kg*m^-1*s^-1"
const g_con = (9.80600 ± 0.000005)u"m/s^2" # https://www.metas.ch/metas/de/home/dok/gravitationszonen.html

m_R(velocity) = 6*π*velocity*μ / g_con

mass_radius_quotient = m_R.(M4_average_velocity_T)

## Mass calculation for different R:
m_Laser = mass_radius_quotient * R̄ |> u"kg"
m_Phone1 = mass_radius_quotient * mean_Radius_Phone1 |> u"kg"
m_Phone2 = mass_radius_quotient * mean_Radius_Phone2 |> u"kg"

masses = [m_Laser, m_Phone1, m_Phone2]

# Reynolds number for different R
Re(v, R_val) = @. 2 * ρ * v * R_val * 1/μ

# Laser [1], Phone1 [2] and Phone2 [3]
radii = [R̄, mean_Radius_Phone1, mean_Radius_Phone2] .|> u"m";
# reynolds_numbers = ustrip.(Re(M4_average_velocity_T, radii))
reynolds_numbers = (Re(M4_average_velocity_T, radii))
maximal_reynolds = Measurements.value.(reynolds_numbers) .+ Measurements.uncertainty.(reynolds_numbers)

## Charge calculation
d = (4 ± 0.25)u"mm"
r = (3 ± 0.25)u"mm" # 2r
r = r/2

α_measured_from_plot = 2.5 ± 0.41
α_measured_from_plot = 1.40 ± 0.21
# α_measured_from_plot = 2.20 ± 0.21
# α_fem = 0.0003? -> Failure

Ec(V) = α_measured_from_plot * V / (d|>u"m")
electric_field_vals = Ec.(mass_charge_distr_Voltages_DC_units)

g = (9.80600 ± 0.000005)u"m/s^2"

q_m_vals = g./electric_field_vals .|> u"C/g"
q_m_vals_mu = g./electric_field_vals .|> u"μC/g"

# using DelimitedFiles
# writedlm("data\\exp_pro\\q_m_vals_mu_alpha_1.4.txt", q_m_vals_mu)


q(m, V) = @. m * g  / Ec(V)

mass_ill_use = Measurements.value(mean(masses)) ± Measurements.uncertainty(std(masses))
mass_charge_distr_Voltages_DC_units = mass_charge_distr_Voltages_DC .*1u"V"

charges = (q(mass_ill_use, mass_charge_distr_Voltages_DC_units)) .|>u"C"

using PhysicalConstants.CODATA2018: e
charges_in_unit_charges = (charges / ustrip(e)) .|> u"kC"
charges_in_1k_unit_charges = (charges ./ 1000e) # charge given in 1000s of unit charges

average_charge_in_1ks_of_unit_charges = mean(charges / e) / 1000

## Printing in latex table format
begin
    for (i, e1, e2) in zip(1:15, ustrip.(charges .|> u"fC"), charges_in_1k_unit_charges)
        println(i, raw" & $", e1, raw"$ & $", e2, raw"$ Å")
        println(raw"\hline")
    end
    println(raw"Mean & $", mean(ustrip.(charges .|> u"fC")), raw"$ & $", mean(charges_in_1k_unit_charges), raw"$ Å")
end



## Histogram fitting only new charge / mass data
begin
    histogram(Measurements.value.(charges_in_1k_unit_charges),
        # bins = 4,
        bins=range(minimum(Measurements.value.(charges_in_1k_unit_charges)), stop=maximum(Measurements.value.(charges_in_1k_unit_charges)), length=6),
        weights = ones(length(charges_in_1k_unit_charges)) / (10*length(charges_in_1k_unit_charges)),
        label = "1000's of unit charges \n",
        ylims = (0, 0.055), legend=:topright,
        # xlims = (145, 455),
        yticks =([0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07], string.([0, 10, 20, 30, 40, 50 , 60, 70])),
        xlabel = "Amount of unit charges"
    )

    fitted_distribution = Distributions.fit(Distributions.Normal, Measurements.value.(charges_in_1k_unit_charges))
    μ_fitted = round(fitted_distribution.μ; digits=2)
    σ_fitted = round(fitted_distribution.σ; digits=2)
    fc(x) = 12.5 * pdf(fitted_distribution, x)
    plot!(0:1000, fc.(0:1000), ylabel = "Percentage of data set", label = "Fit: μ=$(μ_fitted), σ=$(σ_fitted)")
    savefig("plots/q_distr_newdata_alpha=1_4.pdf")
end
# begin
#     histogram(Measurements.value.(charges_in_1k_unit_charges),
#         bins = 5,
#         # bins=range(minimum(Measurements.value.(charges_in_1k_unit_charges)), stop=maximum(Measurements.value.(charges_in_1k_unit_charges)), length=6),
#         weights = ones(length(charges_in_1k_unit_charges)) / (10*length(charges_in_1k_unit_charges)),
#         label = "Charge distribution in 1000's of unit charges \n",
#         # ylims = (0, 0.055),
#         # xlims = (145, 455),
#         # yticks =([0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07], string.([0, 10, 20, 30, 40, 50 , 60, 70])),
#         xlabel = "Amount of unit charges",
#         dpi = 500
#     )
#
#     fitted_distribution = Distributions.fit(Distributions.Normal, Measurements.value.(charges_in_1k_unit_charges))
#     μ_fitted = round(fitted_distribution.μ; digits=2)
#     σ_fitted = round(fitted_distribution.σ; digits=2)
#     plot!(fitted_distribution, ylabel = "Percentage of data set", label = "Gaussian Fit: μ=$(μ_fitted) & σ=$(σ_fitted)")
#     # plot(fitted_distribution, seriestype=:line)
# end
