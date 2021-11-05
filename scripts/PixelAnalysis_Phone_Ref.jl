using CSV
using Measurements
using Unitful
using UnitfulUS
using DataFrames
using StatsBase
using Statistics

cd("C:/Users/marcu/OneDrive/Desktop/PraktikumIII/QuadropoleIonTrap")

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
data_set1, errors1 = analyze_data_phone_screen(PATH1)
println("Mean & standard deviation of the first data set: $(mean(data_set1)) ± ($(std(data_set1))")
println("Relative error: $(std(data_set1) / mean(data_set1) * 100) %")
println("List of relative errors of the pure data: $(round.(errors1; digits=3))")

PATH2 = "data//exp_raw//PhoneScreen2_first2_calib_1_pixel_onPhone.txt";
data_set2, error2 = analyze_data_phone_screen(PATH2)
println("Mean of the second data set: $(mean(data_set2)) with standard deviation: $(std(data_set2))")

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
M4_average_velocity_T = mean(data_M4[1])

## Some quick gaussian boy action

using StatsPlots
M4_T_Cleaned = ustrip.(Measurements.value.(data_M4[1]))
min_val, max_val = extrema(M4_T_Cleaned);
# bins_M4 = range(min_val, max_val; length=4);
histogram_M4_T = histogram(ustrip.(Measurements.value.(M4_T_Cleaned)), bins = 3)

## Stokes law and things
# Stokes law - > calculating m/R
ρ = 1.2u"kg^3/m" # (at 1 ATM)
μ = 1.8e-5u"kg/m/s"
g_con = (9.80600 ± 0.000005)u"m/s^2" # https://www.metas.ch/metas/de/home/dok/gravitationszonen.html

m_R(velocity) = 6*π*velocity*μ / g_con

prop_values = m_R.(M4_average_velocity)

# Reynolds number
Re(v, R) = 2 * ρ * v * R / μ
