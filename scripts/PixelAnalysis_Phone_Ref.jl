using CSV
using Measurements
using Unitful
using UnitfulUS
using DataFrames
using StatsBase

function analyze_data(path; PPI = 458, pixels_calibration = 1)
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

    return (values .* conversion_factor, average_calib)
end;


PATH1 = "data//exp_raw//PhoneScreen1_first2_calib_1_pixel_onPhone.txt";
PATH2 = "data//exp_raw//PhoneScreen2_first2_calib_1_pixel_onPhone.txt";

data_set1, calib1 = analyze_data(PATH1)
println("Mean of the first data set: $(mean(data_set1))")

data_set2 = 0

data_set_same_img = 0
