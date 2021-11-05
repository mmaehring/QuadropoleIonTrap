using Images, TestImages, Plots, ImageView, ImageDraw

imgt = testimage("blobs");
mosaic(imgt);
img_gray = convert(Array{Float64,2},Gray.(imgt));

blobs_log = blob_LoG(img_gray, [10]);

heatmap(img_gray, color=:grays, aspect_ratio=1)

for blobs in blobs_log
    coord = blobs.location
    draw!(imgt, Ellipse(CirclePointRadius(coord[2], coord[1], 4; thickness = 2, fill = false)),
          RGB{N0f8}(0,1,0))
    # scatter!(coord[2],coord[1], alpha=1, facecolors="none", edgecolors="r")
end
mosaic(imgt)

 # s=blobs.amplitude*5000

img_path = "C:/Users/marcu/OneDrive/Desktop/PraktikumIII/QuadropoleIonTrap/data/exp_raw/Power-Screen/Powder2.png"
img = load(img_path)
img_gray = convert(Array{Float64,2}, Gray.(img));

blobs_log = blob_LoG(img_gray, [5]);

for blobs in blobs_log
    coord = blobs.location
    draw!(img, Ellipse(CirclePointRadius(coord[2], coord[1], 3; thickness = 1, fill = false)))
    # scatter!(coord[2],coord[1], alpha=1, facecolors="none", edgecolors="r")
end
mosaic(img)
