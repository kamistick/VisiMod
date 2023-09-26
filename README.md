VisiMod: Map Modeled Visibility Index Across Wildland Landscapes
================
Katherine Mistick
2023-09-26

- <a href="#1-intoduction" id="toc-1-intoduction"><strong>1
  Intoduction</strong></a>
- <a href="#2-download-lidar-data-to-prepare-for-visimod-workflow"
  id="toc-2-download-lidar-data-to-prepare-for-visimod-workflow"><strong>2
  Download lidar data to prepare for VisiMod workflow</strong></a>
- <a href="#3-the-visimod-workflow"
  id="toc-3-the-visimod-workflow"><strong>3 The VisiMod
  Workflow</strong></a>
- <a href="#4-using-visimod-to-generate-an-rgb-vi-map"
  id="toc-4-using-visimod-to-generate-an-rgb-vi-map"><strong>4 Using
  VisiMod() to generate an RGB VI map</strong></a>

------------------------------------------------------------------------

## **1 Intoduction**

The VisiMod package was built to facilitate visibility analyses. It
provides a suite of tools that allows users to build predictive models
that can map visibility (in terms of visibility index (VI) where VI =
(visible area)/(total area)) across entire landscapes. [Mistick et.
al. 2023 (Using airborne lidar and machine learning to predict
visibility across diverse vegetation and terrain
conditions)](https://content.csbs.utah.edu/~pdennison/reprints/denn/2023_Mistick_etal_IJGIS.pdf)
is the original research that led to the development of VisiMod.

Instead of calculating VI by running viewsheds at every pixel in a
rasterized study area, VisiMod allows users to build predictive models
based on a sample of points (at which viewsheds are run using
[terra::viewshed()](https://rdrr.io/cran/terra/man/viewshed.html)) to
predict VI at every pixel.

VisiMod is designed to be used in wildland (i.e. undeveloped) areas. It
relies on a digital surface model (DSM) that is inclusive of vegetation.
If VisiMod were deployed for developed areas any aboveground
obstructions would be treated as vegetation. This is not recommended.

The VisiMod package can be installed by calling
`devtools::install_github("kamistick/VisiMod")`in the RStudio console.
To run the following examples the following packages must be loaded:
`VisiMod`, `terra`, `httr2`, `httr`, `lidR`, `future`, `ggplot2`

## **2 Download lidar data to prepare for VisiMod workflow**

------------------------------------------------------------------------

VisiMod only requires users to provide a digital terrain model (DTM) and
digital surface model (DSM) to run the full workflow. While fine-scale
resolution DTMs are available to download (ex.), fine-scale DSMs most
typically are derived directly from lidar. The example below walks a
user through downloading and processing four USGS 3DEP lidar tiles for
use in the VisiMod workflow.

We’ll use The National Map (TNM) API, provided by the USGS, to access
3DEP lidar data. For more information on how to use the THM API, [click
here](https://apps.nationalmap.gov/tnmaccess/#/). To interact with the
TNM GUI and download data, [click
here](https://apps.nationalmap.gov/downloader/).

For this example we’re going to look at a 4km<sup>2</sup> area in
southern Utah, USA:

![image](images/leaflet.png) \*The above figure is a static image of an
interactive leaflet map that does not render in this GitHub README.md

### **2.1 Set Up Directories**

First, we will set up our directories. In this example the main
directory `lidar_dir` is a folder with four subfolders: `noise`, `dtms`,
`dsms`, and `out_files`. We will download, prepare, and generate new
files in these directories throughout these examples. The code below
will check if the subfolders in your directory exist, and create them if
they do not.

``` r
lidar_dir <- "/path/to/your/Data"
```

``` r
for (folder in c("noise", "dtms", "dsms", "out_files")){
  if (!dir.exists(file.path(lidar_dir, folder))){
  dir.create(file.path(lidar_dir, folder))
  } else {
      print("Directory already exists.")
  }
}
```

    ## [1] "Directory already exists."
    ## [1] "Directory already exists."
    ## [1] "Directory already exists."
    ## [1] "Directory already exists."

### **2.2 Download lidar using The National Map API**

Next we will search for available lidar products within our bounding
box:

``` r
bounding_box <- "-113.796282,38.008962,-113.785245,38.018234"
product_format <- "LAS,LAZ"
url <- paste0("https://tnmaccess.nationalmap.gov/api/v1/products?bbox=", bounding_box, "&prodFormats=", product_format)
req <- request(url)
resp <- req %>% req_perform()
```

Next we’ll grab the download links for the LAZ point clouds, of which
there are 4. We’ll also need the filenames so we can save each file.

``` r
for (x in 1:length(resp_body_json(resp)[2]$items)){
  # get download link
  download_link <- resp_body_json(resp)[2]$items[[x]]$urls$LAZ
  
  #get file name without white spaces
  filename <- gsub(" ", "", resp_body_json(resp)[2]$items[[x]]$title)
  
  #use httr::GET() to download each laz tile
  print(paste0("Downloading ", filename))
  httr::GET(download_link, write_disk(paste0(lidar_dir, "/", filename, ".laz"), overwrite=TRUE))
}
```

### **2.3 Process lidar with lidR in R**

Next we’ll use the [lidR package](https://r-lidar.github.io/lidRbook/)
to process our tiles into a DTM and DSM. First we will read the files in
as a lascatalog, and plot them. It is recommended that you index your
lidar tiles before proceeding using `lidR:::catalog_laxindex(cat)`.
Indexing allows for faster processing and is especially advantageous
when dealing with large LAScatalogs. For more on LAScatalogs check on
this vignette: [LAScatalog formal
class](https://cran.r-project.org/web/packages/lidR/vignettes/lidR-LAScatalog-class.html)

``` r
# read in files as a LAScatalog
cat <- readLAScatalog(lidar_dir) 
lidR:::catalog_laxindex(cat)
```

``` r
plot(cat)
```

![](README_files/figure-gfm/plot-las-cat-1.png)<!-- -->

#### **2.3.1 Remove Noise**

First we’ll clean up the lidar data to remove any noise. This step can
look very different depending on your lidar data. Some datasets are
fairly clean to begin with, others may have extremely high or low noise
flagged as classification 7 or 18. These points are removed by
filtering. However, there may be certain high or low points that were
not classified as noise, but should be. We run an additional noise
filter and remove the noise flagged points a second time to ensure all
extraneous points have been removed. If you have multiple cores
available to work on, we suggest using them by implementing
`plan(multisession, workers = 4L)` and `set_lidr_threads(4L)`, which
would implement 4 cores for parallel processing.

``` r
# set up parallel processing, if available
plan(multisession, workers = 4L)
set_lidr_threads(4L)

# define directoy where clean laz files will be written
clean_dir <- file.path(lidar_dir, "noise")

# filter points already classified as noise or withheld
opt_filter(cat) <- "-drop_class 7 18 -drop_withheld"

# perform an additional noise filter
opt_output_files(cat) <- file.path(clean_dir, "{ORIGINALFILENAME}")
opt_laz_compression(cat) <- TRUE
classify_noise(cat, sor())

# read back in cleaned catalog and filter points classified as noise or withheld
cat_clean <- readLAScatalog(clean_dir)
opt_filter(cat_clean) <- "-drop_class 7 18 -drop_withheld"
```

#### **2.3.2 Generate Rasters**

Next we’ll generate our 1 meter DTM and DSM using `rasterize_terrain()`
and `rasterize_canopy()`. This step will genereate and save a raster (in
TIF format) in the `lidar_dir` directory (in this example they are saved
in subfolders `dtms` and `dsms`).

``` r
# new directories 
dtm_dir <- file.path(lidar_dir, "dtms")
dsm_dir <- file.path(lidar_dir, "dsms")

# generate dtms
opt_output_files(cat_clean) <- file.path(dtm_dir, "dtm_{XLEFT}_{YBOTTOM}")
rasterize_terrain(cat_clean, res= 1, algorithm = tin(), overwrite = TRUE)

# generate dsms
opt_output_files(cat_clean) <- file.path(dsm_dir, "dsm_{XLEFT}_{YBOTTOM}")
rasterize_canopy(cat_clean, res = 1, dsmtin(), overwrite= TRUE)
```

Next we have to mosaic all four DTM tiles together to get a continuous
DTM for the study area. First we’ll read in each individual DTM raster
using `rast()`, then we’ll mosaic them together using `mosaic()`. If
you’d like to save your study-area-wide DTM use `writeRaster()`.

``` r
dtm_dir <- file.path(lidar_dir, "dtms")
dtm_tiles <- list.files(dtm_dir, pattern = "*.tif", full.names = T)
dtm_rasts <- lapply(dtm_tiles, rast)
dtm_sprc <- sprc(dtm_rasts)
dtm <- mosaic(dtm_sprc)
plot(dtm, main = "Digital Terrain Model", las = 1)
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
writeRaster(dtm, file.path(lidar_dir, "out_files", "dtm.tif"), overwrite=TRUE)
```

Next we’ll repeat that process for the DSMs

``` r
dsm_dir <- file.path(lidar_dir, "dsms")
dsm_tiles <- list.files(dsm_dir, pattern = "*.tif", full.names = T)
dsm_rasts <- lapply(dsm_tiles, rast)
dsm_sprc <- sprc(dsm_rasts)
dsm <- mosaic(dsm_sprc)
plot(dsm, main = "Digital Surface Model", las = 1)
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
writeRaster(dsm, file.path(lidar_dir, "out_files", "dsm.tif"), overwrite=TRUE)
```

Those look very similar, so to check that our DSM is inclusive of
vegetation, we can subtract the DTM from the DSM to just look at the
canopy height:

``` r
chm <- dsm-dtm
plot(chm, main = "Canopy Height", las = 1)
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## **3 The VisiMod Workflow**

------------------------------------------------------------------------

### **3.1 Prepare Digital Elevation Models**

The first step in the VisiMod workflow is to run your DTM and DSM
through the `prep_dems()` function. This function checks if the two
input SpatRasters are spatially aligned, in appropriate projected
coordinate systems (where linear units are measured in meters), and if
any interior NA values need to be filled. The function prints useful
update messages, which are suppressed here, which let users know which,
if any, requirements are not met.

In this example, the DSM does have some interior NAs (11, to be exact),
so the `prep_dems()` function will let users know that the
`"dtm and dsm are spatially aligned"`, that the
`"dtm and dsm have suitable coordinate systems"`, and that the
`"dsm has 11 NAs. Removing..."` (among other printed messages). Knowing
that the DSM has interior NAs means we must use the output DSM from
`prep_dems()` as we continue through the VisiMod workflow, because
interior NA values are problematic for subsequent functions.

``` r
dtm <- rast(file.path(lidar_dir, "out_files", "dtm.tif"))
dsm <- rast(file.path(lidar_dir, "out_files", "dsm.tif"))

pdems <- prep_dems(dtm, dsm, file.path(lidar_dir, "out_files","dtm_filled.tif"), file.path(lidar_dir, "out_files","dsm_filled.tif"))
```

### **3.2 Generate Points**

The second step in the VisiMod workflow is to generate x, y locations to
be used as training and validation points for building one or more VI
models. The user-defined number of points are randomly distributed
within the study area at an appropriate distance from the study area
boundary. This distance is controlled by the `max_vi_rad` parameter
which indicates the maximum distance to which VI will be calculated (and
therefore modeled). As `max_vi_rad` increases, processing time will
increase exponentially.

In this example we will generate 200 points with a `max_vi_rad` of 250,
which will allow us to look at smaller radii as well. We’ll convert the
resulting data.frame to a `terra` vect to plot it on top of our DSM.

**NOTE:** We are using the returned DSM (`pdems$dsm`) from `prep_dems()`
because it had interior NAs filled. We could use `pdems$dtm` for our DTM
as well, but since there were no NAs filled `pdems$dtm == dtm`

``` r
gp <- gen_pts(dtm, pdems$dsm, 200, 250)
head(gp)
```

<div class="kable-table">

|        x |       y |
|---------:|--------:|
| 255115.4 | 4210464 |
| 254364.2 | 4210458 |
| 255523.9 | 4211325 |
| 254305.3 | 4211711 |
| 254941.6 | 4211289 |
| 255155.5 | 4210404 |

</div>

While the output of `gen_pts()` is a dataframe, we can vectorize it
using `v <- vect(cbind(gp$x, gp$y))` to plot the output, overlaid on our
DSM:

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### **3.3 Calculate VI**

The next step in the VisiMod workflow is `calc_vi()` which calculates VI
at each point location. VI represents the proportion of area visible
from a given point relative to the total area within a viewing distance
of interest, ranging from 0 (no visibility) to 1 (complete visibility).
This function can be used to calculate omnidirectional VI (in 360
degrees surrounding each point) or directional VI (within a wedge of
specified field of view, and in a specified view-direction (azimuth)).

For this example we will look at omnidirectional VI and set `vi_type`
equal to `omnidir`. We will look at multiple distances, not to exceed
the original `max_vi_rad` set in `gen_pts()`. We will also set
`cores = 4` to run this using just 4 cores.

``` r
cv <- calc_vi(dtm, pdems$dsm, gp, "omnidir", c(125, 250), cores = 4)
```

While the output of `calc_vi()` is a dataframe, we can vectorize it
using `v <- vect(cbind(cv$x, cv$y))` to plot the output, and color the
points according to VI:

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

### **3.4 Generate Model Predictors**

Following the calculation of VI at each point, next in the VisiMod
workflow is the generation of predictors at each point. The
`gen_preds()` function generates and extracts a suite of visibility
predictor pixel values at each input point location. The vi_type should
match that which was used for `calc_vi()`. In this case we will continue
with `"omnidir"` and use our output dataframe from `calc_vi()` as the
input points. Since our input DTM and DSM are at a 1 m spatial
resolution we set our aggregation factor (`agg_fact`) equal to 10L. This
means the predictor rasters will be generated at a spatial resolution
10x the input (in this case, 10 m). We recommend using this combination
of spatial resolutions, but other combinations can be explored by
changing the `agg_fact` parameter or the resolution of data used.

More details on `gen_preds()` can be found in the function
documentation. For more details see Table 1 in [Mistick et.
al. 2023](https://content.csbs.utah.edu/~pdennison/reprints/denn/2023_Mistick_etal_IJGIS.pdf).
*Only local and focal variables from Mistick et. al. 2023 are included
in VisiMod.*

**NOTE:** The maximum radius we are looking to in this example is 250 m.
Our full study area is 2000 m x 2000 m therefore the area in which our
points were randomly sampled is 1750 m x 1750 m. Some points may be
exactly 250 m from the study area boundary, which is *less than* the
largest radius to which focal predictors are generated (which is 32
pixels or 320 m, at 10 m spatial resolution). `gen_preds()` will return
NA values for points that are \< 32 pixels from the study area edge.

``` r
preds <- gen_preds(dtm, pdems$dsm, cv, "omnidir", agg_fact = 10L, save_dir = file.path(lidar_dir, "out_files"))

# show the first 5 rows and 10 cols of the output data.frame
head(preds$pred_pts)
```

<div class="kable-table">

|        x |       y | elevation |     slope | slope_derivative |  curvature | curvature_plan | curvature_prof | aspect_sin | aspect_cos | slope_aspect_sin | slope_aspect_cos |       ch |        cc | elevation_mean_2 | slope_mean_2 | slope_derivative_mean_2 | curvature_mean_2 | curvature_plan_mean_2 | curvature_prof_mean_2 | aspect_sin_mean_2 | aspect_cos_mean_2 | slope_aspect_sin_mean_2 | slope_aspect_cos_mean_2 | ch_mean_2 | cc_mean_2 | elevation_sd_2 | slope_sd_2 | slope_derivative_sd_2 | curvature_sd_2 | curvature_plan_sd_2 | curvature_prof_sd_2 | aspect_sin_sd_2 | aspect_cos_sd_2 | slope_aspect_sin_sd_2 | slope_aspect_cos_sd_2 |   ch_sd_2 |   cc_sd_2 |      tpi_2 | elevation_mean_4 | slope_mean_4 | slope_derivative_mean_4 | curvature_mean_4 | curvature_plan_mean_4 | curvature_prof_mean_4 | aspect_sin_mean_4 | aspect_cos_mean_4 | slope_aspect_sin_mean_4 | slope_aspect_cos_mean_4 | ch_mean_4 | cc_mean_4 | elevation_sd_4 | slope_sd_4 | slope_derivative_sd_4 | curvature_sd_4 | curvature_plan_sd_4 | curvature_prof_sd_4 | aspect_sin_sd_4 | aspect_cos_sd_4 | slope_aspect_sin_sd_4 | slope_aspect_cos_sd_4 |   ch_sd_4 |   cc_sd_4 |      tpi_4 | elevation_mean_6 | slope_mean_6 | slope_derivative_mean_6 | curvature_mean_6 | curvature_plan_mean_6 | curvature_prof_mean_6 | aspect_sin_mean_6 | aspect_cos_mean_6 | slope_aspect_sin_mean_6 | slope_aspect_cos_mean_6 | ch_mean_6 | cc_mean_6 | elevation_sd_6 | slope_sd_6 | slope_derivative_sd_6 | curvature_sd_6 | curvature_plan_sd_6 | curvature_prof_sd_6 | aspect_sin_sd_6 | aspect_cos_sd_6 | slope_aspect_sin_sd_6 | slope_aspect_cos_sd_6 |   ch_sd_6 |   cc_sd_6 |      tpi_6 | elevation_mean_8 | slope_mean_8 | slope_derivative_mean_8 | curvature_mean_8 | curvature_plan_mean_8 | curvature_prof_mean_8 | aspect_sin_mean_8 | aspect_cos_mean_8 | slope_aspect_sin_mean_8 | slope_aspect_cos_mean_8 | ch_mean_8 | cc_mean_8 | elevation_sd_8 | slope_sd_8 | slope_derivative_sd_8 | curvature_sd_8 | curvature_plan_sd_8 | curvature_prof_sd_8 | aspect_sin_sd_8 | aspect_cos_sd_8 | slope_aspect_sin_sd_8 | slope_aspect_cos_sd_8 |   ch_sd_8 |   cc_sd_8 |      tpi_8 | elevation_mean_16 | slope_mean_16 | slope_derivative_mean_16 | curvature_mean_16 | curvature_plan_mean_16 | curvature_prof_mean_16 | aspect_sin_mean_16 | aspect_cos_mean_16 | slope_aspect_sin_mean_16 | slope_aspect_cos_mean_16 | ch_mean_16 | cc_mean_16 | elevation_sd_16 | slope_sd_16 | slope_derivative_sd_16 | curvature_sd_16 | curvature_plan_sd_16 | curvature_prof_sd_16 | aspect_sin_sd_16 | aspect_cos_sd_16 | slope_aspect_sin_sd_16 | slope_aspect_cos_sd_16 |  ch_sd_16 |  cc_sd_16 |    tpi_16 | elevation_mean_32 | slope_mean_32 | slope_derivative_mean_32 | curvature_mean_32 | curvature_plan_mean_32 | curvature_prof_mean_32 | aspect_sin_mean_32 | aspect_cos_mean_32 | slope_aspect_sin_mean_32 | slope_aspect_cos_mean_32 | ch_mean_32 | cc_mean_32 | elevation_sd_32 | slope_sd_32 | slope_derivative_sd_32 | curvature_sd_32 | curvature_plan_sd_32 | curvature_prof_sd_32 | aspect_sin_sd_32 | aspect_cos_sd_32 | slope_aspect_sin_sd_32 | slope_aspect_cos_sd_32 |  ch_sd_32 |  cc_sd_32 |    tpi_32 |    vi_125 |    vi_250 |
|---------:|--------:|----------:|----------:|-----------------:|-----------:|---------------:|---------------:|-----------:|-----------:|-----------------:|-----------------:|---------:|----------:|-----------------:|-------------:|------------------------:|-----------------:|----------------------:|----------------------:|------------------:|------------------:|------------------------:|------------------------:|----------:|----------:|---------------:|-----------:|----------------------:|---------------:|--------------------:|--------------------:|----------------:|----------------:|----------------------:|----------------------:|----------:|----------:|-----------:|-----------------:|-------------:|------------------------:|-----------------:|----------------------:|----------------------:|------------------:|------------------:|------------------------:|------------------------:|----------:|----------:|---------------:|-----------:|----------------------:|---------------:|--------------------:|--------------------:|----------------:|----------------:|----------------------:|----------------------:|----------:|----------:|-----------:|-----------------:|-------------:|------------------------:|-----------------:|----------------------:|----------------------:|------------------:|------------------:|------------------------:|------------------------:|----------:|----------:|---------------:|-----------:|----------------------:|---------------:|--------------------:|--------------------:|----------------:|----------------:|----------------------:|----------------------:|----------:|----------:|-----------:|-----------------:|-------------:|------------------------:|-----------------:|----------------------:|----------------------:|------------------:|------------------:|------------------------:|------------------------:|----------:|----------:|---------------:|-----------:|----------------------:|---------------:|--------------------:|--------------------:|----------------:|----------------:|----------------------:|----------------------:|----------:|----------:|-----------:|------------------:|--------------:|-------------------------:|------------------:|-----------------------:|-----------------------:|-------------------:|-------------------:|-------------------------:|-------------------------:|-----------:|-----------:|----------------:|------------:|-----------------------:|----------------:|---------------------:|---------------------:|-----------------:|-----------------:|-----------------------:|-----------------------:|----------:|----------:|----------:|------------------:|--------------:|-------------------------:|------------------:|-----------------------:|-----------------------:|-------------------:|-------------------:|-------------------------:|-------------------------:|-----------:|-----------:|----------------:|------------:|-----------------------:|----------------:|---------------------:|---------------------:|-----------------:|-----------------:|-----------------------:|-----------------------:|----------:|----------:|----------:|----------:|----------:|
| 254263.8 | 4210661 |  1951.833 | 10.359193 |         7.324886 |  0.0031312 |      0.0019268 |      0.0012044 |  0.8731326 | -0.4821221 |         9.060610 |        -4.969800 | 1.198852 | 0.2923046 |         1951.416 |    10.220500 |                4.429552 |        0.0016456 |             0.0016961 |            -0.0000505 |         0.8347233 |        -0.5062324 |                8.644056 |               -5.001758 |  1.336159 | 0.3097341 |       2.362466 |  1.2073568 |              2.266899 |      0.0010968 |           0.0007257 |           0.0008902 |       0.1226684 |       0.1718581 |             2.0940004 |             1.2622672 | 0.3422813 | 0.0985891 |  0.5074891 |         1950.785 |     9.959498 |                4.218534 |        0.0014986 |             0.0017671 |            -0.0002685 |         0.7593082 |        -0.5244276 |                7.716042 |              -5.0106380 | 1.3477078 | 0.3138467 |       4.094175 |   1.264998 |              2.014707 |      0.0013582 |           0.0009914 |           0.0008639 |       0.2680173 |       0.2714172 |              3.226645 |             2.3379604 | 0.4376403 | 0.1301319 |  1.3241342 |         1949.947 |     9.743723 |                4.844677 |        0.0015970 |             0.0016916 |            -0.0000946 |         0.6701967 |        -0.4892141 |                6.590922 |              -4.5876464 | 1.3808618 | 0.3270973 |       5.426010 |   1.535112 |              2.489922 |      0.0016902 |           0.0012264 |           0.0009824 |       0.4087229 |       0.3766030 |              4.426467 |              3.598010 | 0.4864766 | 0.1406098 |  2.3455662 |         1948.859 |     9.433940 |                5.712109 |        0.0009347 |             0.0010984 |            -0.0001637 |         0.6016622 |        -0.4438090 |                5.635824 |              -4.0129092 | 1.5287503 | 0.3628469 |       6.415300 |   1.914557 |              3.230090 |      0.0024453 |           0.0019508 |           0.0013119 |       0.4749188 |       0.4640459 |              4.919154 |              4.538318 | 0.6328553 | 0.1657834 |  3.7263379 |          1947.039 |      8.301325 |                 6.649167 |        -0.0003172 |              0.0000799 |             -0.0003971 |          0.5803333 |         -0.4388264 |                 4.870579 |                -3.018843 |  1.6193043 |  0.3823686 |        9.407798 |    2.921781 |               5.126100 |       0.0038717 |            0.0032485 |            0.0019527 |        0.4290552 |        0.5356469 |               4.322650 |               5.094846 | 0.8889679 | 0.2124571 |  5.445934 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |       NaN |       NaN |       NaN | 0.0060097 | 0.0040947 |
| 254288.3 | 4210352 |  1928.460 |  6.370967 |         3.686711 |  0.0004039 |      0.0007203 |     -0.0003164 |  0.9847680 | -0.1573422 |         6.272604 |        -1.008496 | 1.692592 | 0.4236675 |         1928.422 |     6.276399 |                3.568702 |       -0.0001252 |             0.0004093 |            -0.0005345 |         0.9521462 |        -0.2480639 |                5.975206 |               -1.542394 |  1.097103 | 0.2599005 |       1.464444 |  0.8548997 |              1.161014 |      0.0014041 |           0.0014396 |           0.0005817 |       0.0572972 |       0.1662728 |             0.8715543 |             1.1031210 | 0.9393680 | 0.2236805 |  0.0444034 |         1928.470 |     6.279349 |                3.607016 |       -0.0001772 |             0.0003338 |            -0.0005111 |         0.8479140 |        -0.4013857 |                5.296367 |              -2.5226420 | 0.8900159 | 0.2083157 |       2.510198 |   1.283039 |              1.706525 |      0.0017593 |           0.0019474 |           0.0007576 |       0.1807543 |       0.2976103 |              1.483960 |             2.1250106 | 0.8526127 | 0.2063038 | -0.0307316 |         1928.542 |     6.563533 |                3.705579 |        0.0001291 |             0.0004936 |            -0.0003645 |         0.7787038 |        -0.4848907 |                5.055928 |              -3.2206460 | 0.8439041 | 0.1983430 |       3.530409 |   1.486693 |              1.767980 |      0.0022244 |           0.0021303 |           0.0008775 |       0.2303704 |       0.3259391 |              1.752454 |              2.511955 | 0.8396843 | 0.2047369 | -0.1246971 |         1928.544 |     6.673140 |                3.861244 |       -0.0001156 |             0.0001962 |            -0.0003118 |         0.7267362 |        -0.5294928 |                4.843651 |              -3.5421359 | 0.8880172 | 0.2064167 |       4.597423 |   1.611867 |              2.054921 |      0.0028783 |           0.0025255 |           0.0010960 |       0.2733109 |       0.3414683 |              2.046490 |              2.632245 | 0.8934676 | 0.2140930 | -0.1142173 |          1929.334 |      6.470890 |                 3.726801 |        -0.0001698 |              0.0001083 |             -0.0002781 |          0.7556785 |         -0.4943258 |                 5.015428 |                -2.989332 |  0.8940628 |  0.2072569 |        8.844855 |    2.551044 |               2.870407 |       0.0028062 |            0.0024304 |            0.0011158 |        0.2666490 |        0.3371165 |               2.817011 |               2.522330 | 0.8948724 | 0.2183528 | -1.138112 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |       NaN |       NaN |       NaN | 0.0011612 | 0.0002903 |
| 254291.3 | 4211399 |  1999.444 | 15.784904 |         3.601424 | -0.0003703 |     -0.0002037 |     -0.0001666 | -0.9780587 | -0.2079809 |       -15.439310 |        -3.279440 | 1.107763 | 0.2237061 |         1999.455 |    14.402943 |                9.727591 |        0.0012762 |             0.0005032 |             0.0007729 |        -0.9578941 |        -0.2563556 |              -13.940149 |               -3.413714 |  1.166653 | 0.2570632 |       3.577170 |  2.4808947 |              9.416227 |      0.0036184 |           0.0017809 |           0.0024714 |       0.0489640 |       0.0964109 |             2.6501729 |             0.4381207 | 0.3993606 | 0.1335385 | -0.0090739 |         1998.852 |    13.075500 |               11.321369 |        0.0021940 |             0.0013588 |             0.0008352 |        -0.6774701 |        -0.3376798 |               -9.834587 |              -3.9141722 | 1.3450151 | 0.3077934 |       5.483619 |   3.287305 |              8.254641 |      0.0056810 |           0.0036090 |           0.0031386 |       0.6150628 |       0.1838810 |              8.138157 |             0.9315933 | 0.4818301 | 0.1477967 |  0.8550405 |         1997.776 |    12.458296 |               11.219459 |        0.0012406 |             0.0007122 |             0.0005284 |        -0.4628885 |        -0.3882535 |               -6.509739 |              -4.2856271 | 1.5333190 | 0.3628450 |       6.493644 |   3.414697 |              7.732490 |      0.0056726 |           0.0035297 |           0.0030677 |       0.7731859 |       0.1939598 |             10.257471 |              0.963379 | 0.5542117 | 0.1630398 |  2.2313450 |         1996.786 |    11.734696 |               10.620705 |        0.0001554 |            -0.0000657 |             0.0002211 |        -0.3013558 |        -0.4705843 |               -4.498329 |              -4.7528074 | 1.6712435 | 0.3986053 |       7.081897 |   3.819942 |              7.621822 |      0.0062386 |           0.0041442 |           0.0031445 |       0.7966090 |       0.2368535 |             10.417678 |              1.165531 | 0.5975976 | 0.1630718 |  3.4657175 |          1996.212 |     10.856684 |                 9.822775 |        -0.0000071 |             -0.0000032 |             -0.0000039 |         -0.2669345 |         -0.5646897 |                -4.084263 |                -5.516887 |  1.8461593 |  0.4334607 |       11.179761 |    3.924785 |               6.973664 |       0.0062906 |            0.0041871 |            0.0031896 |        0.7358785 |        0.2626167 |               9.026758 |               2.179994 | 0.7939511 | 0.1900046 |  3.456285 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |       NaN |       NaN |       NaN | 0.0299873 | 0.0362058 |
| 254305.3 | 4211711 |  2017.532 | 15.289977 |         8.165973 |  0.0022115 |      0.0003536 |      0.0018579 |  0.9887375 |  0.1459603 |        15.114559 |         2.254831 | 1.439401 | 0.3286644 |         2017.261 |    13.960216 |                8.540890 |        0.0007697 |             0.0005306 |             0.0002391 |         0.9853590 |         0.0786335 |               13.760377 |                1.333835 |  1.495542 | 0.3591694 |       3.388481 |  1.9442130 |              3.146895 |      0.0021294 |           0.0010886 |           0.0014967 |       0.0112876 |       0.1490709 |             1.9345604 |             1.9383094 | 0.3837156 | 0.1462836 |  0.3245516 |         2016.948 |    13.008895 |                7.585759 |        0.0006807 |             0.0007066 |            -0.0000258 |         0.9489996 |        -0.0131656 |               12.479671 |               0.3396027 | 1.5162086 | 0.3632880 |       5.608843 |   2.678873 |              3.873775 |      0.0022619 |           0.0019805 |           0.0017006 |       0.0881799 |       0.2992644 |              3.169266 |             3.2166845 | 0.4450099 | 0.1520276 |  0.7209452 |         2016.540 |    12.207375 |                7.840884 |        0.0004377 |             0.0002727 |             0.0001649 |         0.8290569 |        -0.0611060 |               10.645984 |              -0.1333946 | 1.5951975 | 0.3868907 |       7.387960 |   3.377986 |              4.327071 |      0.0034662 |           0.0033318 |           0.0018414 |       0.2856843 |       0.4744198 |              5.086698 |              4.582389 | 0.5317587 | 0.1640366 |  1.2119924 |         2016.171 |    11.914661 |                7.920002 |        0.0000647 |            -0.0001936 |             0.0002583 |         0.5997756 |        -0.0627569 |                8.005327 |              -0.2412763 | 1.6278178 | 0.3976113 |       8.483287 |   3.580241 |              4.327323 |      0.0038655 |           0.0035953 |           0.0019041 |       0.5533372 |       0.5744436 |              7.556163 |              5.786571 | 0.5735300 | 0.1659095 |  1.6656992 |          2016.075 |     12.370166 |                 7.715543 |         0.0002648 |              0.0001544 |              0.0001104 |          0.0741323 |         -0.1232886 |                 1.399256 |                -1.113328 |  1.5951932 |  0.3830913 |       11.472190 |    3.418378 |               5.014583 |       0.0041723 |            0.0035764 |            0.0020001 |        0.7254079 |        0.6736731 |               9.756757 |               8.150265 | 0.5793813 | 0.1621555 |  1.499574 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |       NaN |       NaN |       NaN | 0.0022613 | 0.0125083 |
| 254325.0 | 4210904 |  1958.530 | 11.595578 |        13.142299 |  0.0005380 |     -0.0027065 |      0.0032445 |  0.8763575 |  0.4812325 |        10.156502 |         5.590870 | 1.928525 | 0.5123694 |         1958.245 |    10.538185 |               11.777655 |        0.0020184 |             0.0007371 |             0.0012813 |         0.9098608 |         0.3644158 |                9.460029 |                4.276494 |  2.269774 | 0.5644329 |       2.600482 |  2.6037301 |              4.652311 |      0.0027849 |           0.0015198 |           0.0026753 |       0.0551018 |       0.1890925 |             1.8986609 |             2.5395444 | 0.6116671 | 0.1478378 |  0.3651330 |         1957.633 |     8.415659 |               11.865016 |        0.0007020 |             0.0005634 |             0.0001386 |         0.8755558 |         0.0949689 |                7.427119 |               2.0696829 | 2.4465145 | 0.5701533 |       3.880380 |   3.458194 |              6.550676 |      0.0042429 |           0.0021953 |           0.0036662 |       0.1407855 |       0.4550455 |              3.063932 |             3.7535365 | 0.8768229 | 0.1592610 |  1.1643981 |         1957.188 |     7.226276 |               10.517052 |        0.0003797 |             0.0001415 |             0.0002382 |         0.6851043 |        -0.1570844 |                5.410114 |               0.2278188 | 2.3768329 | 0.5416952 |       4.351192 |   3.386609 |              6.489480 |      0.0044978 |           0.0025590 |           0.0034126 |       0.4214690 |       0.5761680 |              4.134831 |              4.175509 | 0.9054584 | 0.1711174 |  1.6270246 |         1956.807 |     6.955333 |                9.252425 |        0.0003313 |             0.0000636 |             0.0002677 |         0.5373638 |        -0.3109314 |                4.031821 |              -1.1909850 | 2.3680869 | 0.5445369 |       4.479870 |   3.040157 |              5.957042 |      0.0042272 |           0.0024768 |           0.0031412 |       0.5057544 |       0.6012456 |              4.472426 |              4.480537 | 0.9004181 | 0.1746309 |  2.0492871 |          1955.424 |      6.683509 |                 6.287687 |        -0.0004036 |             -0.0002341 |             -0.0001695 |          0.4334018 |         -0.5420433 |                 2.906334 |                -3.286553 |  1.9237188 |  0.4461700 |        6.682045 |    2.584912 |               4.715198 |       0.0034137 |            0.0026223 |            0.0020545 |        0.4851536 |        0.5324835 |               3.907253 |               4.106309 | 1.1054034 | 0.2450034 |  3.582296 |          1960.275 |           NaN |                      NaN |        -0.0000065 |              0.0002245 |             -0.0002310 |                NaN |                NaN |                      NaN |                      NaN |  1.5133413 |  0.3571470 |        15.69190 |         NaN |                    NaN |       0.0032893 |            0.0027241 |            0.0016490 |              NaN |              NaN |                    NaN |                    NaN | 0.9777541 | 0.2317009 | -3.377805 | 0.0005093 | 0.0009473 |
| 254333.4 | 4210496 |  1931.688 |  4.359625 |         3.156425 |  0.0013720 |      0.0017351 |     -0.0003631 |  0.6837411 | -0.7232397 |         2.986383 |        -3.148468 | 1.673918 | 0.3437439 |         1931.535 |     4.737705 |                3.536437 |        0.0002240 |             0.0008019 |            -0.0005779 |         0.6555355 |        -0.7044130 |                3.170801 |               -3.275211 |  1.671049 | 0.3907874 |       1.056822 |  0.7546908 |              1.275338 |      0.0016784 |           0.0014858 |           0.0006172 |       0.1961437 |       0.1955229 |             1.2506239 |             0.8467805 | 0.5738219 | 0.1351505 |  0.1833067 |         1931.546 |     5.358494 |                5.109109 |       -0.0005833 |             0.0000799 |            -0.0006632 |         0.7186743 |        -0.5505743 |                3.990154 |              -2.7336485 | 1.5938637 | 0.3717571 |       1.988212 |   1.565894 |              3.540596 |      0.0029778 |           0.0027395 |           0.0012707 |       0.2793159 |       0.3226591 |              2.092523 |             1.8509498 | 0.8679999 | 0.1994957 |  0.1373521 |         1931.802 |     5.877242 |                5.171654 |       -0.0003329 |             0.0001633 |            -0.0004962 |         0.7410074 |        -0.4622340 |                4.536247 |              -2.3555826 | 1.4760641 | 0.3400438 |       2.983695 |   2.026385 |              4.138360 |      0.0032659 |           0.0027884 |           0.0014056 |       0.3379701 |       0.3528857 |              2.746296 |              2.242713 | 0.9227017 | 0.2134133 | -0.2309053 |         1932.071 |     6.055916 |                4.996059 |       -0.0000658 |             0.0003556 |            -0.0004214 |         0.7380349 |        -0.4466984 |                4.631198 |              -2.2773411 | 1.3424925 | 0.3098183 |       3.945495 |   2.163860 |              4.103073 |      0.0033849 |           0.0026995 |           0.0015622 |       0.3542082 |       0.3624177 |              2.975735 |              2.430288 | 0.9345290 | 0.2195344 | -0.5919962 |          1932.766 |      6.190685 |                 4.165531 |        -0.0002426 |              0.0001006 |             -0.0003433 |          0.6753452 |         -0.5321210 |                 4.276819 |                -2.815579 |  1.0676110 |  0.2463010 |        7.529263 |    2.505027 |               4.147017 |       0.0031597 |            0.0027366 |            0.0014471 |        0.3360303 |        0.3847634 |               3.059549 |               3.005112 | 0.9791103 | 0.2338113 | -1.320014 |          1935.928 |      6.515153 |                      NaN |        -0.0001218 |              0.0000835 |             -0.0002053 |          0.6097426 |         -0.5560847 |                 4.291691 |                -2.915192 |  0.8602406 |  0.1982415 |        15.45976 |    3.248733 |                    NaN |       0.0029028 |            0.0025174 |            0.0012874 |        0.3567888 |         0.437928 |               3.689489 |               3.532265 | 0.9322108 | 0.2284307 | -5.310427 | 0.0070690 | 0.0017673 |

</div>

``` r
# plot an example predictor raster
plot(preds$pred_rast[["slope_aspect_sin_mean_8"]], main = "Predictor Raster: Focal Mean Slope Aspect Sine r = 8 pixels ", las = 1)
```

![](README_files/figure-gfm/gen-preds-1.png)<!-- -->

### **3.5 Build Random Forest Model**

Now that we have points with VI and predictor values we can move to the
next step in VisiMod which is to build a predictive model using
`mod_vi()`. This function will build a random forest model to predict VI
based on the predictor values at each point. We will use the output
data.frame from `gen_preds()` as input into `mod_vi()` because it
contains all our points, inclusive of VI and predictor values.

We will set `cross_validate = TRUE` and `tune = FALSE`. These parameters
can be toggled depending on the user’s goal. In this case we’d like to
perform a 4-fold spatial cross-validation procedure to assess model
performance, but we don’t want the hyperparameter tuning to slow down
our processing significantly.

Since we have calculated VI at two different radii we will call
`mod_vi()` twice, once for each radii.

``` r
mod125 <- mod_vi(preds$pred_pts, 125, cross_validate = TRUE, tune = FALSE, num_cores = 4)
mod250<- mod_vi(preds$pred_pts, 250, cross_validate = TRUE, tune = FALSE, num_cores = 4)
```

Then we can look at model performance:

<div class="kable-table">

|        r2 |      rmse |     nrmse | VI Radius (m) |
|----------:|----------:|----------:|--------------:|
| 0.8079203 | 0.0700826 | 0.1481233 |           125 |
| 0.6833429 | 0.0558025 | 0.1512293 |           250 |

</div>

### **3.6 Map VI**

The final step in the VisiMod workflow is to map VI across the entire
study area, using the random forest model built in `mod_vi()` and the
predictor raster stack generated with `gen_preds()`.

We will run this function twice, once per radius, to map VI at both our
radii (125 m and 250 m).

``` r
map125 <- map_vi(mod125$ranger_mod, preds$pred_rast, 4, FALSE)
plot(map125, main = "Omnidirectional VI at a 125 m radius", las=1)
```

![](README_files/figure-gfm/map-vi-1.png)<!-- -->

``` r
map250 <- map_vi(mod250$ranger_mod, preds$pred_rast, 4, FALSE)
plot(map250, main = "Omnidirectional VI at a 250 m radius", las=1)
```

![](README_files/figure-gfm/map-vi-2.png)<!-- -->

## **4 Using VisiMod() to generate an RGB VI map**

------------------------------------------------------------------------

### **4.1 Run the Function**

The all-encompassing `VisiMod()` function runs through the entire
VisiMod workflow and returns a list of rasters with mapped VI values.
All users must do is provide an input DTM and DSM, then they can run
`VisiMod()`. In this example we are going to generate three VI maps, at
azimuths of 0 degrees (looking north), 120 degrees (looking southeast),
and 240 degrees (looking southwest) (`vi_azi = c(0,120,240)`) each with
a field of view of 120 degrees (`vi_fov = 120`).

We recommend running `VisiMod()` if you are only interested in the
output VI map(s). If you are interested in saving intermediary outputs
(the points used in modeling, the ranger model objects, etc.) we
recommend running the workflow described in section 3.

`VisiMod()` makes assumptions for input parameters in two of the
internal functions: `gen_preds` and `mod_vi()`. For `gen_preds()` a
default aggregation factor of 10 is assumed. For `mod_vi()`
`cross_validate` and `tune` are assumed to be `FALSE`. If users wish to
toggle these parameters they should run the workflow described in
section 3.

Running this example with 4 cores resulted in a total run time of 3
minutes and 25 seconds. Many messages will print indicating progress.
When working with larger datasets it is not unusual for `calc_vi()` and
`gen_preds()` to perform slowly. Also notable, the step that generates
predictors at 32 pixel radii tends to be slow..

``` r
# define your input params
dtm <- rast(file.path(lidar_dir, "out_files","dtm.tif"))
dsm <- rast(file.path(lidar_dir, "out_files","dsm.tif"))

visi_out <- VisiMod(dtm=dtm, dsm=dsm,
                    num_pts = 200,dist =125,
                    vi_type = "directional_single",
                    vi_fov = 120, vi_azi = c(0,120,240),
                    save_dir =paste0(lidar_dir, "/out_files"), cores = 4)
```

To visualize the area in which VI and predictors are being generated, we
can use the wedge() function to draw wedges. In the example below three
wedges will be drawn, one for each azimuth in our previous call to
VisiMod() (0, 120, 240). We will color the wedges according to the RGB
channel with which they will be visualized in the next step.

``` r
# plot the dsm as our basemap (zoom in a bit)
plot(dsm, ext = ext(dsm)-750, main = "Example of Wedge() Function")

# use the midpoint as the point from which our wedge will be drawn 
x <- (ext(dtm)[2]+ext(dtm)[1])/2
y <- (ext(dtm)[4]+ext(dtm)[3])/2

w <- wedge(x = x[[1]], y = y[[1]], 125, 0, 120)
plot(w, col = "red", add=TRUE)

w <- wedge(x = x[[1]], y = y[[1]], 125, 120, 120)
plot(w, col = "green", add=TRUE)

w <- wedge(x = x[[1]], y = y[[1]], 125, 240, 120)
plot(w, col = "blue", add=TRUE)
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### **4.2 Using plotRGB() to plot VI**

The first option for plotting VI is using the plotRGB() function from
the `terra` library. This is an easy way to plot 3 band rasters

Since `VisiMod()` returns a list of three single-band rasters, we must
combine them into a multiband raster before we can plot an RGB map. We
apply a scale of 1 to tell the `plotRGB()` function that our values
range from 0-1 (the default assumes 0-255). We have also applied a
stretch to the plot, since VI values tend to not cover the full 0-1
range (100% visibiliy, in this study area, is rare), we stretch the
values so the pattern of colors is more easily discerned. The 32 pixel
buffer is also left in place so users can see that predictions are
returned with buffer equal to the largest focal predictor radius
considered.

``` r
vi <- c(visi_out[[1]],visi_out[[2]],visi_out[[3]])
plotRGB(vi, r=1,g=2,b=3,
        scale = 1, stretch="lin",
        axes = TRUE, mar = c(3,3), 
        xlab = "UTM Easting", ylab = "UTM Northing",
        main = "Directional VI")
```

![](README_files/figure-gfm/plot-visimod-output-1.png)<!-- -->

### **4.3 Using ggplot2 to plot RGB VI map**

`ggplot2` is a powerful and popular library for creating graphics. While
not built explicitly to plot raster data, it can easily plot rasters
that have been converted to a data.frame. For smaller datasets this
option is beneficial in that it allows for more control over the plot,
and users may already be more familiar with this plotting library. The
example below excludes the original border, unlike the `plotRGB()`
example. This example also does not apply a stretch to the values,
meaning the color ramp is not stretch from 0 to 1 but reflects 0 to the
maximum VI for that color channel.

``` r
vi <- c(visi_out[[1]],visi_out[[2]],visi_out[[3]])
vi_df <- as.data.frame(vi, xy=TRUE, na.rm) 
vi_df <- na.omit(vi_df)
p <- ggplot(data=vi_df, aes(x=x,y=y)) + 
  geom_raster(fill = rgb(r = vi_df$vi_d125_f120_a0,
                         g=vi_df$vi_d125_f120_a120,
                         b=vi_df$vi_d125_f120_a240, maxColorValue =1)) +
  ylab("UTM Northing")+
  xlab("UTM Easting") + 
  ggtitle("Directional VI") +
  scale_x_continuous(expand = expansion(c(0,0)))+
  scale_y_continuous(expand = expansion(c(0,0)))+
  theme(panel.background = element_rect(color = NA, fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        text = element_text(size = 14),
        aspect.ratio = 1) +
  guides(fill = guide_legend(override.aes = list(size = 5)))+
  labs(fill = "Radius (m)", size = "Field of View \n(degrees)") 
p
```

![](README_files/figure-gfm/ggplot-rgb-1.png)<!-- -->
