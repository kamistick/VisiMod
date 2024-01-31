VisiMod: Map Modeled Visibility Index Across Wildland Landscapes
================
Katherine Mistick
2024-01-31

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

The [VisiMod package](https://github.com/kamistick/VisiMod/) was built
to facilitate visibility analyses. It provides a suite of tools that
allows users to build predictive models that can map visibility (in
terms of visibility index (VI) where VI = (visible area)/(total area))
across entire landscapes. [Mistick et.
al. 2023](https://content.csbs.utah.edu/~pdennison/reprints/denn/2023_Mistick_etal_IJGIS.pdf)
is the original research that led to the development of VisiMod.

Instead of calculating VI by running viewsheds from every pixel in a
rasterized study area, VisiMod allows users to build predictive models
based on a sample of points (at which viewsheds are run using
[terra::viewshed()](https://rdrr.io/cran/terra/man/viewshed.html) to
predict VI at every pixel.

VisiMod is designed to be used in wildland (i.e. undeveloped) areas. It
relies on a digital surface model (DSM) that represents elevations of
the ground and any aboveground features present on the landscape,
inclusive of vegetation. If VisiMod were deployed for developed areas
any aboveground obstructions would be treated as vegetation. We have not
tested the performance of VisiMod in developed areas.

The VisiMod package can be installed by calling
`devtools::install_github("kamistick/VisiMod")`in the RStudio console.
To run the following examples the following packages must be loaded:
`VisiMod`, `terra`, `httr2`, `lidR`, `future`, `ggplot2`

## **2 Download lidar data to prepare for VisiMod workflow**

------------------------------------------------------------------------

VisiMod only requires users to provide a digital terrain model (DTM) and
digital surface model (DSM) to run the full workflow. While fine-scale
resolution DTMs are available to download download (e.g., from the USGS
3D Elevation Program or OpenTopography), fine-scale DSMs are not widely
available, and must be derived from some other data source, such as
airborne lidar. The example below walks a user through downloading and
processing four USGS 3DEP lidar tiles for use in the VisiMod workflow.

We’ll use The National Map (TNM) API, provided by the USGS, to access
3DEP lidar data. For more information on how to use the TNM API, [click
here](https://apps.nationalmap.gov/tnmaccess/#/). To interact with the
TNM GUI and download data, [click
here](https://apps.nationalmap.gov/downloader/).

For this example we’re going to look at a 4km<sup>2</sup> area in
southern Utah, USA:

![image](images/leaflet.png) \*The above figure is a static image of an
interactive leaflet map that does not render in this GitHub README.md

### **2.1 Set Up Directories**

First, we will set up our directories. In this example the main
directory `lidar_dir` is a folder with four subfolders: `clean`, `dtms`,
`dsms`, and `out_files`. We will download, prepare, and generate new
files in these directories throughout these examples.

``` r
lidar_dir <- "/path/to/your/Data"
```

``` r
for (folder in c("clean", "dtms", "dsms", "out_files")){
  dir.create(file.path(lidar_dir, folder))
}
```

### **2.2 Download lidar using The National Map API**

Next we will search for available lidar products within our bounding
box. The API requires a bounding box in the format
\[xMin,yMin,xMax,yMax\]. If you have a polygon that defines your study
area you can read that in using `my_vector <- terra:vect()` and then
call `ext(my_vector)` to get the extent in (xmin, xmax, ymin, ymax)
order. TNM API requires the bounding box be defined in decimal degrees,
if your polygon is in a different coordinate reference system you can
use `project(my_vector, "+proj=longlat +datum=WGS84")` to get the extent
in decimal degrees.

``` r
bounding_box <- "-113.796282,38.008962,-113.785245,38.018234"
product_format <- "LAS,LAZ"
url <- paste0("https://tnmaccess.nationalmap.gov/api/v1/products?bbox=", bounding_box, "&prodFormats=", product_format)
req <- request(url)
resp <- req %>% req_perform()
```

Next we’ll grab the download links for the LAZ point clouds, of which
there are 4. We’ll also need the filenames so we can save each file. We
have found that when downloading tiles using TNM API there is a 50 file
limit per request. While that limit does not affect this example, we
recommend for larger jobs to use a smaller bounding box until each API
request is \< 50 tiles

``` r
for (x in 1:length(resp_body_json(resp)[2]$items)){
  # get download link
  download_link <- resp_body_json(resp)[2]$items[[x]]$urls$LAZ
  
  #get file name without white spaces
  filename <- gsub(" ", "", resp_body_json(resp)[2]$items[[x]]$title)
  
  #use download.file() to download each laz tile
  print(paste0("Downloading ", filename))
  download.file(download_link, paste0(lidar_dir, "/", filename, ".laz"), mode = "wb")

}
```

### **2.3 Process lidar with lidR in R**

Next we’ll use the [lidR package](https://r-lidar.github.io/lidRbook/)
to process our point cloud tiles to derive a DTM and DSM. First, we will
read the files in as a LAScatalog, and plot them. It is recommended that
you index your lidar tiles before proceeding using
`lidR:::catalog_laxindex(cat)`. Indexing allows for faster processing
and is especially advantageous when dealing with large LAScatalogs. For
more on LAScatalogs check on this vignette: [LAScatalog formal
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

First we’ll clean up the lidar data to remove any noise. This step isn’t
always strictly necessary and depends heavily upon the quality of the
data source. Some datasets are fairly clean to begin with, others may
have extremely high or low noise flagged as classification 7 or 18.
These points can be removed by filtering. However, there may be certain
high or low points that were not classified as noise by the data vendor,
but should be. In this example, we will run an additional noise filter
and remove the noise flagged points to ensure all extraneous points have
been removed. If you have a CPU with multiple cores available to work
on, we suggest using them by implementing
`plan(multisession, workers = 4L)` and `set_lidr_threads(4L)`, which
would implement 4 cores for parallel processing.

``` r
# set up parallel processing, if available
plan(multisession, workers = 4L)
set_lidr_threads(4L)

# define directoy where clean laz files will be written
clean_dir <- file.path(lidar_dir, "clean")

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
rasterize_terrain(cat_clean, res= 1, algorithm = tin())

# generate dsms
opt_output_files(cat_clean) <- file.path(dsm_dir, "dsm_{XLEFT}_{YBOTTOM}")
rasterize_canopy(cat_clean, res = 1, dsmtin())
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

pdems <- prep_dems(dtm, dsm,
                   file.path(lidar_dir, "out_files","dtm_filled.tif"),
                   file.path(lidar_dir, "out_files","dsm_filled.tif"))
```

### **3.2 Generate Points**

The second step in the VisiMod workflow is to generate x, y locations to
be used as training and validation points for building one or more VI
models. The user-defined number of points are randomly distributed
within the study area at an appropriate distance from the study area
boundary. This distance is controlled by the `max_vi_rad` parameter
which indicates the maximum distance to which VI will be calculated (and
therefore modeled). As `max_vi_rad` increases, processing time later in
the workflow will increase exponentially.

In this example we will generate 200 points with a `max_vi_rad` of 250.
We’ll convert the resulting data.frame to a `terra` vect to plot it on
top of our DSM.

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
| 254467.0 | 4210811 |
| 255059.6 | 4210596 |
| 254269.1 | 4211645 |
| 254534.8 | 4211648 |
| 255416.5 | 4211592 |
| 254809.5 | 4210351 |

</div>

While the output of `gen_pts()` is a dataframe, we can vectorize it
using `v <- vect(cbind(gp$x, gp$y))` to plot the output, overlaid on our
DSM. You can see how `max_vi_rad` ensured that no points felt within 250
m of the study area’s edge:

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### **3.3 Calculate VI**

The next step in the VisiMod workflow is `calc_vi()` which calculates VI
at each point location. VI represents the proportion of area visible
from a given point relative to the total area within a viewing radius of
interest, ranging from 0 (no visibility) to 1 (complete visibility).
This function can be used to calculate omnidirectional VI (in 360
degrees surrounding each point) or directional VI (within a wedge of
specified field of view, and in a specified view-direction (azimuth)).

For this example we will look at omnidirectional VI and set `vi_type`
equal to `omnidir`. We will look at multiple distances (125 and 250 m),
not to exceed the original `max_vi_rad` set in `gen_pts()`. We will also
set `cores = 4` to run this using just 4 cores.

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
predictor pixel values at each input point location. The `vi_type`
should match that which was used for `calc_vi()`. In this case we will
continue with `"omnidir"` and use our output dataframe from `calc_vi()`
as the input points. Since our input DTM and DSM are at a 1 m spatial
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
pixels or 320 m at 10 m spatial resolution). `gen_preds()` will return
NA values for points that are \< 32 pixels from the study area edge.

``` r
preds <- gen_preds(dtm, pdems$dsm, cv, "omnidir", agg_fact = 10L, save_dir = file.path(lidar_dir, "out_files"))

# show the first 5 rows and 10 cols of the output data.frame
head(preds$pred_pts)
```

<div class="kable-table">

|        x |       y | elevation |     slope | slope_derivative |  curvature | curvature_plan | curvature_prof | aspect_sin | aspect_cos | slope_aspect_sin | slope_aspect_cos |       ch |        cc | elevation_mean_2 | slope_mean_2 | slope_derivative_mean_2 | curvature_mean_2 | curvature_plan_mean_2 | curvature_prof_mean_2 | aspect_sin_mean_2 | aspect_cos_mean_2 | slope_aspect_sin_mean_2 | slope_aspect_cos_mean_2 | ch_mean_2 | cc_mean_2 | elevation_sd_2 | slope_sd_2 | slope_derivative_sd_2 | curvature_sd_2 | curvature_plan_sd_2 | curvature_prof_sd_2 | aspect_sin_sd_2 | aspect_cos_sd_2 | slope_aspect_sin_sd_2 | slope_aspect_cos_sd_2 |   ch_sd_2 |   cc_sd_2 |      tpi_2 | elevation_mean_4 | slope_mean_4 | slope_derivative_mean_4 | curvature_mean_4 | curvature_plan_mean_4 | curvature_prof_mean_4 | aspect_sin_mean_4 | aspect_cos_mean_4 | slope_aspect_sin_mean_4 | slope_aspect_cos_mean_4 | ch_mean_4 | cc_mean_4 | elevation_sd_4 | slope_sd_4 | slope_derivative_sd_4 | curvature_sd_4 | curvature_plan_sd_4 | curvature_prof_sd_4 | aspect_sin_sd_4 | aspect_cos_sd_4 | slope_aspect_sin_sd_4 | slope_aspect_cos_sd_4 |   ch_sd_4 |   cc_sd_4 |      tpi_4 | elevation_mean_6 | slope_mean_6 | slope_derivative_mean_6 | curvature_mean_6 | curvature_plan_mean_6 | curvature_prof_mean_6 | aspect_sin_mean_6 | aspect_cos_mean_6 | slope_aspect_sin_mean_6 | slope_aspect_cos_mean_6 | ch_mean_6 | cc_mean_6 | elevation_sd_6 | slope_sd_6 | slope_derivative_sd_6 | curvature_sd_6 | curvature_plan_sd_6 | curvature_prof_sd_6 | aspect_sin_sd_6 | aspect_cos_sd_6 | slope_aspect_sin_sd_6 | slope_aspect_cos_sd_6 |   ch_sd_6 |   cc_sd_6 |      tpi_6 | elevation_mean_8 | slope_mean_8 | slope_derivative_mean_8 | curvature_mean_8 | curvature_plan_mean_8 | curvature_prof_mean_8 | aspect_sin_mean_8 | aspect_cos_mean_8 | slope_aspect_sin_mean_8 | slope_aspect_cos_mean_8 | ch_mean_8 | cc_mean_8 | elevation_sd_8 | slope_sd_8 | slope_derivative_sd_8 | curvature_sd_8 | curvature_plan_sd_8 | curvature_prof_sd_8 | aspect_sin_sd_8 | aspect_cos_sd_8 | slope_aspect_sin_sd_8 | slope_aspect_cos_sd_8 |   ch_sd_8 |   cc_sd_8 |      tpi_8 | elevation_mean_16 | slope_mean_16 | slope_derivative_mean_16 | curvature_mean_16 | curvature_plan_mean_16 | curvature_prof_mean_16 | aspect_sin_mean_16 | aspect_cos_mean_16 | slope_aspect_sin_mean_16 | slope_aspect_cos_mean_16 | ch_mean_16 | cc_mean_16 | elevation_sd_16 | slope_sd_16 | slope_derivative_sd_16 | curvature_sd_16 | curvature_plan_sd_16 | curvature_prof_sd_16 | aspect_sin_sd_16 | aspect_cos_sd_16 | slope_aspect_sin_sd_16 | slope_aspect_cos_sd_16 |  ch_sd_16 |  cc_sd_16 |     tpi_16 | elevation_mean_32 | slope_mean_32 | slope_derivative_mean_32 | curvature_mean_32 | curvature_plan_mean_32 | curvature_prof_mean_32 | aspect_sin_mean_32 | aspect_cos_mean_32 | slope_aspect_sin_mean_32 | slope_aspect_cos_mean_32 | ch_mean_32 | cc_mean_32 | elevation_sd_32 | slope_sd_32 | slope_derivative_sd_32 | curvature_sd_32 | curvature_plan_sd_32 | curvature_prof_sd_32 | aspect_sin_sd_32 | aspect_cos_sd_32 | slope_aspect_sin_sd_32 | slope_aspect_cos_sd_32 | ch_sd_32 | cc_sd_32 | tpi_32 |    vi_125 |    vi_250 |
|---------:|--------:|----------:|----------:|-----------------:|-----------:|---------------:|---------------:|-----------:|-----------:|-----------------:|-----------------:|---------:|----------:|-----------------:|-------------:|------------------------:|-----------------:|----------------------:|----------------------:|------------------:|------------------:|------------------------:|------------------------:|----------:|----------:|---------------:|-----------:|----------------------:|---------------:|--------------------:|--------------------:|----------------:|----------------:|----------------------:|----------------------:|----------:|----------:|-----------:|-----------------:|-------------:|------------------------:|-----------------:|----------------------:|----------------------:|------------------:|------------------:|------------------------:|------------------------:|----------:|----------:|---------------:|-----------:|----------------------:|---------------:|--------------------:|--------------------:|----------------:|----------------:|----------------------:|----------------------:|----------:|----------:|-----------:|-----------------:|-------------:|------------------------:|-----------------:|----------------------:|----------------------:|------------------:|------------------:|------------------------:|------------------------:|----------:|----------:|---------------:|-----------:|----------------------:|---------------:|--------------------:|--------------------:|----------------:|----------------:|----------------------:|----------------------:|----------:|----------:|-----------:|-----------------:|-------------:|------------------------:|-----------------:|----------------------:|----------------------:|------------------:|------------------:|------------------------:|------------------------:|----------:|----------:|---------------:|-----------:|----------------------:|---------------:|--------------------:|--------------------:|----------------:|----------------:|----------------------:|----------------------:|----------:|----------:|-----------:|------------------:|--------------:|-------------------------:|------------------:|-----------------------:|-----------------------:|-------------------:|-------------------:|-------------------------:|-------------------------:|-----------:|-----------:|----------------:|------------:|-----------------------:|----------------:|---------------------:|---------------------:|-----------------:|-----------------:|-----------------------:|-----------------------:|----------:|----------:|-----------:|------------------:|--------------:|-------------------------:|------------------:|-----------------------:|-----------------------:|-------------------:|-------------------:|-------------------------:|-------------------------:|-----------:|-----------:|----------------:|------------:|-----------------------:|----------------:|---------------------:|---------------------:|-----------------:|-----------------:|-----------------------:|-----------------------:|---------:|---------:|-------:|----------:|----------:|
| 254261.6 | 4211270 |  1981.610 | 13.713031 |         4.836152 | -0.0006532 |     -0.0001182 |     -0.0005350 | -0.8709745 | -0.4901150 |      -11.9496506 |        -6.710705 | 1.872068 | 0.4429025 |         1981.683 |    13.278369 |                5.000895 |       -0.0000972 |             0.0002403 |            -0.0003375 |        -0.8590981 |        -0.5059774 |             -11.4257028 |               -6.691791 |  2.064680 | 0.4998173 |       3.130173 |  1.0335127 |              3.086304 |      0.0012735 |           0.0007820 |           0.0009248 |       0.0383560 |       0.0632483 |             1.1485910 |             0.8027909 | 0.6654974 | 0.1760111 | -0.0872438 |         1981.785 |    12.146732 |                7.700300 |       -0.0005357 |             0.0004232 |            -0.0009589 |        -0.7898716 |        -0.5664155 |              -9.9675833 |               -6.582414 |  2.333037 | 0.5448615 |       5.367796 |   2.629777 |              6.024163 |      0.0034881 |           0.0018882 |           0.0022296 |       0.1662535 |       0.1497663 |              3.093170 |              1.417825 | 0.8723594 | 0.1758332 | -0.2201196 |         1982.022 |    10.965063 |                8.366605 |       -0.0002665 |             0.0006450 |            -0.0009114 |        -0.5406073 |        -0.6087223 |              -7.1179222 |               -6.269690 |  2.244161 | 0.5090738 |       6.797439 |   3.379757 |              6.225783 |      0.0057246 |           0.0040034 |           0.0028513 |       0.5349155 |       0.2245785 |              6.087985 |              2.131548 | 0.9610642 | 0.1844364 | -0.5299986 |         1982.224 |    10.712748 |                8.559836 |        0.0000906 |             0.0007594 |            -0.0006688 |        -0.2760086 |        -0.5725025 |              -4.0612781 |               -5.849623 |  2.196398 | 0.4977898 |       7.394139 |   3.219312 |              6.008173 |      0.0057802 |           0.0041163 |           0.0029022 |       0.7207704 |       0.2799355 |              8.170098 |              2.805619 | 1.0006546 | 0.1986052 | -0.7874794 |          1983.219 |     10.193937 |                 8.283773 |        -0.0006227 |              0.0001124 |             -0.0007351 |          0.0194719 |         -0.4808405 |               -0.7604854 |                -4.385452 |   2.166422 |  0.5077887 |        9.185544 |    3.579978 |               6.040657 |       0.0055198 |            0.0038237 |            0.0027848 |        0.7942277 |        0.3710734 |               9.204818 |               3.487461 | 0.8736256 | 0.1920967 | -1.9473222 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |      NaN |      NaN |    NaN | 0.0012631 | 0.0003158 |
| 254262.5 | 4210646 |  1950.470 |  9.912198 |         5.537601 |  0.0026232 |      0.0025528 |      0.0000703 |  0.7729588 | -0.6268630 |        7.6910022 |        -6.181036 | 0.950337 | 0.1648604 |         1950.089 |     9.895582 |                3.867720 |        0.0017617 |             0.0021006 |            -0.0003389 |         0.7293449 |        -0.6304937 |               7.3567637 |               -6.079370 |  1.338211 | 0.3129493 |       2.268675 |  0.9818932 |              1.911031 |      0.0010410 |           0.0007097 |           0.0007018 |       0.1826580 |       0.1835575 |             2.3987532 |             1.3157071 | 0.3150602 | 0.1023118 |  0.4678890 |         1949.432 |     9.738675 |                3.976037 |        0.0013135 |             0.0017280 |            -0.0004145 |         0.6526975 |        -0.6290130 |               6.4656302 |               -5.973328 |  1.364365 | 0.3200382 |       3.934688 |   1.096268 |              1.873544 |      0.0012670 |           0.0010539 |           0.0008206 |       0.3312658 |       0.2573812 |              3.625255 |              2.275827 | 0.4513263 | 0.1332793 |  1.3252382 |         1948.666 |     9.604479 |                4.493677 |        0.0012742 |             0.0014916 |            -0.0002174 |         0.5814076 |        -0.6018455 |               5.6023839 |               -5.651684 |  1.347132 | 0.3159693 |       5.296040 |   1.399406 |              2.198534 |      0.0015949 |           0.0013567 |           0.0009141 |       0.4484710 |       0.3115568 |              4.684904 |              2.960944 | 0.4558434 | 0.1326649 |  2.2430519 |         1947.757 |     9.216747 |                5.426389 |        0.0006815 |             0.0009291 |            -0.0002476 |         0.5497256 |        -0.5487505 |               5.0905380 |               -4.888559 |  1.480344 | 0.3505203 |       6.369873 |   1.822887 |              2.752121 |      0.0029016 |           0.0024247 |           0.0011867 |       0.4877797 |       0.3995020 |              4.950153 |              3.746727 | 0.6004774 | 0.1635679 |  3.3728323 |          1946.364 |      8.267451 |                 6.788742 |        -0.0002585 |              0.0001347 |             -0.0003932 |          0.6002271 |         -0.4232695 |                5.0177675 |                -2.804953 |   1.610836 |  0.3795503 |        9.446887 |    2.988803 |               5.244349 |       0.0039799 |            0.0033390 |            0.0020126 |        0.4112044 |        0.5402399 |               4.220496 |               5.143875 | 0.8914608 | 0.2121961 |  4.6071992 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |      NaN |      NaN |    NaN | 0.0057449 | 0.0075172 |
| 254269.1 | 4211645 |  2015.906 | 11.360659 |         7.602951 |  0.0003760 |      0.0002771 |      0.0000989 |  0.0612487 | -0.9980820 |        0.6965038 |       -11.338843 | 1.143346 | 0.2750463 |         2015.854 |    10.495452 |                7.301916 |       -0.0000382 |            -0.0002666 |             0.0002284 |        -0.0108234 |        -0.9721438 |              -0.0601397 |              -10.249870 |  1.534582 | 0.3756965 |       2.541465 |  1.7074367 |              3.023076 |      0.0019520 |           0.0022463 |           0.0011599 |       0.2269458 |       0.0418870 |             2.0704124 |             1.9055884 | 0.5495349 | 0.1651123 |  0.0641070 |         2015.878 |     9.692884 |                6.727752 |        0.0001912 |             0.0000025 |             0.0001886 |        -0.0860654 |        -0.7998923 |              -0.9062279 |               -8.147747 |  1.464120 | 0.3600917 |       4.062199 |   2.482098 |              3.858645 |      0.0025398 |           0.0034538 |           0.0018810 |       0.5073005 |       0.2957152 |              4.238674 |              3.835309 | 0.4805955 | 0.1458878 |  0.0181719 |         2015.781 |     9.758282 |                7.342618 |        0.0006537 |             0.0002348 |             0.0004189 |        -0.0735984 |        -0.5794520 |              -0.9397725 |               -6.274251 |  1.464234 | 0.3557911 |       5.141474 |   2.675841 |              4.442736 |      0.0028261 |           0.0032081 |           0.0019896 |       0.6360680 |       0.5010540 |              6.059131 |              5.010459 | 0.4863672 | 0.1503251 |  0.1512862 |         2015.428 |    10.110068 |                7.557088 |        0.0005121 |             0.0001230 |             0.0003891 |        -0.0819901 |        -0.4188874 |              -0.8800428 |               -4.753003 |  1.495424 | 0.3602871 |       6.083861 |   2.906325 |              4.423169 |      0.0029230 |           0.0030535 |           0.0018389 |       0.6763182 |       0.5999668 |              7.353033 |              5.758430 | 0.4880002 | 0.1481220 |  0.6424035 |          2013.974 |     11.067215 |                 8.054791 |         0.0000719 |              0.0000917 |             -0.0000198 |         -0.1056094 |         -0.2939010 |               -1.3918060 |                -2.573399 |   1.653608 |  0.3950393 |        9.996587 |    3.765331 |               5.434985 |       0.0045584 |            0.0035310 |            0.0024847 |        0.7324326 |        0.6053891 |               9.266498 |               6.505241 | 0.6301132 | 0.1664494 |  2.4341684 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |      NaN |      NaN |    NaN | 0.0125898 | 0.0057296 |
| 254288.0 | 4210645 |  1946.073 | 10.564654 |         4.417495 | -0.0002761 |      0.0009774 |     -0.0012535 |  0.9076828 | -0.4194873 |        9.5877494 |        -4.435286 | 1.200804 | 0.2794998 |         1946.097 |    10.259158 |                3.670922 |        0.0004342 |             0.0010890 |            -0.0006548 |         0.8969835 |        -0.4292127 |               9.2096136 |               -4.389629 |  1.214710 | 0.2846959 |       2.398183 |  0.8063123 |              1.466709 |      0.0010973 |           0.0006197 |           0.0006473 |       0.0476514 |       0.0936500 |             0.9661455 |             0.9369128 | 0.2918260 | 0.0906148 | -0.0275695 |         1945.917 |     9.762910 |                4.014085 |        0.0006554 |             0.0011192 |            -0.0004639 |         0.8539999 |        -0.4566873 |               8.4097021 |               -4.356122 |  1.284320 | 0.3003475 |       4.141065 |   1.082522 |              1.915924 |      0.0012422 |           0.0011130 |           0.0007877 |       0.1485868 |       0.1959145 |              1.989985 |              1.653142 | 0.3399680 | 0.1015795 |  0.2353106 |         1945.585 |     9.253997 |                4.661162 |        0.0005960 |             0.0010588 |            -0.0004628 |         0.7867989 |        -0.4695633 |               7.3468761 |               -4.234184 |  1.361019 | 0.3219995 |       5.492108 |   1.541400 |              2.166626 |      0.0016508 |           0.0014778 |           0.0008275 |       0.2812076 |       0.2829367 |              3.118125 |              2.517851 | 0.4277034 | 0.1275102 |  0.6621808 |         1945.154 |     8.789118 |                5.021721 |        0.0006977 |             0.0010748 |            -0.0003771 |         0.7056091 |        -0.4723078 |               6.1699694 |               -4.018459 |  1.454046 | 0.3431825 |       6.514370 |   1.994746 |              2.516838 |      0.0018984 |           0.0015493 |           0.0008953 |       0.3874964 |       0.3581299 |              3.988227 |              3.332455 | 0.5159743 | 0.1446458 |  1.2162529 |          1944.157 |      7.840729 |                 6.149886 |        -0.0002144 |              0.0001803 |             -0.0003947 |          0.5661268 |         -0.4745410 |                4.6103570 |                -2.988381 |   1.481068 |  0.3470639 |        9.231529 |    3.045412 |               5.114635 |       0.0035347 |            0.0029392 |            0.0018642 |        0.4257046 |        0.5229101 |               4.245648 |               4.750717 | 0.9616770 | 0.2299496 |  2.2678369 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |      NaN |      NaN |    NaN | 0.0044411 | 0.0011103 |
| 254288.8 | 4211215 |  1979.024 | 13.033809 |         4.717597 |  0.0023390 |      0.0027178 |     -0.0003788 | -0.6587729 | -0.7473569 |       -8.6120304 |        -9.719215 | 1.532461 | 0.3927215 |         1978.698 |    12.919998 |                4.554495 |        0.0015331 |             0.0021250 |            -0.0005919 |        -0.6188551 |        -0.7514454 |              -8.1392068 |               -9.596307 |  1.435098 | 0.3539592 |       2.983899 |  1.0425038 |              1.857948 |      0.0008667 |           0.0010100 |           0.0006633 |       0.1776541 |       0.1271919 |             2.7871990 |             1.0985602 | 0.4540820 | 0.1501357 |  0.3981930 |         1978.136 |    12.109388 |                7.117271 |        0.0009778 |             0.0019806 |            -0.0010027 |        -0.4886927 |        -0.7727529 |              -6.2189171 |               -9.179339 |  1.804128 | 0.4215947 |       5.174129 |   2.011347 |              6.117043 |      0.0037841 |           0.0020642 |           0.0022813 |       0.3659859 |       0.1642128 |              4.914267 |              1.814923 | 0.8287974 | 0.1701015 |  1.1338561 |         1977.588 |    10.945342 |                8.302665 |        0.0003929 |             0.0015603 |            -0.0011673 |        -0.2550540 |        -0.7463884 |              -3.4677430 |               -7.996499 |  2.000389 | 0.4545434 |       6.399788 |   3.024382 |              5.958027 |      0.0046163 |           0.0029605 |           0.0025629 |       0.5806652 |       0.2067502 |              6.809877 |              2.596817 | 0.9628345 | 0.1784885 |  1.7685225 |         1977.226 |    10.139661 |                8.500556 |       -0.0003400 |             0.0007486 |            -0.0010886 |        -0.1297924 |        -0.7009772 |              -1.9701576 |               -6.839830 |  2.121444 | 0.4838760 |       7.074201 |   3.313769 |              5.790195 |      0.0053343 |           0.0038859 |           0.0026216 |       0.6527983 |       0.2591712 |              7.421959 |              2.869355 | 1.0117353 | 0.1944909 |  2.1665952 |          1979.303 |     10.107605 |                 7.605051 |        -0.0004125 |              0.0002886 |             -0.0007011 |         -0.0572587 |         -0.5309603 |               -1.4387831 |                -4.835391 |   2.085570 |  0.4892479 |        9.830509 |    3.523873 |               5.763515 |       0.0048081 |            0.0033624 |            0.0024597 |        0.7792537 |        0.3281500 |               8.975127 |               2.930723 | 0.8727765 | 0.1893673 | -0.9602529 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |      NaN |      NaN |    NaN | 0.0042577 | 0.0142399 |
| 254290.0 | 4210850 |  1959.791 |  7.278730 |         9.417492 |  0.0031509 |      0.0013961 |      0.0017548 | -0.2686590 | -0.9492942 |       -2.1088601 |        -6.873158 | 1.530910 | 0.3593948 |         1959.313 |     7.209230 |                7.046494 |        0.0024009 |             0.0013282 |             0.0010727 |        -0.1172325 |        -0.8838886 |              -1.4680066 |               -6.396176 |  1.593134 | 0.3841811 |       1.629686 |  1.8667068 |              3.062235 |      0.0014627 |           0.0016500 |           0.0011087 |       0.4211348 |       0.0956179 |             2.8818207 |             1.8068683 | 0.5706633 | 0.1482188 |  0.5868906 |         1958.422 |     7.646554 |                6.611842 |        0.0017697 |             0.0009048 |             0.0008648 |         0.0370689 |        -0.7874292 |              -0.3649251 |               -6.152959 |  1.735174 | 0.4234442 |       2.839513 |   2.177274 |              2.788832 |      0.0019870 |           0.0019687 |           0.0009802 |       0.5627181 |       0.2351393 |              4.132931 |              2.765202 | 0.6056337 | 0.1527924 |  1.7586532 |         1957.424 |     7.948945 |                7.168607 |        0.0010527 |             0.0003916 |             0.0006611 |         0.1897994 |        -0.6936279 |               1.1662865 |               -5.519507 |  1.831219 | 0.4478300 |       3.905047 |   2.337499 |              3.698836 |      0.0032549 |           0.0027700 |           0.0014511 |       0.5781731 |       0.3842246 |              4.747839 |              3.760410 | 0.6573314 | 0.1621630 |  2.9426691 |         1956.507 |     8.242990 |                7.442571 |        0.0005428 |             0.0000145 |             0.0005283 |         0.3521928 |        -0.5519230 |               2.8900951 |               -4.250480 |  1.918723 | 0.4688302 |       4.688016 |   2.645136 |              4.564167 |      0.0036736 |           0.0031504 |           0.0019145 |       0.5625148 |       0.5068311 |              5.108430 |              4.749298 | 0.7141760 | 0.1753271 |  4.0422615 |          1955.361 |      7.422135 |                 6.909944 |        -0.0001914 |             -0.0000880 |             -0.0001034 |          0.5901128 |         -0.3596705 |                4.6057681 |                -2.105288 |   1.899388 |  0.4449421 |        7.790299 |    2.955246 |               4.737072 |       0.0034215 |            0.0026107 |            0.0021289 |        0.4526467 |        0.5638004 |               4.350798 |               4.389752 | 1.0125595 | 0.2264864 |  4.8458056 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |      NaN |      NaN |    NaN | 0.0038706 | 0.0009677 |

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
`mod_vi()` twice, once for each radius

``` r
mod125 <- mod_vi(preds$pred_pts, 125, cross_validate = TRUE, tune = FALSE, num_cores = 4)
mod250<- mod_vi(preds$pred_pts, 250, cross_validate = TRUE, tune = FALSE, num_cores = 4)
```

Then we can look at model performance using the `perf_mets` object
returned by `mod_vi()`:

``` r
mod125$perf_mets
mod250$perf_mets
```

<div class="kable-table">

|        r2 |      rmse |     nrmse | VI Radius (m) |
|----------:|----------:|----------:|--------------:|
| 0.7480373 | 0.0873732 | 0.1586089 |           125 |
| 0.5639866 | 0.0765572 | 0.1795530 |           250 |

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
predictors at 32 pixel radii tends to be slow.

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
