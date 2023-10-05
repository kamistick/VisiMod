VisiMod: Map Modeled Visibility Index Across Wildland Landscapes
================
Katherine Mistick
2023-10-05

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
al. 2023](https://content.csbs.utah.edu/~pdennison/reprints/denn/2023_Mistick_etal_IJGIS.pdf)
is the original research that led to the development of VisiMod.

Instead of calculating VI by running viewsheds from every pixel in a
rasterized study area, VisiMod allows users to build predictive models
based on a sample of points (at which viewsheds are run using
[terra::viewshed()](https://rdrr.io/cran/terra/man/viewshed.html)) to
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
| 254582.5 | 4211050 |
| 255036.3 | 4210476 |
| 255212.7 | 4211446 |
| 255056.1 | 4211712 |
| 254713.1 | 4210792 |
| 255372.8 | 4210609 |

</div>

While the output of `gen_pts()` is a dataframe, we can vectorize it
using `v <- vect(cbind(gp$x, gp$y))` to plot the output, overlaid on our
DSM. You can see how max_vi_rad ensured that no points felt within 250 m
of the study area’s edge:

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

|        x |       y | elevation |     slope | slope_derivative |  curvature | curvature_plan | curvature_prof | aspect_sin | aspect_cos | slope_aspect_sin | slope_aspect_cos |        ch |        cc | elevation_mean_2 | slope_mean_2 | slope_derivative_mean_2 | curvature_mean_2 | curvature_plan_mean_2 | curvature_prof_mean_2 | aspect_sin_mean_2 | aspect_cos_mean_2 | slope_aspect_sin_mean_2 | slope_aspect_cos_mean_2 | ch_mean_2 | cc_mean_2 | elevation_sd_2 | slope_sd_2 | slope_derivative_sd_2 | curvature_sd_2 | curvature_plan_sd_2 | curvature_prof_sd_2 | aspect_sin_sd_2 | aspect_cos_sd_2 | slope_aspect_sin_sd_2 | slope_aspect_cos_sd_2 |   ch_sd_2 |   cc_sd_2 |      tpi_2 | elevation_mean_4 | slope_mean_4 | slope_derivative_mean_4 | curvature_mean_4 | curvature_plan_mean_4 | curvature_prof_mean_4 | aspect_sin_mean_4 | aspect_cos_mean_4 | slope_aspect_sin_mean_4 | slope_aspect_cos_mean_4 | ch_mean_4 | cc_mean_4 | elevation_sd_4 | slope_sd_4 | slope_derivative_sd_4 | curvature_sd_4 | curvature_plan_sd_4 | curvature_prof_sd_4 | aspect_sin_sd_4 | aspect_cos_sd_4 | slope_aspect_sin_sd_4 | slope_aspect_cos_sd_4 |   ch_sd_4 |   cc_sd_4 |      tpi_4 | elevation_mean_6 | slope_mean_6 | slope_derivative_mean_6 | curvature_mean_6 | curvature_plan_mean_6 | curvature_prof_mean_6 | aspect_sin_mean_6 | aspect_cos_mean_6 | slope_aspect_sin_mean_6 | slope_aspect_cos_mean_6 | ch_mean_6 | cc_mean_6 | elevation_sd_6 | slope_sd_6 | slope_derivative_sd_6 | curvature_sd_6 | curvature_plan_sd_6 | curvature_prof_sd_6 | aspect_sin_sd_6 | aspect_cos_sd_6 | slope_aspect_sin_sd_6 | slope_aspect_cos_sd_6 |   ch_sd_6 |   cc_sd_6 |      tpi_6 | elevation_mean_8 | slope_mean_8 | slope_derivative_mean_8 | curvature_mean_8 | curvature_plan_mean_8 | curvature_prof_mean_8 | aspect_sin_mean_8 | aspect_cos_mean_8 | slope_aspect_sin_mean_8 | slope_aspect_cos_mean_8 | ch_mean_8 | cc_mean_8 | elevation_sd_8 | slope_sd_8 | slope_derivative_sd_8 | curvature_sd_8 | curvature_plan_sd_8 | curvature_prof_sd_8 | aspect_sin_sd_8 | aspect_cos_sd_8 | slope_aspect_sin_sd_8 | slope_aspect_cos_sd_8 |   ch_sd_8 |   cc_sd_8 |      tpi_8 | elevation_mean_16 | slope_mean_16 | slope_derivative_mean_16 | curvature_mean_16 | curvature_plan_mean_16 | curvature_prof_mean_16 | aspect_sin_mean_16 | aspect_cos_mean_16 | slope_aspect_sin_mean_16 | slope_aspect_cos_mean_16 | ch_mean_16 | cc_mean_16 | elevation_sd_16 | slope_sd_16 | slope_derivative_sd_16 | curvature_sd_16 | curvature_plan_sd_16 | curvature_prof_sd_16 | aspect_sin_sd_16 | aspect_cos_sd_16 | slope_aspect_sin_sd_16 | slope_aspect_cos_sd_16 |  ch_sd_16 |  cc_sd_16 |     tpi_16 | elevation_mean_32 | slope_mean_32 | slope_derivative_mean_32 | curvature_mean_32 | curvature_plan_mean_32 | curvature_prof_mean_32 | aspect_sin_mean_32 | aspect_cos_mean_32 | slope_aspect_sin_mean_32 | slope_aspect_cos_mean_32 | ch_mean_32 | cc_mean_32 | elevation_sd_32 | slope_sd_32 | slope_derivative_sd_32 | curvature_sd_32 | curvature_plan_sd_32 | curvature_prof_sd_32 | aspect_sin_sd_32 | aspect_cos_sd_32 | slope_aspect_sin_sd_32 | slope_aspect_cos_sd_32 | ch_sd_32 | cc_sd_32 | tpi_32 |    vi_125 |    vi_250 |
|---------:|--------:|----------:|----------:|-----------------:|-----------:|---------------:|---------------:|-----------:|-----------:|-----------------:|-----------------:|----------:|----------:|-----------------:|-------------:|------------------------:|-----------------:|----------------------:|----------------------:|------------------:|------------------:|------------------------:|------------------------:|----------:|----------:|---------------:|-----------:|----------------------:|---------------:|--------------------:|--------------------:|----------------:|----------------:|----------------------:|----------------------:|----------:|----------:|-----------:|-----------------:|-------------:|------------------------:|-----------------:|----------------------:|----------------------:|------------------:|------------------:|------------------------:|------------------------:|----------:|----------:|---------------:|-----------:|----------------------:|---------------:|--------------------:|--------------------:|----------------:|----------------:|----------------------:|----------------------:|----------:|----------:|-----------:|-----------------:|-------------:|------------------------:|-----------------:|----------------------:|----------------------:|------------------:|------------------:|------------------------:|------------------------:|----------:|----------:|---------------:|-----------:|----------------------:|---------------:|--------------------:|--------------------:|----------------:|----------------:|----------------------:|----------------------:|----------:|----------:|-----------:|-----------------:|-------------:|------------------------:|-----------------:|----------------------:|----------------------:|------------------:|------------------:|------------------------:|------------------------:|----------:|----------:|---------------:|-----------:|----------------------:|---------------:|--------------------:|--------------------:|----------------:|----------------:|----------------------:|----------------------:|----------:|----------:|-----------:|------------------:|--------------:|-------------------------:|------------------:|-----------------------:|-----------------------:|-------------------:|-------------------:|-------------------------:|-------------------------:|-----------:|-----------:|----------------:|------------:|-----------------------:|----------------:|---------------------:|---------------------:|-----------------:|-----------------:|-----------------------:|-----------------------:|----------:|----------:|-----------:|------------------:|--------------:|-------------------------:|------------------:|-----------------------:|-----------------------:|-------------------:|-------------------:|-------------------------:|-------------------------:|-----------:|-----------:|----------------:|------------:|-----------------------:|----------------:|---------------------:|---------------------:|-----------------:|-----------------:|-----------------------:|-----------------------:|---------:|---------:|-------:|----------:|----------:|
| 254271.7 | 4211030 |  1959.638 |  4.221355 |         4.122323 | -0.0011012 |     -0.0008054 |     -0.0002958 |  0.6889122 | -0.7143504 |         2.946805 |        -2.975908 | 3.2876549 | 0.5800823 |         1959.871 |     4.593181 |                4.147725 |       -0.0016207 |            -0.0012238 |            -0.0003969 |         0.5617132 |        -0.6830530 |                2.724503 |               -2.990823 |  2.939082 | 0.5938527 |      0.9927259 |  0.9341392 |              2.270315 |      0.0018823 |           0.0016625 |           0.0010640 |       0.3638534 |       0.2638702 |              2.034086 |             1.0438598 | 0.9568660 | 0.1830849 | -0.2914080 |         1960.470 |     5.132715 |                4.457108 |       -0.0012059 |            -0.0008326 |            -0.0003733 |         0.4976504 |        -0.5666009 |               2.7960534 |               -2.449015 |  3.015147 | 0.6312185 |       1.786613 |   1.371932 |              3.601225 |      0.0028182 |           0.0021415 |           0.0014566 |       0.4647220 |       0.4581859 |              2.854801 |              2.456786 | 0.9060480 | 0.1666553 | -1.0938664 |         1961.160 |     5.526597 |                4.671471 |       -0.0008657 |            -0.0006275 |            -0.0002381 |         0.4838668 |        -0.4921754 |                3.008975 |               -2.065476 |  2.820787 | 0.6099957 |       2.519823 |   1.682495 |              4.218806 |      0.0030875 |           0.0021605 |           0.0018757 |       0.4690483 |       0.5494793 |              3.162432 |              3.160342 | 0.8899469 | 0.1600639 | -1.9195281 |         1961.866 |     5.690128 |                4.955098 |       -0.0007318 |            -0.0005928 |            -0.0001390 |         0.4802286 |        -0.4526504 |                3.091072 |               -1.949622 |  2.697694 | 0.5989059 |       3.238225 |   1.862014 |              4.704021 |      0.0033157 |           0.0023113 |           0.0019999 |       0.4793075 |       0.5778699 |              3.355122 |              3.349350 | 0.8394245 | 0.1542461 | -2.7782842 |          1965.133 |      6.999870 |                 6.084065 |        -0.0004429 |             -0.0000707 |             -0.0003722 |          0.3715443 |         -0.5226269 |                2.5875140 |                -3.253435 |   2.293302 |  0.5371017 |        6.609104 |    2.702185 |               4.642562 |       0.0032528 |            0.0024584 |            0.0019512 |        0.5901458 |        0.4907261 |               5.066602 |               3.651840 | 0.8329517 | 0.1700450 | -6.6103535 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |      NaN |      NaN |    NaN | 0.0006315 | 0.0001579 |
| 254281.9 | 4210598 |  1943.241 |  7.115165 |         6.637297 |  0.0010836 |      0.0021174 |     -0.0010339 |  0.8437862 | -0.5136694 |         6.034255 |        -3.615499 | 1.3444394 | 0.3245889 |         1943.053 |     7.606914 |                5.331801 |        0.0012376 |             0.0019330 |            -0.0006954 |         0.6742495 |        -0.6046552 |                5.230100 |               -4.503368 |  1.215339 | 0.2690123 |      1.6688722 |  1.2800515 |              1.669930 |      0.0016736 |           0.0021607 |           0.0009840 |       0.3056530 |       0.2741293 |              2.710507 |             1.9524138 | 0.2872554 | 0.0960258 |  0.2326378 |         1942.625 |     8.200677 |                4.229713 |        0.0008470 |             0.0012030 |            -0.0003560 |         0.5485944 |        -0.6427355 |               4.5805332 |               -5.216003 |  1.343436 | 0.3084576 |       3.092116 |   1.346481 |              1.841670 |      0.0019797 |           0.0020037 |           0.0008397 |       0.4208672 |       0.3275000 |              3.694603 |              2.664752 | 0.3961876 | 0.1294511 |  0.8034049 |         1942.169 |     8.220738 |                4.337358 |        0.0003569 |             0.0007845 |            -0.0004276 |         0.5409534 |        -0.6429815 |                4.525409 |               -5.232211 |  1.461047 | 0.3395254 |       4.387478 |   1.507595 |              2.309549 |      0.0023041 |           0.0021714 |           0.0009081 |       0.4250948 |       0.3379219 |              3.777733 |              2.788801 | 0.5786360 | 0.1591472 |  1.3364323 |         1941.833 |     8.184162 |                5.456676 |       -0.0001652 |             0.0005059 |            -0.0006711 |         0.5789958 |        -0.5794389 |                4.770854 |               -4.590801 |  1.641242 | 0.3812361 |       5.490966 |   1.849464 |              4.188937 |      0.0028290 |           0.0023681 |           0.0015487 |       0.4257354 |       0.3849993 |              3.902569 |              3.365993 | 0.7465936 | 0.1803469 |  1.7208462 |          1942.285 |      7.794143 |                 6.132432 |        -0.0001600 |              0.0002042 |             -0.0003642 |          0.6261750 |         -0.4662025 |                4.9558210 |                -2.886799 |   1.492119 |  0.3492344 |        9.258957 |    3.083126 |               5.188091 |       0.0038262 |            0.0031902 |            0.0018894 |        0.3831918 |        0.4939757 |               3.926416 |               4.687242 | 0.9301047 | 0.2221841 |  0.8143859 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |      NaN |      NaN |    NaN | 0.0014464 | 0.0003616 |
| 254295.1 | 4210821 |  1955.830 | 10.442064 |         8.608534 |  0.0013048 |     -0.0003600 |      0.0016647 | -0.4266753 | -0.9042800 |        -4.461195 |        -9.439907 | 1.3551765 | 0.3328561 |         1955.554 |     9.625867 |                6.732131 |        0.0019128 |             0.0009341 |             0.0009787 |        -0.2752223 |        -0.9325361 |               -2.982879 |               -8.945882 |  1.692290 | 0.4092447 |      2.3166826 |  1.8233872 |              3.248226 |      0.0022162 |           0.0018630 |           0.0010227 |       0.2288881 |       0.0372436 |              2.151502 |             1.5601698 | 0.6033248 | 0.1533761 |  0.3448873 |         1954.831 |     8.946964 |                6.710862 |        0.0013544 |             0.0011419 |             0.0002126 |        -0.0405565 |        -0.8677622 |              -0.8551166 |               -7.848357 |  1.736091 | 0.4115175 |       3.703642 |   1.963248 |              3.828527 |      0.0029093 |           0.0025023 |           0.0014620 |       0.4733368 |       0.1557835 |              3.994909 |              2.413214 | 0.7572968 | 0.1808614 |  1.3142379 |         1954.129 |     8.533456 |                6.679651 |       -0.0000234 |             0.0001715 |            -0.0001949 |         0.1924761 |        -0.7324317 |                1.272993 |               -6.332207 |  1.788169 | 0.4233890 |       4.566920 |   2.117523 |              4.333696 |      0.0039771 |           0.0035707 |           0.0017959 |       0.5601217 |       0.3377943 |              4.804978 |              3.550784 | 0.8017725 | 0.1873394 |  2.1370170 |         1953.822 |     8.709974 |                6.849710 |       -0.0000761 |             0.0000700 |            -0.0001461 |         0.3747270 |        -0.5718823 |                3.200628 |               -4.779197 |  1.881083 | 0.4481595 |       5.159433 |   2.199482 |              4.428872 |      0.0037788 |           0.0033042 |           0.0018060 |       0.5556884 |       0.4752982 |              5.099155 |              4.667396 | 0.7938012 | 0.1851004 |  2.4238975 |          1953.815 |      7.437386 |                 6.835801 |        -0.0001038 |              0.0000004 |             -0.0001042 |          0.5855522 |         -0.3916203 |                4.6250591 |                -2.305782 |   1.751915 |  0.4116467 |        8.137242 |    2.953063 |               4.775697 |       0.0033153 |            0.0025144 |            0.0020974 |        0.4582636 |        0.5424991 |               4.370755 |               4.272697 | 1.0249498 | 0.2339000 |  2.0244643 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |      NaN |      NaN |    NaN | 0.0207385 | 0.0073644 |
| 254299.4 | 4211153 |  1969.582 |  7.421292 |         9.116500 | -0.0016240 |     -0.0000038 |     -0.0016203 | -0.1108791 | -0.9931425 |        -0.837068 |        -7.369687 | 1.4654485 | 0.3628580 |         1969.817 |     7.142751 |                7.042369 |       -0.0003698 |             0.0007182 |            -0.0010880 |        -0.1983374 |        -0.9535375 |               -1.461875 |               -6.804049 |  1.666637 | 0.4192372 |      1.6714354 |  1.8102503 |              2.131571 |      0.0020924 |           0.0018929 |           0.0007682 |       0.2048906 |       0.0755407 |              1.482961 |             1.8503119 | 0.4372035 | 0.1043697 | -0.2882648 |         1969.906 |     7.257073 |                6.279982 |        0.0000084 |             0.0008819 |            -0.0008736 |        -0.2564781 |        -0.8930339 |              -2.0457797 |               -6.378325 |  1.865898 | 0.4479359 |       2.954722 |   2.515262 |              3.909337 |      0.0033314 |           0.0024075 |           0.0017782 |       0.3351422 |       0.1560748 |              2.935499 |              2.344660 | 0.6747768 | 0.1439927 | -0.3635382 |         1969.963 |     7.398166 |                6.750455 |       -0.0007850 |             0.0002056 |            -0.0009906 |        -0.1856729 |        -0.8649181 |               -1.723625 |               -6.263411 |  2.059583 | 0.4772508 |       4.097122 |   2.931293 |              4.710007 |      0.0038639 |           0.0030386 |           0.0018343 |       0.4359745 |       0.1646217 |              3.819190 |              2.563528 | 0.8760061 | 0.1722810 | -0.4079831 |         1970.369 |     7.885889 |                6.862033 |       -0.0008793 |             0.0001057 |            -0.0009850 |        -0.1428348 |        -0.7950964 |               -1.637201 |               -6.022931 |  2.174966 | 0.5051670 |       5.160552 |   3.153047 |              4.906298 |      0.0038603 |           0.0030755 |           0.0018254 |       0.5468169 |       0.2184484 |              5.108544 |              2.651356 | 0.9067271 | 0.1824048 | -0.9569505 |          1973.962 |      9.323542 |                 6.373364 |        -0.0004350 |              0.0002730 |             -0.0007080 |          0.0083324 |         -0.6043115 |               -0.7621786 |                -5.230570 |   2.146256 |  0.5012109 |        9.596091 |    3.482380 |               5.013113 |       0.0039593 |            0.0028962 |            0.0019832 |        0.7327844 |        0.3126769 |               7.906498 |               2.930277 | 0.8988101 | 0.1867564 | -5.5974319 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |      NaN |      NaN |    NaN | 0.0036262 | 0.0009065 |
| 254307.2 | 4210779 |  1948.948 | 11.168335 |         4.629034 |  0.0008950 |      0.0013322 |     -0.0004372 | -0.2377718 | -0.9696935 |        -2.658609 |       -10.829047 | 0.8642901 | 0.1732117 |         1948.843 |     9.911809 |                7.006697 |       -0.0003037 |             0.0008001 |            -0.0011038 |        -0.0941628 |        -0.9435989 |               -1.292122 |               -9.489318 |  1.585555 | 0.3598417 |      2.4579332 |  1.7843311 |              5.378580 |      0.0033175 |           0.0024057 |           0.0017357 |       0.2704044 |       0.1072127 |              2.251121 |             1.9901485 | 0.8504192 | 0.2079411 |  0.1267247 |         1949.021 |     9.401572 |                7.864414 |       -0.0007788 |             0.0003715 |            -0.0011504 |         0.2070630 |        -0.7199053 |               1.5664048 |               -6.793855 |  1.872268 | 0.4194728 |       3.653154 |   2.183039 |              5.796325 |      0.0042737 |           0.0039821 |           0.0020884 |       0.4741528 |       0.4593620 |              4.493780 |              4.884450 | 0.9978422 | 0.2155005 | -0.1507315 |         1949.408 |     9.375200 |                7.515446 |       -0.0004216 |             0.0004134 |            -0.0008350 |         0.3849648 |        -0.5565665 |                3.482946 |               -4.909971 |  1.871771 | 0.4279193 |       4.499274 |   2.250949 |              5.115348 |      0.0037510 |           0.0033950 |           0.0018814 |       0.5015713 |       0.5396218 |              5.108297 |              5.541101 | 0.9377625 | 0.2090851 | -0.6663880 |         1949.793 |     9.076773 |                7.019383 |       -0.0003964 |             0.0002670 |            -0.0006634 |         0.4822833 |        -0.4769100 |                4.350608 |               -3.882665 |  1.760068 | 0.4096334 |       5.353524 |   2.318753 |              4.513125 |      0.0033855 |           0.0029940 |           0.0017027 |       0.5009484 |       0.5385455 |              5.096427 |              5.282957 | 0.9336918 | 0.2161355 | -1.1511951 |          1951.011 |      7.424007 |                 6.056871 |         0.0001109 |              0.0001476 |             -0.0000366 |          0.5434428 |         -0.4601410 |                4.3117011 |                -2.832592 |   1.532983 |  0.3626837 |        8.510274 |    2.849928 |               4.532063 |       0.0029436 |            0.0023642 |            0.0018023 |        0.4708753 |        0.5212007 |               4.345496 |               4.214967 | 1.0085331 | 0.2392857 | -2.4843554 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |      NaN |      NaN |    NaN | 0.0214719 | 0.0068093 |
| 254307.6 | 4211353 |  2000.463 | 11.771687 |        22.659692 |  0.0061763 |      0.0037515 |      0.0024248 | -0.7626857 | -0.5368374 |        -9.920223 |        -5.443765 | 1.7624618 | 0.4399743 |         1999.344 |    12.044555 |               13.795276 |        0.0055454 |             0.0039763 |             0.0015691 |        -0.4116892 |        -0.5079884 |               -6.284386 |               -5.116242 |  1.339874 | 0.3031138 |      2.6024050 |  3.9959312 |              9.711564 |      0.0066433 |           0.0056738 |           0.0029759 |       0.6952011 |       0.2694369 |              9.519159 |             0.6841081 | 0.5411269 | 0.1602291 |  1.3881510 |         1997.345 |    12.787950 |               10.886973 |        0.0031753 |             0.0022011 |             0.0009742 |        -0.2450999 |        -0.4420487 |              -3.7709664 |               -5.063902 |  1.403264 | 0.3226844 |       4.062632 |   3.128470 |              8.350836 |      0.0056835 |           0.0045423 |           0.0026718 |       0.8353931 |       0.2127229 |             11.468116 |              1.105572 | 0.5164313 | 0.1524034 |  3.9922669 |         1995.447 |    12.663985 |               10.289995 |        0.0015954 |             0.0012112 |             0.0003842 |        -0.2044990 |        -0.4279544 |               -2.984862 |               -4.980911 |  1.510444 | 0.3535340 |       5.321037 |   2.791682 |              7.471505 |      0.0057554 |           0.0039485 |           0.0029645 |       0.8615095 |       0.1892434 |             11.558836 |              1.161404 | 0.5598524 | 0.1605119 |  6.1799014 |         1993.972 |    12.204257 |               10.812641 |        0.0000866 |             0.0003479 |            -0.0002614 |        -0.2753796 |        -0.4566493 |               -3.706850 |               -5.097639 |  1.639130 | 0.3889270 |       6.350941 |   2.964521 |              7.134843 |      0.0063087 |           0.0041456 |           0.0031151 |       0.8240367 |       0.1990698 |             10.813554 |              1.254026 | 0.5884512 | 0.1633286 |  7.8501190 |          1992.774 |     11.184769 |                10.154253 |         0.0000868 |              0.0001751 |             -0.0000883 |         -0.2314306 |         -0.5487814 |               -3.6488962 |                -5.594023 |   1.834749 |  0.4275424 |       11.023378 |    3.815090 |               7.301304 |       0.0062344 |            0.0041846 |            0.0032121 |        0.7587021 |        0.2652723 |               9.459487 |               2.374737 | 0.8341993 | 0.1949997 |  8.1370693 |               NaN |           NaN |                      NaN |               NaN |                    NaN |                    NaN |                NaN |                NaN |                      NaN |                      NaN |        NaN |        NaN |             NaN |         NaN |                    NaN |             NaN |                  NaN |                  NaN |              NaN |              NaN |                    NaN |                    NaN |      NaN |      NaN |    NaN | 0.0003463 | 0.0000866 |

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
| 0.8545271 | 0.0667785 | 0.1112801 |           125 |
| 0.7767998 | 0.0508414 | 0.1158244 |           250 |

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
