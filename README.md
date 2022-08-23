# README

## Code

The files here may be used to reproduce the analysis performed in this project (although this analysis -- particularly calculation of connectivity metrics -- will take some time).

### Library code

`factored.R`: This contains helper functions used in the calculation of connectivity metrics.

`extern/`: This folder must be populated by a file containing the `corvif` function from Mixed effects models and extensions in ecology with R. (2009). Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer. This is used to perform colinearity testing.

### Raster and data outputs

`background-metrics.R`: This produces the primary connectivity metrics for Western and Central Europe from the CSIRO rasters.  
`build-eu-rasters.R`: This produces the rasters of secondary connectivity metrics for Western and Central Europe after running `background-metrics.R`.

### Core analyses

`conn-data-eu.R`: This calculates the primary connectivity metrics for the sites used in the analysis, from the PREDICTS data (to get the locations) and the CSIRO rasters (to calculate the connectivity metrics).  
`get-abund-and-comp-diss.R`: This uses the PREDICTS data to calculate community abundance and compositional similarity metrics.  
`europe.R`: This uses the outputs above to calculate derived connectivity metrics, and fit the models and produce outputs as presented in Sections 2.6 and 3.

`connectivity-analysis.R`: This performs a standalone analysis on the appropriate threshold to consider areas as patches, as described in Section 2.4.

`metrics-figure.R`, `new-forest-metrics.R`: These calculate connectivity metrics for the New Forest from the PREDICTS data and CSIRO rasters, in order to produce the images presented in Figure 2.

## Data

The CSIRO raster data can be accessed directly from https://data.csiro.au/collection/csiro:15276v3. The PREDICTS data that I used was from an unpublished extract which is not available for general use. The analysis may be performed on a public PREDICTS version (e.g. https://data.nhm.ac.uk/dataset/the-2016-release-of-the-predicts-database), however the results presented here incorporate more recent data and will not be precisely replicable.

## Results

This contains the output CSV files that are used in the final report. These are produced as part of running the analyses described above.
