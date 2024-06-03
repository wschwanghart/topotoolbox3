# ToDo-List

## GRIDobj

- [ ] GRIDOBJ              : Create instance of a GRIDobj
Changed refmat to wf (worldfile). Referencing matrices are no longer supported by the mapping toolbox.
- [ ] GRIDOBJ2ASCII        : write/export GRIDobj to ESRI ArcGIS ASCII file
- [ ] GRIDOBJ2GEOTIFF      : Exports an instance of GRIDobj to a geotiff file
- [ ] GRIDOBJ2IM           : GRIDOBJ2IM Create image from GRIDobj
- [ ] GRIDOBJ2MAT          : convert GRIDobj to matrix and coordinate vectors
- [ ] GRIDOBJ2PM           : combine several GRIDobj into a predictor matrix
- [ ] GRIDOBJ2POLYGON      : Conversion from drainage basin grid to polygon or polyline
- [ ] GRIDOBJ2RGB          : Convert GRIDobj to RGB image
- [x] ACV                  : Anisotropic coefficient of variation (ACV) 
- [x] AGGREGATE            : resampling a GRIDobj using aggregation
- [ ] ARCSLOPE             : mean gradient from a digital elevation model sensu ArcGIS
- [ ] ASPECT               : angle of exposition from a digital elevation model (GRIDobj)
- [x] CASTSHADOW           : cast shadow
- [ ] CELLAREA             : calculate cell areas of a GRIDobj in geographic coordinate system
- [ ] CLIP                 : clip a GRIDobj with a polygon or another GRIDobj
- [x] CONTOUR              : contour plot of an instance of GRIDobj
- [ ] COORD2IND            : convert x and y coordinates to linear index
- [ ] COORD2SUB            : convert x and y coordinates to subscripts into a GRIDobj
- [ ] CREATEMASK           : create a binary mask using polygon mapping
- [ ] CREATERECTMASK       : create a binary mask using rectangle mapping
- [ ] CROP                 : crop an instance of GRIDobj with axis-aligned minimum bounding box
- [ ] CURVATURE            : 8-connected neighborhood curvature of a digital elevation model 
- [ ] DEMAREA              : Calculate the corrected surface area of a DEM
- [ ] DEMPROFILE           : get profile along path
- [ ] DIFFUSION            : Solve the diffusion equation
- [ ] DILATE               : morphological dilation
- [ ] DIST2CURVE           : labels pixels in a GRIDobj by their directed distance to a curved line
- [ ] DIST2LINE            : labels pixels in a GRIDobj by their distance to a straight line
- [ ] DISTANCE             : distance transform
- [ ] ELEVATEMINIMA        : elevate regional minima in a DEM to their lowest neighbor
- [ ] ERODE                : morphological erosion
- [ ] EVANSSLOPE           : Calculate surface slope using Evans method
- [ ] EXCESSTOPOGRAPHY     : reconstruct surface with threshold-slope surface
- [ ] FILLSINKS            : fill/remove pits, sinks or topographic depressions
- [ ] FILTER               : 2D-filtering of DEMs with different kernels 
- [ ] FIND                 : Find indices of nonzero elements in GRIDobj
- [ ] FINDCOORD            : Find coordinates of nonzero elements in GRIDobj
- [x] GETCOORDINATES       : get coordinate vectors of an instance of GRIDobj
- [ ] GETEXTENT            : return extent of a GRIDobj
- [ ] GETOUTLINE           : get or plot extent of GRIDobj
- [ ] GRADIENT8            : 8-connected neighborhood gradient of a digital elevation model
- [x] GRIDDEDCONTOUR       : plot contours on grid
- [ ] HEXGRID              : creates an array of haxagonal points
- [ ] HILLSHADE            : create hillshading from a digital elevation model (GRIDobj)
- [ ] HISTOGRAM            : Plot frequency distribution of values in GRIDobj
- [ ] HYPSCURVE            : plot hypsometric curve of a digital elevation model
- [ ] IDENTIFYFLATS        : identify flat terrain in a digital elevation model
- [ ] IMAGESC              : Scale data in GRIDobj and display as image object
- [ ] IMAGESCHS            : plot hillshade image with overlay
- [ ] IND2COORD            : convert linear index to x and y coordinates
- [ ] INFO                 : detailed information on GRIDobj instance
- [ ] INPAINTNANS          : Interpolate or fill missing values in a grid (GRIDobj)
- [x] INTERP               : Interpolate to query locations
- [ ] INTERP2GRIDOBJ       : Interpolate scattered data to GRIDobj
- [ ] INTERPWITHBARRIERS   : Laplace interplation with barriers
- [ ] ISNAN                : returns array elements that are NaNs as logical grid
- [ ] KSDENSITY            : kernel density estimator for GRIDobj
- [ ] LARGESTINSCRIBEDGRID : Find and crop the largest grid with no nans
- [ ] LINE2GRIDOBJ         : convert line to a grid
- [ ] LOCALTOPOGRAPHY      : Local topography
- [x] MEASURE              : take interactive measurements along a polyline
Additional work to be done is to replace impoly with drawpolygon.
- [x] MINMAXNORM           : min-max normalization with optional percent clipping
- [ ] MPOWER               : overloaded power for GRIDobj
- [ ] MRDIVIDE             : overloaded right division for GRIDobj
- [ ] MTIMES               : overloaded multiplication for GRIDobj
- [ ] PAD                  : add or remove a border of pixels around a GRIDobj
- [ ] POLYGON2GRIDOBJ      : convert polygon to a grid
- [ ] POSTPROCFLATS        : postprocess flat terrain for visualization purpose
- [x] PRCCLIP              : percentile clipping
- [x] PROJECT              : transforms a GRIDobj between projected coordinate systems
- [ ] PROMINENCE           : Calculate the prominence of mountain peaks
- [x] RAND                 : Compute a GRIDobj with random numbers
- [ ] RANDOMSAMPLE         : Uniform random sampling of a GRIDobj
- [ ] RECLABEL             : labels GRIDobj by rectangular fields
- [ ] RECLASSIFY           : generate univariate class intervals for an instance of GRIDobj
- [x] REPROJECT2UTM        : Reproject DEM with WGS84 coordinate system to UTM-WGS84 
- [ ] RESAMPLE             : change spatial resolution of a GRIDobj
- [ ] ROUGHNESS            : terrain ruggedness, position and roughness indices of DEMs
- [ ] SHUFFLELABEL         : shufflelabel randomly relabels a label matrix
- [ ] SNAP2STREAM          : snap gauges or pour points to stream raster
- [ ] SUB2COORD            : convert subscripts to x and y coordinates
- [ ] SURF                 : surface plot for GRIDobj
- [ ] TANAKACONTOUR        : Relief depiction using Tanaka contours
- [ ] TOPOSHIELDING        : topographic shielding from cosmic rays
- [x] VALIDATEALIGNMENT    : validates whether instances of GRIDobj are spatially aligned
- [x] ZSCORE               : standardized z-scores for GRIDobj

## FLOWobj


## STREAMobj

