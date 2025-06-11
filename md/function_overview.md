# Function overview
 
## GRIDobj
| Name    | Description |
| -------- | ------- |
| GRIDOBJ | Create instance of a GRIDobj | 
| GRIDOBJ2ASCII | write/export GRIDobj to ESRI ArcGIS ASCII file | 
| GRIDOBJ2GEOTABLE | Convert categorical GRIDobj to geotable (polygon)  | 
| GRIDOBJ2GEOTIFF | Exports an instance of GRIDobj to a geotiff file | 
| GRIDOBJ2IM | Create image from GRIDobj | 
| GRIDOBJ2MAT | Convert GRIDobj to matrix and coordinate vectors | 
| GRIDOBJ2PM | Combine several GRIDobjs into a predictor matrix | 
| GRIDOBJ2POLYGON | Conversion from drainage basin grid to polygon or polyline | 
| GRIDOBJ2RGB | Convert GRIDobj to RGB image | 
| ACV | Anisotropic coefficient of variation (ACV)  | 
| AGGREGATE | Resampling a GRIDobj using aggregation/binning | 
| ARCSLOPE | Slope from a digital elevation model sensu ArcGIS | 
| ASPECT | Direction of the steepest slope | 
| CASTSHADOW | Cast shadow calculated from digital terrain | 
| CELLAREA | Calculate cell areas of a GRIDobj in geographic coordinate system | 
| CLIP | clip a GRIDobj with a polygon or another GRIDobj | 
| CONTOUR | Contour plot of an instance of GRIDobj | 
| COORD2IND | convert x and y coordinates to linear index | 
| COORD2SUB | convert x and y coordinates to subscripts into a GRIDobj | 
| CREATEMASK | Create a binary mask using polygon mapping | 
| CROP | Crop an instance of GRIDobj with axis-aligned minimum bounding box | 
| CURVATURE | Curvature of a digital elevation model  | 
| DEMAREA | Calculate the corrected surface area of a DEM | 
| DEMPROFILE | Get profile along path | 
| DIFFUSION | Solve the diffusion equation | 
| DILATE | Morphological dilation | 
| DIST2CURVE | Labels pixels in a GRIDobj by their directed distance to a curved line | 
| DIST2LINE | labels pixels in a GRIDobj by their distance to a straight line | 
| DISTANCE | distance transform | 
| ELEVATEMINIMA | Elevate regional minima in a DEM to their lowest neighbor | 
| ERODE | Morphological erosion | 
| EVANSSLOPE | Calculate surface slope using Evans method | 
| EXCESSTOPOGRAPHY | reconstruct surface with threshold-slope surface | 
| FILLSINKS | Fill/remove pits, sinks or topographic depressions | 
| FILTER | 2D-filtering of DEMs with different kernels  | 
| FIND | Find indices of nonzero elements in GRIDobj | 
| FINDCOORD | Find coordinates of nonzero elements in GRIDobj | 
| GETCOORDINATES | get coordinate vectors of an instance of GRIDobj | 
| GETEXTENT | return extent of a GRIDobj | 
| GETOUTLINE | Get outline of GRIDobj | 
| GRADIENT8 | 8-connected neighborhood gradient of a digital elevation model | 
| GRIDDEDCONTOUR | plot contours on grid | 
| HEXGRID | creates an array of haxagonal points | 
| HILLSHADE | Calculate hillshading from a digital elevation model | 
| HISTOGRAM | Plot frequency distribution of values in GRIDobj | 
| HYDROGRAM | Generate an hydrogram | 
| HYPSCURVE | plot hypsometric curve of a digital elevation model | 
| IDENTIFYFLATS | identify flat terrain in a digital elevation model | 
| IMAGESC | Scale data in GRIDobj and display as image object | 
| IMAGESCHS | Plot hillshade image with overlay | 
| IND2COORD | convert linear index to x and y coordinates | 
| INFO | Detailed information on GRIDobj instance | 
| INPAINTNANS | Interpolate or fill missing values in a grid (GRIDobj) | 
| INTERP | Interpolate to query locations | 
| INTERP2GRIDOBJ | Interpolate scattered data to GRIDobj | 
| INTERPWITHBARRIERS | Laplace interplation with barriers | 
| ISUNDERLYINGINTEGER | Create instance of a GRIDobj | 
| ISUNDERLYINGNUMERIC | Create instance of a GRIDobj | 
| ISUNDERLYINGTYPE | Create instance of a GRIDobj | 
| ISNAN | Returns array elements that are NaNs as logical grid | 
| KSDENSITY | Kernel density estimator for GRIDobj | 
| LARGESTINSCRIBEDGRID | Find and crop the largest grid with no nans | 
| LINE2GRIDOBJ | Convert line to a grid | 
| LOCALTOPOGRAPHY | Local topography | 
| MAX | Maximum value in GRIDobj | 
| MEASURE | Take interactive measurements along a polyline | 
| MIN | Minimum value in GRIDobj | 
| MINMAXNORM | Min-max normalization with optional percent clipping | 
| MOSAIC | Merge multiple GRIDobjs into a larger GRIDobj | 
| MPOWER | overloaded power for GRIDobj | 
| MRDIVIDE | overloaded right division for GRIDobj | 
| MTIMES | overloaded multiplication for GRIDobj | 
| MUSTBEUNDERLYINGINTEGER | Create instance of a GRIDobj | 
| MUSTBEUNDERLYINGNUMERIC | Create instance of a GRIDobj | 
| MUSTBEUNDERLYINGTYPE | Create instance of a GRIDobj | 
| NCOLS | Return the number of columns in a DEM | 
| NROWS | Return the number of rows in a DEM | 
| PAD | Add or remove a border of pixels around a GRIDobj | 
| POLYGON2GRIDOBJ | Convert polygon to a grid | 
| POSTPROCFLATS | Postprocess flat terrain for visualization purpose | 
| PRCCLIP | Percentile clipping | 
| PROJECT | Warps a GRIDobj to a different coordinate system | 
| PROMINENCE | Calculate the prominence of mountain peaks | 
| RAND | Compute a GRIDobj with random, spatially uncorrelated numbers | 
| RANDOMSAMPLE | Uniform random sampling of a GRIDobj | 
| RECLABEL | Labels GRIDobj by rectangular fields | 
| RECLASSIFY | Generate univariate class intervals for an instance of GRIDobj | 
| REPROJECT2UTM | Reproject DEM with WGS84 coordinate system to UTM-WGS84  | 
| RESAMPLE | Change spatial resolution of a GRIDobj | 
| ROUGHNESS | Terrain ruggedness, position and roughness indices of DEMs | 
| SHUFFLELABEL | Shufflelabel randomly relabels a label matrix | 
| SNAP2STREAM | snap gauges or pour points to stream raster | 
| SUB2COORD | Convert subscripts to x and y coordinates | 
| SURF | Surface plot for GRIDobj | 
| TANAKACONTOUR | Relief depiction using Tanaka contours | 
| TOPOSHIELDING | topographic shielding from cosmic rays | 
| UNDERLYINGTYPE | Create instance of a GRIDobj | 
| VALIDATEALIGNMENT | Checks validity that two GRIDobj are spatially aligned | 
| ZSCORE | Standardized z-scores for GRIDobj | 
 
## FLOWobj
| Name    | Description |
| -------- | ------- |
| GRIDOBJ | Create instance of a GRIDobj | 
| GRIDOBJ2ASCII | write/export GRIDobj to ESRI ArcGIS ASCII file | 
| GRIDOBJ2GEOTABLE | Convert categorical GRIDobj to geotable (polygon)  | 
| GRIDOBJ2GEOTIFF | Exports an instance of GRIDobj to a geotiff file | 
| GRIDOBJ2IM | Create image from GRIDobj | 
| GRIDOBJ2MAT | Convert GRIDobj to matrix and coordinate vectors | 
| GRIDOBJ2PM | Combine several GRIDobjs into a predictor matrix | 
| GRIDOBJ2POLYGON | Conversion from drainage basin grid to polygon or polyline | 
| GRIDOBJ2RGB | Convert GRIDobj to RGB image | 
| ACV | Anisotropic coefficient of variation (ACV)  | 
| AGGREGATE | Resampling a GRIDobj using aggregation/binning | 
| ARCSLOPE | Slope from a digital elevation model sensu ArcGIS | 
| ASPECT | Direction of the steepest slope | 
| CASTSHADOW | Cast shadow calculated from digital terrain | 
| CELLAREA | Calculate cell areas of a GRIDobj in geographic coordinate system | 
| CLIP | clip a GRIDobj with a polygon or another GRIDobj | 
| CONTOUR | Contour plot of an instance of GRIDobj | 
| COORD2IND | convert x and y coordinates to linear index | 
| COORD2SUB | convert x and y coordinates to subscripts into a GRIDobj | 
| CREATEMASK | Create a binary mask using polygon mapping | 
| CROP | Crop an instance of GRIDobj with axis-aligned minimum bounding box | 
| CURVATURE | Curvature of a digital elevation model  | 
| DEMAREA | Calculate the corrected surface area of a DEM | 
| DEMPROFILE | Get profile along path | 
| DIFFUSION | Solve the diffusion equation | 
| DILATE | Morphological dilation | 
| DIST2CURVE | Labels pixels in a GRIDobj by their directed distance to a curved line | 
| DIST2LINE | labels pixels in a GRIDobj by their distance to a straight line | 
| DISTANCE | distance transform | 
| ELEVATEMINIMA | Elevate regional minima in a DEM to their lowest neighbor | 
| ERODE | Morphological erosion | 
| EVANSSLOPE | Calculate surface slope using Evans method | 
| EXCESSTOPOGRAPHY | reconstruct surface with threshold-slope surface | 
| FILLSINKS | Fill/remove pits, sinks or topographic depressions | 
| FILTER | 2D-filtering of DEMs with different kernels  | 
| FIND | Find indices of nonzero elements in GRIDobj | 
| FINDCOORD | Find coordinates of nonzero elements in GRIDobj | 
| GETCOORDINATES | get coordinate vectors of an instance of GRIDobj | 
| GETEXTENT | return extent of a GRIDobj | 
| GETOUTLINE | Get outline of GRIDobj | 
| GRADIENT8 | 8-connected neighborhood gradient of a digital elevation model | 
| GRIDDEDCONTOUR | plot contours on grid | 
| HEXGRID | creates an array of haxagonal points | 
| HILLSHADE | Calculate hillshading from a digital elevation model | 
| HISTOGRAM | Plot frequency distribution of values in GRIDobj | 
| HYDROGRAM | Generate an hydrogram | 
| HYPSCURVE | plot hypsometric curve of a digital elevation model | 
| IDENTIFYFLATS | identify flat terrain in a digital elevation model | 
| IMAGESC | Scale data in GRIDobj and display as image object | 
| IMAGESCHS | Plot hillshade image with overlay | 
| IND2COORD | convert linear index to x and y coordinates | 
| INFO | Detailed information on GRIDobj instance | 
| INPAINTNANS | Interpolate or fill missing values in a grid (GRIDobj) | 
| INTERP | Interpolate to query locations | 
| INTERP2GRIDOBJ | Interpolate scattered data to GRIDobj | 
| INTERPWITHBARRIERS | Laplace interplation with barriers | 
| ISUNDERLYINGINTEGER | Create instance of a GRIDobj | 
| ISUNDERLYINGNUMERIC | Create instance of a GRIDobj | 
| ISUNDERLYINGTYPE | Create instance of a GRIDobj | 
| ISNAN | Returns array elements that are NaNs as logical grid | 
| KSDENSITY | Kernel density estimator for GRIDobj | 
| LARGESTINSCRIBEDGRID | Find and crop the largest grid with no nans | 
| LINE2GRIDOBJ | Convert line to a grid | 
| LOCALTOPOGRAPHY | Local topography | 
| MAX | Maximum value in GRIDobj | 
| MEASURE | Take interactive measurements along a polyline | 
| MIN | Minimum value in GRIDobj | 
| MINMAXNORM | Min-max normalization with optional percent clipping | 
| MOSAIC | Merge multiple GRIDobjs into a larger GRIDobj | 
| MPOWER | overloaded power for GRIDobj | 
| MRDIVIDE | overloaded right division for GRIDobj | 
| MTIMES | overloaded multiplication for GRIDobj | 
| MUSTBEUNDERLYINGINTEGER | Create instance of a GRIDobj | 
| MUSTBEUNDERLYINGNUMERIC | Create instance of a GRIDobj | 
| MUSTBEUNDERLYINGTYPE | Create instance of a GRIDobj | 
| NCOLS | Return the number of columns in a DEM | 
| NROWS | Return the number of rows in a DEM | 
| PAD | Add or remove a border of pixels around a GRIDobj | 
| POLYGON2GRIDOBJ | Convert polygon to a grid | 
| POSTPROCFLATS | Postprocess flat terrain for visualization purpose | 
| PRCCLIP | Percentile clipping | 
| PROJECT | Warps a GRIDobj to a different coordinate system | 
| PROMINENCE | Calculate the prominence of mountain peaks | 
| RAND | Compute a GRIDobj with random, spatially uncorrelated numbers | 
| RANDOMSAMPLE | Uniform random sampling of a GRIDobj | 
| RECLABEL | Labels GRIDobj by rectangular fields | 
| RECLASSIFY | Generate univariate class intervals for an instance of GRIDobj | 
| REPROJECT2UTM | Reproject DEM with WGS84 coordinate system to UTM-WGS84  | 
| RESAMPLE | Change spatial resolution of a GRIDobj | 
| ROUGHNESS | Terrain ruggedness, position and roughness indices of DEMs | 
| SHUFFLELABEL | Shufflelabel randomly relabels a label matrix | 
| SNAP2STREAM | snap gauges or pour points to stream raster | 
| SUB2COORD | Convert subscripts to x and y coordinates | 
| SURF | Surface plot for GRIDobj | 
| TANAKACONTOUR | Relief depiction using Tanaka contours | 
| TOPOSHIELDING | topographic shielding from cosmic rays | 
| UNDERLYINGTYPE | Create instance of a GRIDobj | 
| VALIDATEALIGNMENT | Checks validity that two GRIDobj are spatially aligned | 
| ZSCORE | Standardized z-scores for GRIDobj | 
 
## STREAMobj
| Name    | Description |
| -------- | ------- |
| GRIDOBJ | Create instance of a GRIDobj | 
| GRIDOBJ2ASCII | write/export GRIDobj to ESRI ArcGIS ASCII file | 
| GRIDOBJ2GEOTABLE | Convert categorical GRIDobj to geotable (polygon)  | 
| GRIDOBJ2GEOTIFF | Exports an instance of GRIDobj to a geotiff file | 
| GRIDOBJ2IM | Create image from GRIDobj | 
| GRIDOBJ2MAT | Convert GRIDobj to matrix and coordinate vectors | 
| GRIDOBJ2PM | Combine several GRIDobjs into a predictor matrix | 
| GRIDOBJ2POLYGON | Conversion from drainage basin grid to polygon or polyline | 
| GRIDOBJ2RGB | Convert GRIDobj to RGB image | 
| ACV | Anisotropic coefficient of variation (ACV)  | 
| AGGREGATE | Resampling a GRIDobj using aggregation/binning | 
| ARCSLOPE | Slope from a digital elevation model sensu ArcGIS | 
| ASPECT | Direction of the steepest slope | 
| CASTSHADOW | Cast shadow calculated from digital terrain | 
| CELLAREA | Calculate cell areas of a GRIDobj in geographic coordinate system | 
| CLIP | clip a GRIDobj with a polygon or another GRIDobj | 
| CONTOUR | Contour plot of an instance of GRIDobj | 
| COORD2IND | convert x and y coordinates to linear index | 
| COORD2SUB | convert x and y coordinates to subscripts into a GRIDobj | 
| CREATEMASK | Create a binary mask using polygon mapping | 
| CROP | Crop an instance of GRIDobj with axis-aligned minimum bounding box | 
| CURVATURE | Curvature of a digital elevation model  | 
| DEMAREA | Calculate the corrected surface area of a DEM | 
| DEMPROFILE | Get profile along path | 
| DIFFUSION | Solve the diffusion equation | 
| DILATE | Morphological dilation | 
| DIST2CURVE | Labels pixels in a GRIDobj by their directed distance to a curved line | 
| DIST2LINE | labels pixels in a GRIDobj by their distance to a straight line | 
| DISTANCE | distance transform | 
| ELEVATEMINIMA | Elevate regional minima in a DEM to their lowest neighbor | 
| ERODE | Morphological erosion | 
| EVANSSLOPE | Calculate surface slope using Evans method | 
| EXCESSTOPOGRAPHY | reconstruct surface with threshold-slope surface | 
| FILLSINKS | Fill/remove pits, sinks or topographic depressions | 
| FILTER | 2D-filtering of DEMs with different kernels  | 
| FIND | Find indices of nonzero elements in GRIDobj | 
| FINDCOORD | Find coordinates of nonzero elements in GRIDobj | 
| GETCOORDINATES | get coordinate vectors of an instance of GRIDobj | 
| GETEXTENT | return extent of a GRIDobj | 
| GETOUTLINE | Get outline of GRIDobj | 
| GRADIENT8 | 8-connected neighborhood gradient of a digital elevation model | 
| GRIDDEDCONTOUR | plot contours on grid | 
| HEXGRID | creates an array of haxagonal points | 
| HILLSHADE | Calculate hillshading from a digital elevation model | 
| HISTOGRAM | Plot frequency distribution of values in GRIDobj | 
| HYDROGRAM | Generate an hydrogram | 
| HYPSCURVE | plot hypsometric curve of a digital elevation model | 
| IDENTIFYFLATS | identify flat terrain in a digital elevation model | 
| IMAGESC | Scale data in GRIDobj and display as image object | 
| IMAGESCHS | Plot hillshade image with overlay | 
| IND2COORD | convert linear index to x and y coordinates | 
| INFO | Detailed information on GRIDobj instance | 
| INPAINTNANS | Interpolate or fill missing values in a grid (GRIDobj) | 
| INTERP | Interpolate to query locations | 
| INTERP2GRIDOBJ | Interpolate scattered data to GRIDobj | 
| INTERPWITHBARRIERS | Laplace interplation with barriers | 
| ISUNDERLYINGINTEGER | Create instance of a GRIDobj | 
| ISUNDERLYINGNUMERIC | Create instance of a GRIDobj | 
| ISUNDERLYINGTYPE | Create instance of a GRIDobj | 
| ISNAN | Returns array elements that are NaNs as logical grid | 
| KSDENSITY | Kernel density estimator for GRIDobj | 
| LARGESTINSCRIBEDGRID | Find and crop the largest grid with no nans | 
| LINE2GRIDOBJ | Convert line to a grid | 
| LOCALTOPOGRAPHY | Local topography | 
| MAX | Maximum value in GRIDobj | 
| MEASURE | Take interactive measurements along a polyline | 
| MIN | Minimum value in GRIDobj | 
| MINMAXNORM | Min-max normalization with optional percent clipping | 
| MOSAIC | Merge multiple GRIDobjs into a larger GRIDobj | 
| MPOWER | overloaded power for GRIDobj | 
| MRDIVIDE | overloaded right division for GRIDobj | 
| MTIMES | overloaded multiplication for GRIDobj | 
| MUSTBEUNDERLYINGINTEGER | Create instance of a GRIDobj | 
| MUSTBEUNDERLYINGNUMERIC | Create instance of a GRIDobj | 
| MUSTBEUNDERLYINGTYPE | Create instance of a GRIDobj | 
| NCOLS | Return the number of columns in a DEM | 
| NROWS | Return the number of rows in a DEM | 
| PAD | Add or remove a border of pixels around a GRIDobj | 
| POLYGON2GRIDOBJ | Convert polygon to a grid | 
| POSTPROCFLATS | Postprocess flat terrain for visualization purpose | 
| PRCCLIP | Percentile clipping | 
| PROJECT | Warps a GRIDobj to a different coordinate system | 
| PROMINENCE | Calculate the prominence of mountain peaks | 
| RAND | Compute a GRIDobj with random, spatially uncorrelated numbers | 
| RANDOMSAMPLE | Uniform random sampling of a GRIDobj | 
| RECLABEL | Labels GRIDobj by rectangular fields | 
| RECLASSIFY | Generate univariate class intervals for an instance of GRIDobj | 
| REPROJECT2UTM | Reproject DEM with WGS84 coordinate system to UTM-WGS84  | 
| RESAMPLE | Change spatial resolution of a GRIDobj | 
| ROUGHNESS | Terrain ruggedness, position and roughness indices of DEMs | 
| SHUFFLELABEL | Shufflelabel randomly relabels a label matrix | 
| SNAP2STREAM | snap gauges or pour points to stream raster | 
| SUB2COORD | Convert subscripts to x and y coordinates | 
| SURF | Surface plot for GRIDobj | 
| TANAKACONTOUR | Relief depiction using Tanaka contours | 
| TOPOSHIELDING | topographic shielding from cosmic rays | 
| UNDERLYINGTYPE | Create instance of a GRIDobj | 
| VALIDATEALIGNMENT | Checks validity that two GRIDobj are spatially aligned | 
| ZSCORE | Standardized z-scores for GRIDobj | 
 
## SWATHobj
| Name    | Description |
| -------- | ------- |
| GRIDOBJ | Create instance of a GRIDobj | 
| GRIDOBJ2ASCII | write/export GRIDobj to ESRI ArcGIS ASCII file | 
| GRIDOBJ2GEOTABLE | Convert categorical GRIDobj to geotable (polygon)  | 
| GRIDOBJ2GEOTIFF | Exports an instance of GRIDobj to a geotiff file | 
| GRIDOBJ2IM | Create image from GRIDobj | 
| GRIDOBJ2MAT | Convert GRIDobj to matrix and coordinate vectors | 
| GRIDOBJ2PM | Combine several GRIDobjs into a predictor matrix | 
| GRIDOBJ2POLYGON | Conversion from drainage basin grid to polygon or polyline | 
| GRIDOBJ2RGB | Convert GRIDobj to RGB image | 
| ACV | Anisotropic coefficient of variation (ACV)  | 
| AGGREGATE | Resampling a GRIDobj using aggregation/binning | 
| ARCSLOPE | Slope from a digital elevation model sensu ArcGIS | 
| ASPECT | Direction of the steepest slope | 
| CASTSHADOW | Cast shadow calculated from digital terrain | 
| CELLAREA | Calculate cell areas of a GRIDobj in geographic coordinate system | 
| CLIP | clip a GRIDobj with a polygon or another GRIDobj | 
| CONTOUR | Contour plot of an instance of GRIDobj | 
| COORD2IND | convert x and y coordinates to linear index | 
| COORD2SUB | convert x and y coordinates to subscripts into a GRIDobj | 
| CREATEMASK | Create a binary mask using polygon mapping | 
| CROP | Crop an instance of GRIDobj with axis-aligned minimum bounding box | 
| CURVATURE | Curvature of a digital elevation model  | 
| DEMAREA | Calculate the corrected surface area of a DEM | 
| DEMPROFILE | Get profile along path | 
| DIFFUSION | Solve the diffusion equation | 
| DILATE | Morphological dilation | 
| DIST2CURVE | Labels pixels in a GRIDobj by their directed distance to a curved line | 
| DIST2LINE | labels pixels in a GRIDobj by their distance to a straight line | 
| DISTANCE | distance transform | 
| ELEVATEMINIMA | Elevate regional minima in a DEM to their lowest neighbor | 
| ERODE | Morphological erosion | 
| EVANSSLOPE | Calculate surface slope using Evans method | 
| EXCESSTOPOGRAPHY | reconstruct surface with threshold-slope surface | 
| FILLSINKS | Fill/remove pits, sinks or topographic depressions | 
| FILTER | 2D-filtering of DEMs with different kernels  | 
| FIND | Find indices of nonzero elements in GRIDobj | 
| FINDCOORD | Find coordinates of nonzero elements in GRIDobj | 
| GETCOORDINATES | get coordinate vectors of an instance of GRIDobj | 
| GETEXTENT | return extent of a GRIDobj | 
| GETOUTLINE | Get outline of GRIDobj | 
| GRADIENT8 | 8-connected neighborhood gradient of a digital elevation model | 
| GRIDDEDCONTOUR | plot contours on grid | 
| HEXGRID | creates an array of haxagonal points | 
| HILLSHADE | Calculate hillshading from a digital elevation model | 
| HISTOGRAM | Plot frequency distribution of values in GRIDobj | 
| HYDROGRAM | Generate an hydrogram | 
| HYPSCURVE | plot hypsometric curve of a digital elevation model | 
| IDENTIFYFLATS | identify flat terrain in a digital elevation model | 
| IMAGESC | Scale data in GRIDobj and display as image object | 
| IMAGESCHS | Plot hillshade image with overlay | 
| IND2COORD | convert linear index to x and y coordinates | 
| INFO | Detailed information on GRIDobj instance | 
| INPAINTNANS | Interpolate or fill missing values in a grid (GRIDobj) | 
| INTERP | Interpolate to query locations | 
| INTERP2GRIDOBJ | Interpolate scattered data to GRIDobj | 
| INTERPWITHBARRIERS | Laplace interplation with barriers | 
| ISUNDERLYINGINTEGER | Create instance of a GRIDobj | 
| ISUNDERLYINGNUMERIC | Create instance of a GRIDobj | 
| ISUNDERLYINGTYPE | Create instance of a GRIDobj | 
| ISNAN | Returns array elements that are NaNs as logical grid | 
| KSDENSITY | Kernel density estimator for GRIDobj | 
| LARGESTINSCRIBEDGRID | Find and crop the largest grid with no nans | 
| LINE2GRIDOBJ | Convert line to a grid | 
| LOCALTOPOGRAPHY | Local topography | 
| MAX | Maximum value in GRIDobj | 
| MEASURE | Take interactive measurements along a polyline | 
| MIN | Minimum value in GRIDobj | 
| MINMAXNORM | Min-max normalization with optional percent clipping | 
| MOSAIC | Merge multiple GRIDobjs into a larger GRIDobj | 
| MPOWER | overloaded power for GRIDobj | 
| MRDIVIDE | overloaded right division for GRIDobj | 
| MTIMES | overloaded multiplication for GRIDobj | 
| MUSTBEUNDERLYINGINTEGER | Create instance of a GRIDobj | 
| MUSTBEUNDERLYINGNUMERIC | Create instance of a GRIDobj | 
| MUSTBEUNDERLYINGTYPE | Create instance of a GRIDobj | 
| NCOLS | Return the number of columns in a DEM | 
| NROWS | Return the number of rows in a DEM | 
| PAD | Add or remove a border of pixels around a GRIDobj | 
| POLYGON2GRIDOBJ | Convert polygon to a grid | 
| POSTPROCFLATS | Postprocess flat terrain for visualization purpose | 
| PRCCLIP | Percentile clipping | 
| PROJECT | Warps a GRIDobj to a different coordinate system | 
| PROMINENCE | Calculate the prominence of mountain peaks | 
| RAND | Compute a GRIDobj with random, spatially uncorrelated numbers | 
| RANDOMSAMPLE | Uniform random sampling of a GRIDobj | 
| RECLABEL | Labels GRIDobj by rectangular fields | 
| RECLASSIFY | Generate univariate class intervals for an instance of GRIDobj | 
| REPROJECT2UTM | Reproject DEM with WGS84 coordinate system to UTM-WGS84  | 
| RESAMPLE | Change spatial resolution of a GRIDobj | 
| ROUGHNESS | Terrain ruggedness, position and roughness indices of DEMs | 
| SHUFFLELABEL | Shufflelabel randomly relabels a label matrix | 
| SNAP2STREAM | snap gauges or pour points to stream raster | 
| SUB2COORD | Convert subscripts to x and y coordinates | 
| SURF | Surface plot for GRIDobj | 
| TANAKACONTOUR | Relief depiction using Tanaka contours | 
| TOPOSHIELDING | topographic shielding from cosmic rays | 
| UNDERLYINGTYPE | Create instance of a GRIDobj | 
| VALIDATEALIGNMENT | Checks validity that two GRIDobj are spatially aligned | 
| ZSCORE | Standardized z-scores for GRIDobj | 
 
## PPS
| Name    | Description |
| -------- | ------- |
| GRIDOBJ | Create instance of a GRIDobj | 
| GRIDOBJ2ASCII | write/export GRIDobj to ESRI ArcGIS ASCII file | 
| GRIDOBJ2GEOTABLE | Convert categorical GRIDobj to geotable (polygon)  | 
| GRIDOBJ2GEOTIFF | Exports an instance of GRIDobj to a geotiff file | 
| GRIDOBJ2IM | Create image from GRIDobj | 
| GRIDOBJ2MAT | Convert GRIDobj to matrix and coordinate vectors | 
| GRIDOBJ2PM | Combine several GRIDobjs into a predictor matrix | 
| GRIDOBJ2POLYGON | Conversion from drainage basin grid to polygon or polyline | 
| GRIDOBJ2RGB | Convert GRIDobj to RGB image | 
| ACV | Anisotropic coefficient of variation (ACV)  | 
| AGGREGATE | Resampling a GRIDobj using aggregation/binning | 
| ARCSLOPE | Slope from a digital elevation model sensu ArcGIS | 
| ASPECT | Direction of the steepest slope | 
| CASTSHADOW | Cast shadow calculated from digital terrain | 
| CELLAREA | Calculate cell areas of a GRIDobj in geographic coordinate system | 
| CLIP | clip a GRIDobj with a polygon or another GRIDobj | 
| CONTOUR | Contour plot of an instance of GRIDobj | 
| COORD2IND | convert x and y coordinates to linear index | 
| COORD2SUB | convert x and y coordinates to subscripts into a GRIDobj | 
| CREATEMASK | Create a binary mask using polygon mapping | 
| CROP | Crop an instance of GRIDobj with axis-aligned minimum bounding box | 
| CURVATURE | Curvature of a digital elevation model  | 
| DEMAREA | Calculate the corrected surface area of a DEM | 
| DEMPROFILE | Get profile along path | 
| DIFFUSION | Solve the diffusion equation | 
| DILATE | Morphological dilation | 
| DIST2CURVE | Labels pixels in a GRIDobj by their directed distance to a curved line | 
| DIST2LINE | labels pixels in a GRIDobj by their distance to a straight line | 
| DISTANCE | distance transform | 
| ELEVATEMINIMA | Elevate regional minima in a DEM to their lowest neighbor | 
| ERODE | Morphological erosion | 
| EVANSSLOPE | Calculate surface slope using Evans method | 
| EXCESSTOPOGRAPHY | reconstruct surface with threshold-slope surface | 
| FILLSINKS | Fill/remove pits, sinks or topographic depressions | 
| FILTER | 2D-filtering of DEMs with different kernels  | 
| FIND | Find indices of nonzero elements in GRIDobj | 
| FINDCOORD | Find coordinates of nonzero elements in GRIDobj | 
| GETCOORDINATES | get coordinate vectors of an instance of GRIDobj | 
| GETEXTENT | return extent of a GRIDobj | 
| GETOUTLINE | Get outline of GRIDobj | 
| GRADIENT8 | 8-connected neighborhood gradient of a digital elevation model | 
| GRIDDEDCONTOUR | plot contours on grid | 
| HEXGRID | creates an array of haxagonal points | 
| HILLSHADE | Calculate hillshading from a digital elevation model | 
| HISTOGRAM | Plot frequency distribution of values in GRIDobj | 
| HYDROGRAM | Generate an hydrogram | 
| HYPSCURVE | plot hypsometric curve of a digital elevation model | 
| IDENTIFYFLATS | identify flat terrain in a digital elevation model | 
| IMAGESC | Scale data in GRIDobj and display as image object | 
| IMAGESCHS | Plot hillshade image with overlay | 
| IND2COORD | convert linear index to x and y coordinates | 
| INFO | Detailed information on GRIDobj instance | 
| INPAINTNANS | Interpolate or fill missing values in a grid (GRIDobj) | 
| INTERP | Interpolate to query locations | 
| INTERP2GRIDOBJ | Interpolate scattered data to GRIDobj | 
| INTERPWITHBARRIERS | Laplace interplation with barriers | 
| ISUNDERLYINGINTEGER | Create instance of a GRIDobj | 
| ISUNDERLYINGNUMERIC | Create instance of a GRIDobj | 
| ISUNDERLYINGTYPE | Create instance of a GRIDobj | 
| ISNAN | Returns array elements that are NaNs as logical grid | 
| KSDENSITY | Kernel density estimator for GRIDobj | 
| LARGESTINSCRIBEDGRID | Find and crop the largest grid with no nans | 
| LINE2GRIDOBJ | Convert line to a grid | 
| LOCALTOPOGRAPHY | Local topography | 
| MAX | Maximum value in GRIDobj | 
| MEASURE | Take interactive measurements along a polyline | 
| MIN | Minimum value in GRIDobj | 
| MINMAXNORM | Min-max normalization with optional percent clipping | 
| MOSAIC | Merge multiple GRIDobjs into a larger GRIDobj | 
| MPOWER | overloaded power for GRIDobj | 
| MRDIVIDE | overloaded right division for GRIDobj | 
| MTIMES | overloaded multiplication for GRIDobj | 
| MUSTBEUNDERLYINGINTEGER | Create instance of a GRIDobj | 
| MUSTBEUNDERLYINGNUMERIC | Create instance of a GRIDobj | 
| MUSTBEUNDERLYINGTYPE | Create instance of a GRIDobj | 
| NCOLS | Return the number of columns in a DEM | 
| NROWS | Return the number of rows in a DEM | 
| PAD | Add or remove a border of pixels around a GRIDobj | 
| POLYGON2GRIDOBJ | Convert polygon to a grid | 
| POSTPROCFLATS | Postprocess flat terrain for visualization purpose | 
| PRCCLIP | Percentile clipping | 
| PROJECT | Warps a GRIDobj to a different coordinate system | 
| PROMINENCE | Calculate the prominence of mountain peaks | 
| RAND | Compute a GRIDobj with random, spatially uncorrelated numbers | 
| RANDOMSAMPLE | Uniform random sampling of a GRIDobj | 
| RECLABEL | Labels GRIDobj by rectangular fields | 
| RECLASSIFY | Generate univariate class intervals for an instance of GRIDobj | 
| REPROJECT2UTM | Reproject DEM with WGS84 coordinate system to UTM-WGS84  | 
| RESAMPLE | Change spatial resolution of a GRIDobj | 
| ROUGHNESS | Terrain ruggedness, position and roughness indices of DEMs | 
| SHUFFLELABEL | Shufflelabel randomly relabels a label matrix | 
| SNAP2STREAM | snap gauges or pour points to stream raster | 
| SUB2COORD | Convert subscripts to x and y coordinates | 
| SURF | Surface plot for GRIDobj | 
| TANAKACONTOUR | Relief depiction using Tanaka contours | 
| TOPOSHIELDING | topographic shielding from cosmic rays | 
| UNDERLYINGTYPE | Create instance of a GRIDobj | 
| VALIDATEALIGNMENT | Checks validity that two GRIDobj are spatially aligned | 
| ZSCORE | Standardized z-scores for GRIDobj | 
 
## DIVIDEobj
| Name    | Description |
| -------- | ------- |
| GRIDOBJ | Create instance of a GRIDobj | 
| GRIDOBJ2ASCII | write/export GRIDobj to ESRI ArcGIS ASCII file | 
| GRIDOBJ2GEOTABLE | Convert categorical GRIDobj to geotable (polygon)  | 
| GRIDOBJ2GEOTIFF | Exports an instance of GRIDobj to a geotiff file | 
| GRIDOBJ2IM | Create image from GRIDobj | 
| GRIDOBJ2MAT | Convert GRIDobj to matrix and coordinate vectors | 
| GRIDOBJ2PM | Combine several GRIDobjs into a predictor matrix | 
| GRIDOBJ2POLYGON | Conversion from drainage basin grid to polygon or polyline | 
| GRIDOBJ2RGB | Convert GRIDobj to RGB image | 
| ACV | Anisotropic coefficient of variation (ACV)  | 
| AGGREGATE | Resampling a GRIDobj using aggregation/binning | 
| ARCSLOPE | Slope from a digital elevation model sensu ArcGIS | 
| ASPECT | Direction of the steepest slope | 
| CASTSHADOW | Cast shadow calculated from digital terrain | 
| CELLAREA | Calculate cell areas of a GRIDobj in geographic coordinate system | 
| CLIP | clip a GRIDobj with a polygon or another GRIDobj | 
| CONTOUR | Contour plot of an instance of GRIDobj | 
| COORD2IND | convert x and y coordinates to linear index | 
| COORD2SUB | convert x and y coordinates to subscripts into a GRIDobj | 
| CREATEMASK | Create a binary mask using polygon mapping | 
| CROP | Crop an instance of GRIDobj with axis-aligned minimum bounding box | 
| CURVATURE | Curvature of a digital elevation model  | 
| DEMAREA | Calculate the corrected surface area of a DEM | 
| DEMPROFILE | Get profile along path | 
| DIFFUSION | Solve the diffusion equation | 
| DILATE | Morphological dilation | 
| DIST2CURVE | Labels pixels in a GRIDobj by their directed distance to a curved line | 
| DIST2LINE | labels pixels in a GRIDobj by their distance to a straight line | 
| DISTANCE | distance transform | 
| ELEVATEMINIMA | Elevate regional minima in a DEM to their lowest neighbor | 
| ERODE | Morphological erosion | 
| EVANSSLOPE | Calculate surface slope using Evans method | 
| EXCESSTOPOGRAPHY | reconstruct surface with threshold-slope surface | 
| FILLSINKS | Fill/remove pits, sinks or topographic depressions | 
| FILTER | 2D-filtering of DEMs with different kernels  | 
| FIND | Find indices of nonzero elements in GRIDobj | 
| FINDCOORD | Find coordinates of nonzero elements in GRIDobj | 
| GETCOORDINATES | get coordinate vectors of an instance of GRIDobj | 
| GETEXTENT | return extent of a GRIDobj | 
| GETOUTLINE | Get outline of GRIDobj | 
| GRADIENT8 | 8-connected neighborhood gradient of a digital elevation model | 
| GRIDDEDCONTOUR | plot contours on grid | 
| HEXGRID | creates an array of haxagonal points | 
| HILLSHADE | Calculate hillshading from a digital elevation model | 
| HISTOGRAM | Plot frequency distribution of values in GRIDobj | 
| HYDROGRAM | Generate an hydrogram | 
| HYPSCURVE | plot hypsometric curve of a digital elevation model | 
| IDENTIFYFLATS | identify flat terrain in a digital elevation model | 
| IMAGESC | Scale data in GRIDobj and display as image object | 
| IMAGESCHS | Plot hillshade image with overlay | 
| IND2COORD | convert linear index to x and y coordinates | 
| INFO | Detailed information on GRIDobj instance | 
| INPAINTNANS | Interpolate or fill missing values in a grid (GRIDobj) | 
| INTERP | Interpolate to query locations | 
| INTERP2GRIDOBJ | Interpolate scattered data to GRIDobj | 
| INTERPWITHBARRIERS | Laplace interplation with barriers | 
| ISUNDERLYINGINTEGER | Create instance of a GRIDobj | 
| ISUNDERLYINGNUMERIC | Create instance of a GRIDobj | 
| ISUNDERLYINGTYPE | Create instance of a GRIDobj | 
| ISNAN | Returns array elements that are NaNs as logical grid | 
| KSDENSITY | Kernel density estimator for GRIDobj | 
| LARGESTINSCRIBEDGRID | Find and crop the largest grid with no nans | 
| LINE2GRIDOBJ | Convert line to a grid | 
| LOCALTOPOGRAPHY | Local topography | 
| MAX | Maximum value in GRIDobj | 
| MEASURE | Take interactive measurements along a polyline | 
| MIN | Minimum value in GRIDobj | 
| MINMAXNORM | Min-max normalization with optional percent clipping | 
| MOSAIC | Merge multiple GRIDobjs into a larger GRIDobj | 
| MPOWER | overloaded power for GRIDobj | 
| MRDIVIDE | overloaded right division for GRIDobj | 
| MTIMES | overloaded multiplication for GRIDobj | 
| MUSTBEUNDERLYINGINTEGER | Create instance of a GRIDobj | 
| MUSTBEUNDERLYINGNUMERIC | Create instance of a GRIDobj | 
| MUSTBEUNDERLYINGTYPE | Create instance of a GRIDobj | 
| NCOLS | Return the number of columns in a DEM | 
| NROWS | Return the number of rows in a DEM | 
| PAD | Add or remove a border of pixels around a GRIDobj | 
| POLYGON2GRIDOBJ | Convert polygon to a grid | 
| POSTPROCFLATS | Postprocess flat terrain for visualization purpose | 
| PRCCLIP | Percentile clipping | 
| PROJECT | Warps a GRIDobj to a different coordinate system | 
| PROMINENCE | Calculate the prominence of mountain peaks | 
| RAND | Compute a GRIDobj with random, spatially uncorrelated numbers | 
| RANDOMSAMPLE | Uniform random sampling of a GRIDobj | 
| RECLABEL | Labels GRIDobj by rectangular fields | 
| RECLASSIFY | Generate univariate class intervals for an instance of GRIDobj | 
| REPROJECT2UTM | Reproject DEM with WGS84 coordinate system to UTM-WGS84  | 
| RESAMPLE | Change spatial resolution of a GRIDobj | 
| ROUGHNESS | Terrain ruggedness, position and roughness indices of DEMs | 
| SHUFFLELABEL | Shufflelabel randomly relabels a label matrix | 
| SNAP2STREAM | snap gauges or pour points to stream raster | 
| SUB2COORD | Convert subscripts to x and y coordinates | 
| SURF | Surface plot for GRIDobj | 
| TANAKACONTOUR | Relief depiction using Tanaka contours | 
| TOPOSHIELDING | topographic shielding from cosmic rays | 
| UNDERLYINGTYPE | Create instance of a GRIDobj | 
| VALIDATEALIGNMENT | Checks validity that two GRIDobj are spatially aligned | 
| ZSCORE | Standardized z-scores for GRIDobj | 
 
