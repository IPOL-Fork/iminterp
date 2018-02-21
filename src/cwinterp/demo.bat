@echo off
rem Demo script for using cwinterp
rem Pascal Getreuer 2010

rem Image credits: the test image frog-hr.bmp is by
rem   John D. Willson, USGS Amphibian Research and Monitoring Initiative
rem   http://armi.usgs.gov/gallery/detail.php?search=Genus&subsearch=Bufo&id=323

echo.
echo +============================================================================+
echo + First, we coarsen high resolution image "frog-hr.bmp" by factor 4 to make  +
echo + low resolution image "frog-lr.bmp".                                        +
echo +============================================================================+

.\imcoarsen -g topleft -x4 -p0.35 frog-hr.bmp frog-lr.bmp

echo.
echo.
echo +============================================================================+
echo + Now we run the interpolation on the coarsened image "frog-lr.bmp"...       +
echo +============================================================================+

.\cwinterp -g topleft -x4 -p0.35 frog-lr.bmp interp.bmp

echo.
echo.
echo +============================================================================+
echo + The difference between the original and the interpolation is               +
echo +============================================================================+

.\imdiff frog-hr.bmp interp.bmp

echo.

pause
