set ROOT_DIR=%cd%
call setenv.bat

:: we use vanilla so that we do not rely on any config file for R 
%ROOT_DIR%\lib\R\R-2.13.1\bin\i386\Rscript --vanilla segmentation.R parameters.R

pause
