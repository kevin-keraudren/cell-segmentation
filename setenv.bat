if !%ROOT_DIR%==! set ROOT_DIR=%cd%

set PATH=%ROOT_DIR%\dev\Rtools\bin;%ROOT_DIR%\dev\Rtools\MinGW\bin;%ROOT_DIR%\dev\Rtools\MinGW64\bin;%ROOT_DIR%\lib\gtk\bin;%ROOT_DIR%\lib\ImageMagick-6.6.9-10;%ROOT_DIR%\lib\R\R-2.13.1\bin;%ROOT_DIR%\lib\CED;%ROOT_DIR%\lib\EED

:: For local install of R libraries
set R_LIBS=%ROOT_DIR%\lib\R\R_HOME



