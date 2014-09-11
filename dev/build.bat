set ROOT_DIR=%cd%\..
call %ROOT_DIR%\setenv.bat

:: This variables is used in EBImage\src\Makevars.win
set DEV_DIR=%cd%

:: Transform Windows back-slashes into Cygwin forward-slashes
set DEV_DIR=%DEV_DIR:\=^/%


:: we use vanilla so that we do not rely on any config file for R 
R --vanilla CMD INSTALL EBImage
pause
