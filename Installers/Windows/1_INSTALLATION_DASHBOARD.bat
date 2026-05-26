:: Author: Gijs Bakker, James Gray
:: Version 0.1
:: email


:: ======= Check which instance =========
if "%1" == "" goto install_pixi
goto %1


:: ======= Install pixi =========
:install_pixi

:: Check if pixi is installed
pixi --version >nul 2>&1
if errorlevel 1 (
    echo "Install pixi"
    powershell -ExecutionPolicy ByPass -Command "irm -useb https://pixi.sh/install.ps1 | iex" || goto :error
    echo "Installed pixi"
    start "" %0 pixi_installed || goto :error
    goto end
) else (
    echo "Pixi is already installed!"
    goto pixi_installed
)


:pixi_installed
:: ======= Install R through pixi =========
echo "Install R"
pixi init dashboard || goto :error
cd dashboard || goto :error
pixi add "r-base=4.4.3" || goto :error

echo "✅ R version 4.4.3 installed succesfully!"

:: ======= Install libgit2 ======
echo "\tInstalling libgit2"
pixi add libgit2 || goto :error

:: ======= Getting App repository =========
echo "Getting latest dashboard files from repo"
curl -L "https://github.com/jrdalenberg/UnBlender/archive/refs/heads/main.zip" -o main.zip || goto :error

:: Extracting files
echo "Extracting dashboard!"
tar -xf main.zip || goto :error
del main.zip || goto :error

:: Removing reminants
move UnBlender-main/Unblender-app . || goto :error
rmdir /s /q UnBlender-main || goto :error

:: ======= Installing CRAN packages =======
echo "Installing required CRAN packages"
pixi add "r-seurat" "r-devtools" "r-shiny" "r-tidyverse" "r-shinydashboard" "r-shinyjs" "r-shinyjqui" "r-shinyTree" "r-dt" "r-writexl" "r-data.table" "r-nnls" "r-mcmcpack" "r-e1071" "r-locfdr" "r-quadprog" "r-statmod" "r-corpcor" "r-doparallel" "r-ggally" "r-markdown" "r-xfun" "r-promises"


:: ======= Installing BioConductor packages =======
echo "Installing required BioConductor packages"
pixi run Rscript ..\installR_scripts\install_BIOC_packages.R || goto :error

:: ======= Installing Github packages =======
echo "Installing required Github packages"
pixi run Rscript ..\installR_scripts\install_REPO_packages.R || goto :error

:: ======= Installing runtime documents =====
curl -L "https://ironman.tenwisedev.nl/public/so_downs_5000.rds" -o Unblender-app/app/datafiles/so_downs_5000.rds || goto :error
curl -L "https://ironman.tenwisedev.nl/public/example_gse76225.csv" -o Unblender-app/app/datafiles/example_gse76225.csv || goto :error

echo "✅ DASHBOARD INSTALLED SUCCESFULLY!"
goto :end

:end
exit /b

:error
pause
exit /b