# Unblender Standalone Application repo

Design based on https://github.com/derryleng/Shiny_Desktop_App

This version of Unblender is a GUI Shiny app bundled with R 4.4.3 and Chrome portable.

When running from source on Windows:
- Clone repo to your system
- Download [Chrome Portable](https://portableapps.com/apps/internet/google_chrome_portable]) and install in `./chrome` folder. 
- Download [R 4.4.3](https://cloud.r-project.org/bin/windows/base/old/4.4.3/R-4.4.3-win.exe) and install in `./R` folder.
- Place the following required data files into `./app/datafiles`:
    * https://ironman.tenwisedev.nl/public/so_downs_5000.rds
    * https://ironman.tenwisedev.nl/public/example_gse76225.csv
- Run `run.bat`. (Note: on first run, all necessary R libraries are installed. This may take quite some time).
