#!/bin/bash
cd dashboard
~/.pixi/bin/pixi run Rscript -e 'shiny::runApp(appDir="UnBlender-app/app")'
