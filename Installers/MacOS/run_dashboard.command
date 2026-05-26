#!/bin/zsh
cd "$(dirname "$0")/dashboard"
~/.pixi/bin/pixi run Rscript -e 'shiny::runApp(appDir="Unblender-app/app",launch.browser=TRUE)'