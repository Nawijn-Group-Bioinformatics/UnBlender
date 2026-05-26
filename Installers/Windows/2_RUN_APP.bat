cd dashboard || goto :error
pixi run Rscript -e "shiny::runApp(appDir='Unblender-app/app',launch.browser=TRUE)" || goto :error
exit /b

:error
pause
exit /b