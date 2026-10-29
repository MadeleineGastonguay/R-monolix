library(Rsmlx)
library(lixoftConnectors)
initializeLixoftConnectors(software="monolix",
                           path = "/Applications/MonolixSuite2024R1.app/Contents/Resources/monolixSuite/")

# Check that istallation works
checkInstallation <- function(){
  demoPath <- getDemoPath()
  loadProject(file.path(demoPath, "1.creating_and_using_models", "1.1.libraries_of_models", "theophylline_project.mlxtran"))
  runScenario()
  getEstimatedPopulationParameters()  
}

checkInstallation()

## Code for installation:
# install.packages("RJSONIO")
# install.packages("/Applications/MonolixSuite2024R1.app/Contents/Resources/monolixSuite/connectors/lixoftConnectors.tar.gz",
# repos = NULL, type="source", INSTALL_opts ="--no-multiarch")



