#' Perform necessary tasks when the DiNAMICDuo package is loaded
#'
#'
miniconda_installation <- NULL
miniconda_permission <- NULL
numpy_import <- NULL

.onAttach <- function(libname, pkgname)
	{
	miniconda_installation <- utils::askYesNo("Is miniconda installed?")
	
	if (isFALSE(miniconda_installation))
		{
		miniconda_permission <- utils::askYesNo("Install miniconda?  Downloads 50MB and takes time.")

  		if (isTRUE(miniconda_permission))
			{
    			reticulate::install_miniconda()
    			} else{
    				packageStartupMessage("You should install miniconda before using this package")
				}

		numpy_import <- reticulate::import("numpy", delay_load = TRUE)
		}
	}


