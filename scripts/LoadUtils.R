local({
  if ("wta:utils" %in% search()) {
    warning("Environment `wta:utils` already attached and it will be refreshed",
            call. = FALSE)
    utils_env <- .WtaUtilsEnv
    detach("wta:utils")
  } else {
    utils_env <- new.env()
    assign(".WtaUtilsEnv", utils_env, envir = utils_env)
  }

  # LoadUtils.R was generally used in "analysis" folder which is a sibling
  # of "scripts", so here we use "../scripts/" to load all utilities.
  source("../scripts/UtilityFunctions.R", local = utils_env)
  source("../scripts/EnrichmentUtilities.R", local = utils_env)
  source("../scripts/TrajectoryUtilities.R", local = utils_env)
  source("../scripts/AbundanceUtilities.R", local = utils_env)

  args <- list(what = utils_env, name = "wta:utils")
  do.call(base::attach, args)
})
