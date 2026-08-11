library(PioneerMetastaticAE)
# ------------------------------------------------------------------------------
# 0. Environment
# ------------------------------------------------------------------------------
source("../credentials.r")

# ------------------------------------------------------------------------------
# 1. Site-specific parameters  (edit these for every site)
# ------------------------------------------------------------------------------

JDBC <- "/home/a_kiselev/Jdbc_SAT3"

DatabaseConnector::downloadJdbcDrivers(
  "redshift",
  pathToDriver = JDBC,
  method = "auto",
)

connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms        = DBMS,
  user        = USER,
  password    = PASSWORD,
  server      = SERVER,
  port        = DB_PORT,
  pathToDriver = JDBC
)

## -- Schema names -------------------------------------------------------------
cdmDatabaseSchema        <- "marketscan_ccaemdcr_aug2025"# read-only CDM
vocabularyDatabaseSchema <- cdmDatabaseSchema     # usually same as CDM
cohortDatabaseSchema     <- "marketscan_ccaemdcr_aug2025_results"        # write-enabled; needs read/write/delete
cohortTable              <- "cohort_SAT3"

databaseId          <- "MarketScan"           # short identifier, no spaces
databaseName        <- "MarketScan aug2025"
databaseDescription <- "MarketScan aug2025"


## -- Temp table emulation (Oracle / some Spark configs) -----------------------
options(sqlRenderTempEmulationSchema = NULL)
# Maximum number of cores to be used:
maxCores <- parallel::detectCores()

# The folder where the study intermediate and result files will be written:
outputFolder <- file.path("/home/a_kiselev/output", paste0(databaseId, "_SAT3_PioneerMetastaticAE"))

# Specify another folder where the Incidence rate results will be written:
outputFolderIR<- file.path("/home/a_kiselev/output", paste0(databaseId, "_SAT3_PioneerMetastaticAE"))
# Optional: specify where the temporary files (used by the Andromeda package) will be created:
options(andromedaTempFolder = file.path(outputFolder, "andromedaTemp"))



#Age group breaks
breaks_age <- c(18, 55, 75, 120)

# Charlson Comorbidity Index category breaks (interior cut points).
# The default c(0, 2, 4) yields the strata 0, 1-2, 3-4 and 5+.
breaks_charlson <- c(0, 2, 4)

PioneerMetastaticAE::execute(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTable = cohortTable,
  verifyDependencies = TRUE,
  outputFolder = outputFolder,
  databaseId = databaseId,
  databaseName = databaseName,
  databaseDescription = databaseDescription
)


# Execute the incidence rate analysis
execute_IR  (connectionDetails=connectionDetails,
                       cdmDatabaseSchema=cdmDatabaseSchema,
                       cohortDatabaseSchema = cohortDatabaseSchema,
                       cohortTable = cohortTable,
                       tempEmulationSchema = getOption("sqlRenderTempEmulationSchema"),
                       databaseId = databaseId,
		                   outputFolderIR=outputFolderIR,
                       breaks_age = breaks_age
                      ) 

# Execute the incidence rate analysis stratified by the Charlson Comorbidity Index.
# The CCI is computed as in FeatureExtraction v3.2.0 (all condition_era on/before index).
# Set Grouping = FALSE here because execute_IR above has already (re)built the grouping cohorts.
execute_IR_charlson(connectionDetails = connectionDetails,
                    cdmDatabaseSchema = cdmDatabaseSchema,
                    cohortDatabaseSchema = cohortDatabaseSchema,
                    cohortTable = cohortTable,
                    tempEmulationSchema = getOption("sqlRenderTempEmulationSchema"),
                    databaseId = databaseId,
                    outputFolderIR = outputFolderIR,
                    breaks_charlson = breaks_charlson,
                    Grouping = FALSE)


# if you want to view the shiny app locally, uncomment the following section
# CohortDiagnostics::createMergedResultsFile(
#  dataFolder = outputFolder,
#  sqliteDbPath = file.path(outputFolder,
#                           "MergedCohortDiagnosticsData.sqlite")
#)
#CohortDiagnostics::launchDiagnosticsExplorer(dataFolder = outputFolder)

