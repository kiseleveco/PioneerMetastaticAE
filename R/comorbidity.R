# Copyright 2023 PIONEER
#
# This file is part of PioneerMetastaticAE
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Execute the incidence rate analysis stratified by prior comorbidities
#'
#' @details
#' This function stratifies the adverse-event incidence rates of the study (target) cohorts
#' by the presence / absence of prior comorbidities. Each comorbidity is a cohort from the
#' OHDSI PhenotypeLibrary (https://ohdsi.github.io/PhenotypeLibrary/). The list of
#' comorbidities is ALWAYS read from \code{inst/settings/comorbidities.csv} (columns
#' \code{cohortId}, \code{comorbidity_name}); nothing is hard-coded.
#'
#' A patient is counted as having a comorbidity when there is at least one entry into the
#' comorbidity cohort with \code{cohort_start_date} strictly before the study cohort index
#' date (all prior history, no lower bound). Each comorbidity is used as an independent
#' present / absent stratifier (no cross-combinations, no interaction with age or Charlson).
#'
#' Outputs are written so that they flow through the package distribution scripts
#' (\code{extras/merge_csv.r}, \code{extras/merge_tte.r}, \code{extras/merge_mcf.r} with
#' \code{extras/arrange_surv_objects.r} / \code{extras/arrange_mcf_objects.r}) unchanged:
#' incidence rates as long-format CSVs, and one TTE and one MCF RDS file per comorbidity
#' (each a flat list of one object per target/outcome, so the arrange_* matrices are not
#' overwritten).
#'
#' @param connectionDetails                   An object of type \code{connectionDetails}.
#' @param cdmDatabaseSchema                   Schema holding the OMOP CDM data.
#' @param cohortDatabaseSchema                Schema with write privileges for the cohort tables.
#' @param cohortTable                         The study cohort table.
#' @param comorbidityCohortTable              Table that will hold the PhenotypeLibrary comorbidity
#'                                            cohorts. Defaults to \code{<cohortTable>_comorbidity}.
#' @param tempEmulationSchema                 Schema for temp-table emulation, where required.
#' @param databaseId                          A short string identifying the database.
#' @param outputFolderIR                      Local folder to place results.
#' @param generateComorbidityCohorts          Whether to (re)generate the comorbidity cohorts from
#'                                            the PhenotypeLibrary before the analysis. Set to FALSE
#'                                            to reuse an already-populated \code{comorbidityCohortTable}.
#'
#' @export
execute_IR_comorbidity <- function(connectionDetails,
                                   cdmDatabaseSchema,
                                   cohortDatabaseSchema = cdmDatabaseSchema,
                                   cohortTable = "cohort",
                                   comorbidityCohortTable = paste0(cohortTable, "_comorbidity"),
                                   tempEmulationSchema = getOption("sqlRenderTempEmulationSchema"),
                                   databaseId = databaseId,
                                   outputFolderIR = outputFolderIR,
                                   generateComorbidityCohorts = TRUE)
{
if (!file.exists(outputFolderIR)) {
  dir.create(outputFolderIR, recursive = TRUE)
}

# The analysis settings (target/outcome pairs) and the comorbidity list.
IRsettings    = read.csv("./inst/settings/IRsettings.csv")
comorbidities = read.csv("./inst/settings/comorbidities.csv")
comorbidity_ids   = as.integer(comorbidities$cohortId)
comorbidity_names = as.character(comorbidities$comorbidity_name)

##########################################################################################
# 1. Generate the PhenotypeLibrary comorbidity cohorts (optional)
##########################################################################################
if (generateComorbidityCohorts) {
  message("Fetching ", length(comorbidity_ids), " PhenotypeLibrary comorbidity cohorts")
  cohortDefinitionSet <- PhenotypeLibrary::getPlCohortDefinitionSet(cohortIds = comorbidity_ids)
  cohortTableNames    <- CohortGenerator::getCohortTableNames(cohortTable = comorbidityCohortTable)
  CohortGenerator::createCohortTables(connectionDetails    = connectionDetails,
                                      cohortDatabaseSchema = cohortDatabaseSchema,
                                      cohortTableNames     = cohortTableNames,
                                      incremental          = FALSE)
  CohortGenerator::generateCohortSet(connectionDetails    = connectionDetails,
                                     cdmDatabaseSchema     = cdmDatabaseSchema,
                                     tempEmulationSchema   = tempEmulationSchema,
                                     cohortDatabaseSchema  = cohortDatabaseSchema,
                                     cohortTableNames      = cohortTableNames,
                                     cohortDefinitionSet   = cohortDefinitionSet,
                                     incremental           = FALSE)
}

conn = DatabaseConnector::connect(connectionDetails)

##########################################################################################
# SQL scripts (target / outcome extraction - identical to execute_IR)
##########################################################################################
sql_all_outcomes='SELECT T.COHORT_DEFINITION_ID as target_id, O.COHORT_DEFINITION_ID as outcome_id, T.SUBJECT_ID ,T.COHORT_START_DATE,O.COHORT_START_DATE as OUTCOME_START_DATE,T.COHORT_END_DATE, p.YEAR_OF_BIRTH
FROM @cohortDatabaseSchema.@cohortTable T
JOIN @cohortDatabaseSchema.@cohortTable O ON T.subject_id = O.subject_id
JOIN @cdmDatabaseSchema.person p ON T.subject_id = p.person_id
WHERE T.cohort_definition_id = @target_id
AND YEAR(T.COHORT_START_DATE) - p.YEAR_OF_BIRTH > 18
  AND O.cohort_definition_id = @outcome_id
  AND O.cohort_start_date >= T.cohort_start_date
  AND O.cohort_end_date <= T.cohort_end_date'

sql_first_outcome='SELECT T.COHORT_DEFINITION_ID as target_id, O.COHORT_DEFINITION_ID as outcome_id, T.SUBJECT_ID ,T.COHORT_START_DATE,O.COHORT_START_DATE as OUTCOME_START_DATE,T.COHORT_END_DATE, p.YEAR_OF_BIRTH
FROM @cohortDatabaseSchema.@cohortTable T
JOIN @cohortDatabaseSchema.@cohortTable O ON T.subject_id = O.subject_id
JOIN @cdmDatabaseSchema.person p ON T.subject_id = p.person_id
WHERE T.cohort_definition_id = @target_id
  AND O.cohort_definition_id = @outcome_id
  AND O.cohort_start_date >= T.cohort_start_date
  AND YEAR(T.COHORT_START_DATE) - p.YEAR_OF_BIRTH > 18
  AND O.cohort_start_date = (
    SELECT MIN(cohort_start_date)
    FROM @cohortDatabaseSchema.@cohortTable
    WHERE subject_id = O.subject_id
      AND cohort_definition_id = @outcome_id
      AND cohort_start_date >= T.cohort_start_date
      AND cohort_end_date <= T.cohort_end_date
  )'

sql_target ='
SELECT T.COHORT_DEFINITION_ID as target_id, T.SUBJECT_ID, T.COHORT_START_DATE, T.COHORT_END_DATE, p.YEAR_OF_BIRTH
FROM @cohortDatabaseSchema.@cohortTable T
JOIN @cdmDatabaseSchema.person p ON T.subject_id = p.person_id
WHERE T.COHORT_DEFINITION_ID=@target_id AND YEAR(T.COHORT_START_DATE) - p.YEAR_OF_BIRTH > 18'

sql_target_chronic= '
SELECT T.COHORT_DEFINITION_ID as target_id, T.SUBJECT_ID, T.COHORT_START_DATE, T.COHORT_END_DATE, p.YEAR_OF_BIRTH
FROM @cohortDatabaseSchema.@cohortTable T
JOIN @cdmDatabaseSchema.person p ON T.subject_id = p.person_id
WHERE T.COHORT_DEFINITION_ID = @target_id
AND YEAR(T.COHORT_START_DATE) - p.YEAR_OF_BIRTH > 18
AND NOT EXISTS (
    SELECT 1
    FROM @cohortDatabaseSchema.@cohortTable AS outcome
    WHERE outcome.SUBJECT_ID = T.SUBJECT_ID
    AND outcome.COHORT_DEFINITION_ID = @outcome_id
    AND outcome.COHORT_START_DATE < T.COHORT_START_DATE
);'

# Per-subject comorbidity presence flags for one target cohort. One column per comorbidity
# (cm_<id>) = 1 if the subject has >=1 comorbidity-cohort entry strictly before the study
# index date, else 0. CASE columns are built from the CSV-derived ids - nothing hard-coded.
case_cols <- paste(sprintf("  MAX(CASE WHEN c.cohort_definition_id = %d THEN 1 ELSE 0 END) AS cm_%d",
                           comorbidity_ids, comorbidity_ids), collapse = ",\n")
sql_presence <- paste0(
"SELECT t.subject_id AS subject_id,\n", case_cols, "\n",
"FROM @cohortDatabaseSchema.@cohortTable t\n",
"LEFT JOIN @cohortDatabaseSchema.@comorbidityCohortTable c\n",
"  ON t.subject_id = c.subject_id\n",
"  AND c.cohort_start_date < t.cohort_start_date\n",
"WHERE t.cohort_definition_id = @target_id\n",
"GROUP BY t.subject_id;")

##########################################################################################
# SQL helper functions (mirror execute_IR / execute_IR_charlson)
##########################################################################################
sql_translate_ae=function(sql,target_id,outcome_id)
{
sql_rendered= SqlRender::render(sql,
                cdmDatabaseSchema=cdmDatabaseSchema,
                cohortDatabaseSchema=cohortDatabaseSchema,
                cohortTable=cohortTable,
                target_id=target_id,
                outcome_id=outcome_id)
return(SqlRender::translate(sql_rendered,targetDialect = DBMS))
}

sql_translate_target=function(sql,target_id)
{
sql_rendered= SqlRender::render(sql,
                cdmDatabaseSchema=cdmDatabaseSchema,
                cohortDatabaseSchema=cohortDatabaseSchema,
                cohortTable=cohortTable,
                target_id=target_id)
return(SqlRender::translate(sql_rendered,targetDialect = DBMS))
}

get_data=function(sql)
{
suppressMessages(as.data.frame(DatabaseConnector::querySql(conn,sql)))
}

# Per-subject comorbidity flags for a single target cohort (lower-cased column names).
get_comorbidity_flags=function(target_id)
{
  sql_rendered = SqlRender::render(sql_presence,
                                   cohortDatabaseSchema   = cohortDatabaseSchema,
                                   cohortTable            = cohortTable,
                                   comorbidityCohortTable = comorbidityCohortTable,
                                   target_id              = target_id)
  flags = get_data(SqlRender::translate(sql_rendered, targetDialect = DBMS))
  names(flags) = tolower(names(flags))
  return(flags)
}

# Attach the wide flag table to a per-subject data frame; missing -> absent (0).
attach_flags=function(ds, flags)
{
  suppressMessages(ds <- dplyr::left_join(ds, flags, by = "subject_id"))
  cm_cols <- grep("^cm_", names(ds), value = TRUE)
  for (cc in cm_cols) ds[[cc]][is.na(ds[[cc]])] <- 0
  return(ds)
}

##########################################################################################
# Data preparation (per subject; no age / Charlson), mirrors execute_IR
##########################################################################################
data_prep_chronic <- function(ae, target){
  df_ae <- mutate(ae,
                   target_id              = ae$TARGET_ID,
                   outcome_id             = ae$OUTCOME_ID,
                   subject_id             = ae$SUBJECT_ID,
                   time_to_outcome_years  = as.numeric(difftime(ae$OUTCOME_START_DATE, ae$COHORT_START_DATE, units = "days"))/365.25,
                   index_year_target      = lubridate::year(ae$COHORT_START_DATE)) %>%
              group_by(target_id, outcome_id, subject_id, time_to_outcome_years, index_year_target) %>%
              summarise(count = n())
  df_target <-  mutate(target,
                         target_id         = target$TARGET_ID,
                         subject_id        = target$SUBJECT_ID,
                         follow_up_years   = as.numeric(difftime( target$COHORT_END_DATE, target$COHORT_START_DATE, units = "days"))/365.25,
                         index_year_target = lubridate::year(target$COHORT_START_DATE)) %>%
                      group_by(target_id, subject_id, index_year_target) %>%
                  summarise(follow_up_time = sum(follow_up_years))
  out_id  <- unique(df_ae$outcome_id)
  suppressMessages(
  df_comb <- (left_join(df_target, df_ae) %>%
                mutate(TAR = if_else(!is.na(time_to_outcome_years), time_to_outcome_years, follow_up_time),
                       outcome_id = out_id) %>%
                mutate(count = if_else(is.na(count), 0, count))))
  return(as.data.frame(df_comb))
}

data_prep_episodic <- function(ae, target){
  df_ae <- mutate(ae,
                   target_id              = ae$TARGET_ID,
                   outcome_id             = ae$OUTCOME_ID,
                   subject_id             = ae$SUBJECT_ID,
                   time_to_outcome_years  = as.numeric(difftime(ae$OUTCOME_START_DATE, ae$COHORT_START_DATE, units = "days"))/365.25,
                   index_year_target      = lubridate::year(ae$COHORT_START_DATE)) %>%
              group_by(target_id, outcome_id, subject_id, index_year_target) %>%
              summarise(count = n(),
                        time_to_outcome_years = min(time_to_outcome_years))
  df_target <-  mutate(target,
                         target_id         = target$TARGET_ID,
                         subject_id        = target$SUBJECT_ID,
                         follow_up_years   = as.numeric(difftime( target$COHORT_END_DATE, target$COHORT_START_DATE, units = "days"))/365.25,
                         index_year_target = lubridate::year(target$COHORT_START_DATE)) %>%
                      group_by(target_id, subject_id, index_year_target) %>%
                  summarise(follow_up_time = sum(follow_up_years))
  out_id  <- unique(df_ae$outcome_id)
  suppressMessages(
  df_comb <- (left_join(df_target, df_ae) %>%
                mutate(TAR =  follow_up_time,
                                 outcome_id = out_id) %>%
                mutate(count = if_else(is.na(count), 0, count))))
  return(as.data.frame(df_comb))
}

# Episodic data for reda::mcf (one row per target entry (event=0) and per outcome (event=1)).
mcf_episodic_prep=function(ae,target)
{
  df_ae <- mutate(ae,
                   target_id              = ae$TARGET_ID,
                   outcome_id             = ae$OUTCOME_ID,
                   subject_id             = ae$SUBJECT_ID,
                   time_to_outcome_years  = as.numeric(difftime(ae$OUTCOME_START_DATE, ae$COHORT_START_DATE, units = "days"))/365.25,
                   index_year_target      = lubridate::year(ae$COHORT_START_DATE)) %>%
              group_by(target_id, outcome_id, subject_id, time_to_outcome_years, index_year_target) %>%
              summarise(count = n())
  df_target <-  mutate(target,
                         target_id         = target$TARGET_ID,
                         subject_id        = target$SUBJECT_ID,
                         follow_up_years   = as.numeric(difftime( target$COHORT_END_DATE, target$COHORT_START_DATE, units = "days"))/365.25,
                         index_year_target = lubridate::year(target$COHORT_START_DATE)) %>%
                       group_by(target_id, subject_id, index_year_target) %>%
                  summarise(follow_up_time = sum(follow_up_years))

 df_ae_1=data.frame(cbind(  df_ae$target_id,
                            df_ae$outcome_id,
                            df_ae$subject_id,
                            df_ae$index_year_target,
                            df_ae$time_to_outcome_years,
                            event=1))
names(df_ae_1)=c("target_id","outcome_id","subject_id","Year","time","event")

outcome_id=unique(df_ae_1$outcome_id)
 df_o_1=data.frame(cbind(   df_target$target_id ,
                            outcome_id,
                            df_target$subject_id ,
                            df_target$index_year_target,
                            df_target$follow_up_time,
                            event=0 ))
names(df_o_1)=c("target_id","outcome_id","subject_id","Year","time","event")

df=data.frame(rbind(df_o_1, df_ae_1))
return(df)
}

##########################################################################################
# Incidence-rate summaries by comorbidity presence/absence
##########################################################################################
add_poisson_ci=function(res2)
{
  CI=survival::cipoisson(res2$total_events,res2$total_time_at_risk)*1000
  if(nrow(res2)==1)
  {
    res2$lower=CI[1]
    res2$upper=CI[2]
  }else{
    res2$lower=CI[,1]
    res2$upper=CI[,2]
  }
  return(res2)
}

# Aggregate one comorbidity's present(0/1) x index-year cells (filling missing cells with 0).
aggregate_present = function(df, present)
{
  df$present <- present
  out_id     <- unique(df$outcome_id)
  suppressMessages(
  agg <- (group_by(df, target_id, present, index_year_target) %>%
            summarise(nr_events = sum(count, na.rm = TRUE), time_at_risk = sum(TAR), nr_patients = n())))
  res        <- expand.grid(unique(agg$target_id), c(0, 1), unique(agg$index_year_target))
  names(res) <- c("target_id", "present", "index_year_target")
  suppressMessages(
  res        <- (left_join(res, agg) %>%
                   mutate(nr_events    = coalesce(nr_events, rep(0, nrow(res))),
                          time_at_risk = coalesce(time_at_risk, rep(0, nrow(res))),
                          nr_patients  = coalesce(nr_patients, rep(0, nrow(res))))))
  res$outcome_id <- out_id
  return(as.data.frame(res))
}

get_IR_by_comorbidity_only=function(df)
{
  res2 <- df %>%
    group_by(target_id,outcome_id,present )%>%
    summarize(
      total_events        = sum(nr_events),
      total_time_at_risk  = sum(time_at_risk),
      total_nr_patients   = sum(nr_patients))%>%
    mutate(incidence_rate = (total_events / total_time_at_risk) * 1000)
  return(add_poisson_ci(res2))
}

get_IR_by_comorbidity_and_year=function(df)
{
  res2 <- df %>%
    group_by(target_id,outcome_id,present,index_year_target )%>%
    summarize(
      total_events        = sum(nr_events),
      total_time_at_risk  = sum(time_at_risk),
      total_nr_patients   = sum(nr_patients))%>%
    mutate(incidence_rate = (total_events / total_time_at_risk) * 1000)
  return(add_poisson_ci(res2))
}

##########################################################################################
# Writing / censoring helpers (mirror execute_IR)
##########################################################################################
write_table=function(object,csv_file)
{
if (file.exists(csv_file)) {
      write.table(object, csv_file, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
    } else 	{
      write.table(object, csv_file, sep = ",", col.names = TRUE, row.names = FALSE)
        }
}

censor_counts_IR=function(df)
{
df$total_time_at_risk[df$total_nr_patients<5 & df$total_nr_patients!=0]=NA
df$incidence_rate[df$total_nr_patients<5 | df$total_events<5]=NA
df$lower[df$total_nr_patients<5 | df$total_events<5]=NA
df$upper[df$total_nr_patients<5 | df$total_events<5]=NA
df$total_events [df$total_events<5 & df$total_events!=0]="<5"
df$total_nr_patients[df$total_nr_patients<5 & df$total_nr_patients!=0]="<5"
return(df)
}

# Compute, label, censor and append the IR tables for one comorbidity (long format).
write_incidence_rates_csv=function(df, chronic_indicator, comorbidity_id, comorbidity_name)
{
  only = get_IR_by_comorbidity_only(df)
  only$comorbidity_id = comorbidity_id
  only$comorbidity_name = comorbidity_name
  only$databaseId = databaseId
  only = censor_counts_IR(only)

  yr = get_IR_by_comorbidity_and_year(df)
  yr$comorbidity_id = comorbidity_id
  yr$comorbidity_name = comorbidity_name
  yr$databaseId = databaseId
  yr = censor_counts_IR(yr)

  if (chronic_indicator == 1) {
    write_table(only, paste0(outputFolderIR, "/chronic_IR_by_comorbidity_only.csv"))
    write_table(yr,   paste0(outputFolderIR, "/chronic_IR_by_comorbidity_and_year.csv"))
  } else {
    write_table(only, paste0(outputFolderIR, "/episodic_IR_by_comorbidity_only.csv"))
    write_table(yr,   paste0(outputFolderIR, "/episodic_IR_by_comorbidity_and_year.csv"))
  }
}

present_factor = function(x) factor(ifelse(x == 1, "present", "absent"), levels = c("absent", "present"))

# TRUE when both present and absent are observed, so a stratified fit is possible. When a
# comorbidity is all-present or all-absent in a cell, a one-level factor would break
# survfit/reda::mcf ("contrasts can be applied only to factors with 2 or more levels"),
# so those fits fall back to an unstratified (~ 1) model.
has_two_levels = function(x) length(unique(x[!is.na(x)])) >= 2

##########################################################################################
# Main loop
##########################################################################################
# One flat list per comorbidity (keyed by cohortId as character) so that each file holds a
# single object per target/outcome and arrange_tte_objects / arrange_mcf_objects do not clash.
tte_chronic  = list()
tte_episodic = list()
mcf_comorb   = list()
for (cid in as.character(comorbidity_ids)) {
  tte_chronic[[cid]]  = list()
  tte_episodic[[cid]] = list()
  mcf_comorb[[cid]]   = list()
}

target_id <- unique(IRsettings$target_id)

for(i in 1:length(target_id)){
    target <- sql_translate_target(sql = sql_target, target_id = target_id[i]) %>% get_data
    if(nrow(target)<5)
    {
    print(paste("No data for cohort:", target_id[i]))
    next
    }
    # Comorbidity presence for every subject in this target cohort (computed once per target).
    flags <- get_comorbidity_flags(target_id[i])

    OutcomeList=IRsettings[IRsettings$target_id==target_id[i],]
    for(j in 1:nrow(OutcomeList)){
        if(OutcomeList$chronic[j]==1)
        {
            target_chronic= sql_translate_ae(sql = sql_target_chronic,target_id=target_id[i], outcome_id = OutcomeList$outcome_id[j]) %>% get_data
            if(nrow(target_chronic)<5)
            {
                print(paste("Not enough data for target cohort:", target_id[i],"and outcome:",OutcomeList$outcome_id[j]))
                next
            }
        outcome <- sql_translate_ae(sql = sql_first_outcome,target_id=target_id[i], outcome_id = OutcomeList$outcome_id[j]) %>% get_data
        if(nrow(outcome)<5)
            {
                print(paste("Not enough outcomes for target cohort:", target_id[i],"and outcome:",OutcomeList$outcome_id[j]))
                next
            }
        print(paste(" Comorbidity-stratified (chronic) incidence for target:", target_id[i],"and outcome:", OutcomeList$outcome_id[j]))
        ds <- attach_flags(data_prep_chronic(ae = outcome, target = target_chronic), flags)

        for (k in seq_along(comorbidity_ids)) {
          cid <- as.character(comorbidity_ids[k]); nm <- comorbidity_names[k]
          col <- paste0("cm_", comorbidity_ids[k])
          present <- if (!is.null(ds[[col]])) ds[[col]] else rep(0, nrow(ds))

          ds_ir <- aggregate_present(ds, present)
          write_incidence_rates_csv(ds_ir, 1, comorbidity_ids[k], nm)

          grp <- present_factor(present)
          sf <- if (has_two_levels(present)) survfit(Surv(ds$TAR, ds$count) ~ grp, data = ds) else survfit(Surv(ds$TAR, ds$count) ~ 1, data = ds)
          sf$target = target_id[i]; sf$outcome = OutcomeList$outcome_id[j]
          sf$comorbidity_id = comorbidity_ids[k]; sf$comorbidity_name = nm; sf$databaseId = databaseId
          tte_chronic[[cid]][[length(tte_chronic[[cid]]) + 1]] <- sf
        }
        }

        if(OutcomeList$episodic[j]==1)
        {
           outcome <- sql_translate_ae(sql = sql_all_outcomes,target_id=target_id[i], outcome_id = OutcomeList$outcome_id[j]) %>% get_data
           if(nrow(outcome)<5)
               {
               print(paste("Not enough outcomes for cohort:", target_id[i],"and outcome:",OutcomeList$outcome_id[j]))
               next
               }
        print(paste(" Comorbidity-stratified (episodic) incidence for target:", target_id[i],"and outcome:", OutcomeList$outcome_id[j]))
        ds <- attach_flags(data_prep_episodic(ae = outcome, target = target), flags)
        modified_count <- ifelse(ds$count> 0, 1, 0)

        ds_ep <- mcf_episodic_prep(ae = outcome, target = target)
        suppressMessages(ds_ep <- dplyr::left_join(ds_ep, flags, by = "subject_id"))
        cm_cols <- grep("^cm_", names(ds_ep), value = TRUE)
        for (cc in cm_cols) ds_ep[[cc]][is.na(ds_ep[[cc]])] <- 0

        for (k in seq_along(comorbidity_ids)) {
          cid <- as.character(comorbidity_ids[k]); nm <- comorbidity_names[k]
          col <- paste0("cm_", comorbidity_ids[k])
          present <- if (!is.null(ds[[col]])) ds[[col]] else rep(0, nrow(ds))

          ds_ir <- aggregate_present(ds, present)
          write_incidence_rates_csv(ds_ir, 0, comorbidity_ids[k], nm)

          grp <- present_factor(present)
          sf <- if (has_two_levels(present)) survfit(Surv(ds$TAR, modified_count) ~ grp, data = ds) else survfit(Surv(ds$TAR, modified_count) ~ 1, data = ds)
          sf$target = target_id[i]; sf$outcome = OutcomeList$outcome_id[j]
          sf$comorbidity_id = comorbidity_ids[k]; sf$comorbidity_name = nm; sf$databaseId = databaseId
          tte_episodic[[cid]][[length(tte_episodic[[cid]]) + 1]] <- sf

          ds_ep$present <- present_factor(ds_ep[[col]])
          mcf1 <- if (has_two_levels(ds_ep[[col]])) reda::mcf(Recur(id = subject_id, time = time, event = event) ~ present, data = ds_ep) else reda::mcf(Recur(id = subject_id, time = time, event = event) ~ 1, data = ds_ep)
          m <- mcf1@MCF
          m$target = target_id[i]; m$outcome = OutcomeList$outcome_id[j]
          m$comorbidity_id = comorbidity_ids[k]; m$comorbidity_name = nm; m$databaseId = databaseId
          mcf_comorb[[cid]][[length(mcf_comorb[[cid]]) + 1]] <- m
        }
        }
    }
   }

##########################################################################################
# Write one TTE + one MCF file per comorbidity (filenames carry "tte"/"mcf" so the
# distribution scripts pick them up and arrange them into per-comorbidity matrices).
##########################################################################################
for (k in seq_along(comorbidity_ids)) {
  cid <- as.character(comorbidity_ids[k]); id <- comorbidity_ids[k]
  saveRDS(tte_chronic[[cid]],  paste0(outputFolderIR, "/chronic_tte_by_comorbidity_",  id, ".rds"))
  saveRDS(tte_episodic[[cid]], paste0(outputFolderIR, "/episodic_tte_by_comorbidity_", id, ".rds"))
  saveRDS(mcf_comorb[[cid]],   paste0(outputFolderIR, "/mcf_episodic_by_comorbidity_", id, ".rds"))
}

DatabaseConnector::disconnect(conn)
}
