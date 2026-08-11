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

#' Execute the incidence rate analysis stratified by the Charlson Comorbidity Index
#'
#' @details
#' This function extracts the data needed to calculate incidence rates of outcomes across
#' targets, stratified by the patient's Charlson Comorbidity Index (CCI).
#'
#' The Charlson index is computed exactly as in \code{FeatureExtraction} v3.2.0
#' (the 17-category "Romano adaptation"): for every patient the weights of the distinct
#' Charlson categories whose \code{condition_era} starts on or before the cohort index date
#' are summed. Following FeatureExtraction, no lower observation window is applied
#' (all history on or before index is considered) and a patient with no qualifying
#' condition_era receives a score of 0 (\code{missing_means_zero}).
#'
#' The continuous score is then binned with \code{breaks_charlson} and used as the sole
#' stratifier, producing incidence-rate tables (overall and by index year within each
#' Charlson category), Kaplan-Meier / time-to-event fits and recurrent-event (MCF)
#' estimates by Charlson category. No age stratification is applied.
#'
#' @param connectionDetails                   An object of type \code{connectionDetails} as created
#'                                            using the
#'                                            \code{\link[DatabaseConnector]{createConnectionDetails}}
#'                                            function in the DatabaseConnector package.
#' @param cdmDatabaseSchema                   Schema name where your patient-level data in OMOP CDM
#'                                            format resides. Note that for SQL Server, this should
#'                                            include both the database and schema name, for example
#'                                            'cdm_data.dbo'.
#' @param cohortDatabaseSchema                Schema name where intermediate data can be stored. You
#'                                            will need to have write privileges in this schema.
#' @param cohortTable                         The name of the table that holds the exposure and
#'                                            outcome cohorts used in this study.
#' @param tempEmulationSchema                 Some database platforms like Oracle and Impala do not
#'                                            truly support temp tables. To emulate temp tables,
#'                                            provide a schema with write privileges where temp tables
#'                                            can be created.
#' @param breaks_charlson                     Numeric vector of interior break points for the Charlson
#'                                            score. The default \code{c(0, 2, 4)} yields the
#'                                            categories 0, 1-2, 3-4 and 5+.
#' @param databaseId                          A short string for identifying the database (e.g.
#'                                            'Synpuf').
#' @param outputFolderIR                      Name of local folder to place results; make sure to use
#'                                            forward slashes (/).
#' @param Grouping                            Whether to (re)build the grouping cohorts before the
#'                                            analysis. Set to FALSE if \code{\link{execute_IR}} has
#'                                            already been run in the same session/output folder.
#'
#' @export
execute_IR_charlson <- function(connectionDetails,
                                cdmDatabaseSchema,
                                cohortDatabaseSchema = cdmDatabaseSchema,
                                cohortTable = "cohort",
                                tempEmulationSchema = getOption("sqlRenderTempEmulationSchema"),
                                breaks_charlson = c(0, 2, 4),
                                databaseId = databaseId,
                                outputFolderIR = outputFolderIR,
                                Grouping = FALSE)
{
# check if output folder is present and if not create it
if (!file.exists(outputFolderIR)) {
  dir.create(outputFolderIR, recursive = TRUE)
}
# the csv where the analysis settings are
IRsettings = read.csv("./inst/settings/IRsettings.csv")
# Connection details
conn = DatabaseConnector::connect(connectionDetails)

##########################################################################################
# SQL scripts (target / outcome extraction - identical to execute_IR)
##########################################################################################

# All the outcomes during the target time at risk (episodic events)
sql_all_outcomes='SELECT T.COHORT_DEFINITION_ID as target_id, O.COHORT_DEFINITION_ID as outcome_id, T.SUBJECT_ID ,T.COHORT_START_DATE,O.COHORT_START_DATE as OUTCOME_START_DATE,T.COHORT_END_DATE, p.YEAR_OF_BIRTH
FROM @cohortDatabaseSchema.@cohortTable T
JOIN @cohortDatabaseSchema.@cohortTable O ON T.subject_id = O.subject_id
JOIN @cdmDatabaseSchema.person p ON T.subject_id = p.person_id
WHERE T.cohort_definition_id = @target_id
AND YEAR(T.COHORT_START_DATE) - p.YEAR_OF_BIRTH > 18
  AND O.cohort_definition_id = @outcome_id
  AND O.cohort_start_date >= T.cohort_start_date
  AND O.cohort_end_date <= T.cohort_end_date'

# first outcome during/after target (chronic events)
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

# All the target population
sql_target ='
SELECT T.COHORT_DEFINITION_ID as target_id, T.SUBJECT_ID, T.COHORT_START_DATE, T.COHORT_END_DATE, p.YEAR_OF_BIRTH
FROM @cohortDatabaseSchema.@cohortTable T
JOIN @cdmDatabaseSchema.person p ON T.subject_id = p.person_id
WHERE T.COHORT_DEFINITION_ID=@target_id AND YEAR(T.COHORT_START_DATE) - p.YEAR_OF_BIRTH > 18'

# For the chronic population patients having the event before the cohort start date are excluded
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

sql_grouping='
WITH new_cohort AS (
  SELECT
    SUBJECT_ID,
    MIN(COHORT_START_DATE) as min_start_date,
    MAX(COHORT_END_DATE) as max_end_date
  FROM
    @cohortDatabaseSchema.@cohortTable
  WHERE
    COHORT_DEFINITION_ID IN (@cohort_list)
  GROUP BY
    SUBJECT_ID
)
INSERT INTO @cohortDatabaseSchema.@cohortTable (COHORT_DEFINITION_ID, SUBJECT_ID, COHORT_START_DATE, COHORT_END_DATE)
SELECT
  @new_cohort_id,
  SUBJECT_ID,
  min_start_date,
  max_end_date
FROM
  new_cohort;'

sql_drop_cohort_id= "DELETE FROM @cohortDatabaseSchema.@cohortTable
WHERE COHORT_DEFINITION_ID = @cohort_id;"

##########################################################################################
# Charlson Comorbidity Index - FeatureExtraction v3.2.0 ("Romano adaptation")
#
# The 17 diagnostic categories, their weights and ancestor concept ids are taken verbatim
# from FeatureExtraction v3.2.0 inst/sql/sql_server/CharlsonIndex.sql. As in
# FeatureExtraction, a category contributes its weight once if the patient has any
# condition_era (descendant of the ancestor) starting on or before the cohort index date;
# no lower observation window is applied.
##########################################################################################

# id | name | weight | ancestor_concept_id(s)
charlson_categories <- list(
  list(id = 1,  name = "Myocardial infarction",                weight = 1, ancestors = "4329847"),
  list(id = 2,  name = "Congestive heart failure",             weight = 1, ancestors = "316139"),
  list(id = 3,  name = "Peripheral vascular disease",          weight = 1, ancestors = "321052"),
  list(id = 4,  name = "Cerebrovascular disease",              weight = 1, ancestors = "381591,434056"),
  list(id = 5,  name = "Dementia",                             weight = 1, ancestors = "4182210"),
  list(id = 6,  name = "Chronic pulmonary disease",            weight = 1, ancestors = "4063381"),
  list(id = 7,  name = "Rheumatologic disease",                weight = 1, ancestors = "257628,134442,80800,80809,256197,255348"),
  list(id = 8,  name = "Peptic ulcer disease",                 weight = 1, ancestors = "4247120"),
  list(id = 9,  name = "Mild liver disease",                   weight = 1, ancestors = "4064161,4212540"),
  list(id = 10, name = "Diabetes (mild to moderate)",          weight = 1, ancestors = "201820"),
  list(id = 11, name = "Diabetes with chronic complications",  weight = 2, ancestors = "443767,442793"),
  list(id = 12, name = "Hemoplegia or paralegia",              weight = 2, ancestors = "192606,374022"),
  list(id = 13, name = "Renal disease",                        weight = 2, ancestors = "4030518"),
  list(id = 14, name = "Any malignancy",                       weight = 2, ancestors = "443392"),
  list(id = 15, name = "Moderate to severe liver disease",     weight = 3, ancestors = "4245975,4029488,192680,24966"),
  list(id = 16, name = "Metastatic solid tumor",               weight = 6, ancestors = "432851"),
  list(id = 17, name = "AIDS",                                 weight = 6, ancestors = "439727")
)

# Build the SQL that (re)creates and populates the #charlson_scoring and #charlson_concepts
# temp tables. Executed once per connection.
scoring_values <- paste(
  vapply(charlson_categories, function(c)
    sprintf("(%d, '%s', %d)", c$id, gsub("'", "''", c$name), c$weight),
    character(1)),
  collapse = ",\n\t")

concepts_selects <- paste(
  vapply(charlson_categories, function(c)
    sprintf("SELECT %d AS diag_category_id, descendant_concept_id\nFROM @cdmDatabaseSchema.concept_ancestor\nWHERE ancestor_concept_id IN (%s)",
            c$id, c$ancestors),
    character(1)),
  collapse = "\nUNION ALL\n")

sql_charlson_setup <- paste0(
"IF OBJECT_ID('tempdb..#charlson_concepts', 'U') IS NOT NULL DROP TABLE #charlson_concepts;
CREATE TABLE #charlson_concepts (diag_category_id INT, concept_id INT);

IF OBJECT_ID('tempdb..#charlson_scoring', 'U') IS NOT NULL DROP TABLE #charlson_scoring;
CREATE TABLE #charlson_scoring (diag_category_id INT, diag_category_name VARCHAR(255), weight INT);

INSERT INTO #charlson_scoring (diag_category_id, diag_category_name, weight)
VALUES\n\t", scoring_values, ";

INSERT INTO #charlson_concepts (diag_category_id, concept_id)\n", concepts_selects, ";")

# Per-subject Charlson score for one target cohort. Missing means zero, so subjects with no
# qualifying condition_era simply do not appear here and are coalesced to 0 downstream.
sql_charlson_scores <-
"SELECT temp.subject_id AS subject_id, SUM(temp.weight) AS charlson_score
FROM (
  SELECT DISTINCT charlson_scoring.diag_category_id,
    charlson_scoring.weight,
    cohort.subject_id AS subject_id
  FROM @cohortDatabaseSchema.@cohortTable cohort
  INNER JOIN @cdmDatabaseSchema.condition_era condition_era
    ON cohort.subject_id = condition_era.person_id
  INNER JOIN #charlson_concepts charlson_concepts
    ON condition_era.condition_concept_id = charlson_concepts.concept_id
  INNER JOIN #charlson_scoring charlson_scoring
    ON charlson_concepts.diag_category_id = charlson_scoring.diag_category_id
  WHERE condition_era.condition_era_start_date <= cohort.cohort_start_date
    AND cohort.cohort_definition_id = @target_id
) temp
GROUP BY temp.subject_id;"

sql_charlson_cleanup <-
"IF OBJECT_ID('tempdb..#charlson_concepts', 'U') IS NOT NULL DROP TABLE #charlson_concepts;
IF OBJECT_ID('tempdb..#charlson_scoring', 'U') IS NOT NULL DROP TABLE #charlson_scoring;"

##########################################################################################
# Utility / SQL helper functions (mirrors execute_IR)
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

sql_execute_grouping=function(sql,cohortDatabaseSchema,cohortTable,cohort_list,new_cohort_id)
{
sql_rendered= SqlRender::render(sql,
          cohortDatabaseSchema=cohortDatabaseSchema,
          cohortTable=cohortTable,
          cohort_list=cohort_list,
          new_cohort_id=new_cohort_id)
 DatabaseConnector::executeSql(conn, SqlRender::translate(sql_rendered,targetDialect = DBMS))
}

sql_execute_drop_cohort_id=function(cohortDatabaseSchema,cohortTable,cohort_id)
{
sql_rendered=SqlRender::render(sql = sql_drop_cohort_id,
                              cohortDatabaseSchema=cohortDatabaseSchema,
                              cohortTable=cohortTable,
                              cohort_id=cohort_id)
  DatabaseConnector::executeSql(conn, SqlRender::translate(sql_rendered,targetDialect = DBMS))
}

get_data=function(sql)
{
suppressMessages(as.data.frame(DatabaseConnector::querySql(conn,sql)))
}

# Build the Charlson scoring / concept temp tables once.
build_charlson_tables=function()
{
  sql_rendered = SqlRender::render(sql_charlson_setup, cdmDatabaseSchema = cdmDatabaseSchema)
  DatabaseConnector::executeSql(conn, SqlRender::translate(sql_rendered, targetDialect = DBMS))
}

# Drop the Charlson temp tables.
drop_charlson_tables=function()
{
  DatabaseConnector::executeSql(conn, SqlRender::translate(sql_charlson_cleanup, targetDialect = DBMS))
}

# Per-subject Charlson score for a single target cohort.
get_charlson_scores=function(target_id)
{
  sql_rendered = SqlRender::render(sql_charlson_scores,
                                   cdmDatabaseSchema = cdmDatabaseSchema,
                                   cohortDatabaseSchema = cohortDatabaseSchema,
                                   cohortTable = cohortTable,
                                   target_id = target_id)
  get_data(SqlRender::translate(sql_rendered, targetDialect = DBMS))
}

# Attach the Charlson score (coalesced to 0) to a raw cohort data frame.
attach_charlson=function(df, charlson)
{
  suppressMessages(df <- dplyr::left_join(df, charlson, by = "SUBJECT_ID"))
  df$CHARLSON_SCORE[is.na(df$CHARLSON_SCORE)] <- 0
  return(df)
}

##########################################################################################
# Data preparation (chronic / episodic) - carries the Charlson score as a grouping key
##########################################################################################

data_prep_chronic <- function(ae, target, charlson){
  ae     <- attach_charlson(ae, charlson)
  target <- attach_charlson(target, charlson)

  df_ae <- mutate(ae,
                   target_id              = ae$TARGET_ID,
                   outcome_id             = ae$OUTCOME_ID,
                   subject_id             = ae$SUBJECT_ID,
                   charlson_score         = ae$CHARLSON_SCORE,
                   time_to_outcome_years  = as.numeric(difftime(ae$OUTCOME_START_DATE, ae$COHORT_START_DATE, units = "days"))/365.25,
                   index_year_target      = lubridate::year(ae$COHORT_START_DATE)) %>%
              group_by(target_id, outcome_id, subject_id, charlson_score, time_to_outcome_years, index_year_target) %>%
              summarise(count = n())

  df_target <-  mutate(target,
                         target_id         = target$TARGET_ID,
                         subject_id        = target$SUBJECT_ID,
                         charlson_score    = target$CHARLSON_SCORE,
                         follow_up_years   = as.numeric(difftime( target$COHORT_END_DATE, target$COHORT_START_DATE, units = "days"))/365.25,
                         index_year_target = lubridate::year(target$COHORT_START_DATE)) %>%
                      group_by(target_id, subject_id, charlson_score, index_year_target) %>%
                  summarise(follow_up_time = sum(follow_up_years))

  out_id  <- unique(df_ae$outcome_id)
  suppressMessages(
  df_comb <- (left_join(df_target, df_ae) %>%
                mutate(TAR = if_else(!is.na(time_to_outcome_years), time_to_outcome_years, follow_up_time),
                       outcome_id = out_id) %>%
                mutate(count = if_else(is.na(count), 0, count))))

  return(as.data.frame(df_comb))
}

data_prep_episodic <- function(ae, target, charlson){
  ae     <- attach_charlson(ae, charlson)
  target <- attach_charlson(target, charlson)

  df_ae <- mutate(ae,
                   target_id              = ae$TARGET_ID,
                   outcome_id             = ae$OUTCOME_ID,
                   subject_id             = ae$SUBJECT_ID,
                   charlson_score         = ae$CHARLSON_SCORE,
                   time_to_outcome_years  = as.numeric(difftime(ae$OUTCOME_START_DATE, ae$COHORT_START_DATE, units = "days"))/365.25,
                   index_year_target      = lubridate::year(ae$COHORT_START_DATE)) %>%
              group_by(target_id, outcome_id, subject_id, charlson_score, index_year_target) %>%
              summarise(count = n(),
                        time_to_outcome_years = min(time_to_outcome_years))

  df_target <-  mutate(target,
                         target_id         = target$TARGET_ID,
                         subject_id        = target$SUBJECT_ID,
                         charlson_score    = target$CHARLSON_SCORE,
                         follow_up_years   = as.numeric(difftime( target$COHORT_END_DATE, target$COHORT_START_DATE, units = "days"))/365.25,
                         index_year_target = lubridate::year(target$COHORT_START_DATE)) %>%
                      group_by(target_id, subject_id, charlson_score, index_year_target) %>%
                  summarise(follow_up_time = sum(follow_up_years))

  out_id  <- unique(df_ae$outcome_id)
  suppressMessages(
  df_comb <- (left_join(df_target, df_ae) %>%
                mutate(TAR =  follow_up_time,
                                 outcome_id = out_id) %>%
                mutate(count = if_else(is.na(count), 0, count))))

  return(as.data.frame(df_comb))
}

# Bin the Charlson score. Default breaks_charlson = c(0,2,4) -> "0","(0,2]","(2,4]","(4,Inf]".
charlson_cut = function(score, breaks_charlson)
{
  cut(score, breaks = c(-Inf, breaks_charlson, Inf))
}

##########################################################################################
# Cut / aggregate at the desired Charlson categories (and index year)
##########################################################################################
cut_and_aggregate = function(df)
{
  df                <- mutate(df, charlson_cat = charlson_cut(charlson_score, breaks_charlson))
  out_id            <- unique(df$outcome_id)
  suppressMessages(
  df <- (group_by(df, target_id, charlson_cat, index_year_target) %>%
           summarise(nr_events = sum(count, na.rm = TRUE), time_at_risk = sum(TAR), nr_patients = n())))

  # Include possibly missing groups
  res        <- expand.grid(unique(df$target_id),
                            unique(df$charlson_cat),
                            unique(df$index_year_target))
  names(res) <- c("target_id", "charlson_cat", "index_year_target")
    suppressMessages(
  res        <- (left_join(res, df) %>%
                   mutate(nr_events    = coalesce(nr_events, rep(0, nrow(res))),
                          time_at_risk = coalesce(time_at_risk, rep(0, nrow(res))),
                          nr_patients  = coalesce(nr_patients, rep(0, nrow(res))))))
  res$outcome_id <- out_id
  return(as.data.frame(res))
}

# Episodic data preparation conforming to reda mcf function (adds a Charlson_Group column)
mcf_episodic_prep=function(ae,target,charlson)
{
  ae     <- attach_charlson(ae, charlson)
  target <- attach_charlson(target, charlson)

  df_ae <- mutate(ae,
                   target_id              = ae$TARGET_ID,
                   outcome_id             = ae$OUTCOME_ID,
                   subject_id             = ae$SUBJECT_ID,
                   charlson_score         = ae$CHARLSON_SCORE,
                   time_to_outcome_years  = as.numeric(difftime(ae$OUTCOME_START_DATE, ae$COHORT_START_DATE, units = "days"))/365.25,
                   index_year_target      = lubridate::year(ae$COHORT_START_DATE)) %>%
              group_by(target_id, outcome_id, subject_id, charlson_score, time_to_outcome_years, index_year_target) %>%
              summarise(count = n())

  df_target <-  mutate(target,
                         target_id         = target$TARGET_ID,
                         subject_id        = target$SUBJECT_ID,
                         charlson_score    = target$CHARLSON_SCORE,
                         follow_up_years   = as.numeric(difftime( target$COHORT_END_DATE, target$COHORT_START_DATE, units = "days"))/365.25,
                         index_year_target = lubridate::year(target$COHORT_START_DATE)) %>%
                       group_by(target_id, subject_id, charlson_score, index_year_target) %>%
                  summarise(follow_up_time = sum(follow_up_years))

 df_ae_1=data.frame(cbind(  df_ae$target_id,
                            df_ae$outcome_id,
                            df_ae$subject_id,
                            df_ae$charlson_score,
                            df_ae$index_year_target,
                            df_ae$time_to_outcome_years,
                            event=1))
names(df_ae_1)=c("target_id","outcome_id","subject_id","Charlson_Score","Year","time","event")

outcome_id=unique(df_ae_1$outcome_id)
 df_o_1=data.frame(cbind(   df_target$target_id ,
                            outcome_id,
                            df_target$subject_id ,
                            df_target$charlson_score,
                            df_target$index_year_target,
                            df_target$follow_up_time,
                            event=0 ))
names(df_o_1)=c("target_id","outcome_id","subject_id","Charlson_Score","Year","time","event")

df=data.frame(rbind(df_o_1, df_ae_1))
df$Charlson_Group= charlson_cut(df$Charlson_Score, breaks_charlson)
return(df)
}

##########################################################################################
# Incidence-rate summaries by Charlson category
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

get_IR_by_charlson_only=function(df)
{
  res2 <- df %>%
    group_by(target_id,outcome_id,charlson_cat )%>%
    summarize(
      total_events        = sum(nr_events),
      total_time_at_risk  = sum(time_at_risk),
      total_nr_patients   = sum(nr_patients))%>%
    mutate(incidence_rate = (total_events / total_time_at_risk) * 1000)
  return(add_poisson_ci(res2))
}

get_IR_by_charlson_and_year=function(df)
{
  res2 <- df %>%
    group_by(target_id,outcome_id,charlson_cat,index_year_target )%>%
    summarize(
      total_events        = sum(nr_events),
      total_time_at_risk  = sum(time_at_risk),
      total_nr_patients   = sum(nr_patients))%>%
    mutate(incidence_rate = (total_events / total_time_at_risk) * 1000)
  return(add_poisson_ci(res2))
}

##########################################################################################
# Writing / censoring helpers (mirrors execute_IR)
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

write_incidence_rates_csv=function(df,chronic_indicator)
{
IR_by_charlson_only=get_IR_by_charlson_only(df)
IR_by_charlson_only$databaseId=databaseId
IR_by_charlson_only=censor_counts_IR(IR_by_charlson_only)

IR_by_charlson_and_year=get_IR_by_charlson_and_year(df)
IR_by_charlson_and_year$databaseId=databaseId
IR_by_charlson_and_year=censor_counts_IR(IR_by_charlson_and_year)

    if(chronic_indicator==1)
    {
    write_table(IR_by_charlson_only,     paste0(outputFolderIR,"/chronic_IR_by_charlson_only.csv"))
    write_table(IR_by_charlson_and_year, paste0(outputFolderIR,"/chronic_IR_by_charlson_and_year.csv"))
    }
    if (chronic_indicator==0)
    {
    write_table(IR_by_charlson_only,     paste0(outputFolderIR,"/episodic_IR_by_charlson_only.csv"))
    write_table(IR_by_charlson_and_year, paste0(outputFolderIR,"/episodic_IR_by_charlson_and_year.csv"))
    }
}

create_grouping_cohorts=function()
{
groups=read.csv(file="./inst/settings/groups.csv")
for(i in 1:nrow(groups))
{
groups$cohort_list[i]=gsub(";",",",groups$cohort_list[i])
sql_execute_drop_cohort_id(cohortDatabaseSchema,
                            cohortTable,
                            groups$new_cohort_id[i])
sql_execute_grouping(sql_grouping,
                      cohortDatabaseSchema,
                      cohortTable,
                      groups$cohort_list[i],
                      groups$new_cohort_id[i])
}
}

##########################################################################################
# Main loop
##########################################################################################
if(Grouping==TRUE)
{
create_grouping_cohorts()
}

# Build the Charlson concept / scoring temp tables once for the whole run.
build_charlson_tables()

target_id  <- unique(IRsettings$target_id)
counter1=0   # counter for chronic events
counter2=0   # counter for episodic events
s1_charlson=NULL
s2_charlson=NULL
mcf_episodic_by_charlson=NULL

for(i in 1:length(target_id)){
    target <- sql_translate_target(sql = sql_target, target_id = target_id[i]) %>% get_data
    if(nrow(target)<5)
    {
    print(paste("No data for cohort:", target_id[i]))
    next
    }
    # Charlson score for every subject in this target cohort (computed once per target)
    charlson <- get_charlson_scores(target_id[i])

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
        counter1=counter1+1
        print(paste( " Calculating Charlson-stratified incidence for target:", target_id[i],"and outcome:", OutcomeList$outcome_id[j]))
        ds <- data_prep_chronic(ae = outcome, target = target_chronic, charlson = charlson)
        Charlson_Group= charlson_cut(ds$charlson_score, breaks_charlson)

        s1_charlson[[counter1]] <- survfit(Surv(ds$TAR, ds$count) ~  Charlson_Group,   data = ds)
        s1_charlson[[counter1]]$target      =target_id[i]
        s1_charlson[[counter1]]$outcome     =OutcomeList$outcome_id[j]
        s1_charlson[[counter1]]$databaseId  =databaseId

        ds_ir     <- cbind(cut_and_aggregate(ds),databaseId)
        write_incidence_rates_csv(ds_ir,1)
        }

        if(OutcomeList$episodic[j]==1)
        {
           outcome <- sql_translate_ae(sql = sql_all_outcomes,target_id=target_id[i], outcome_id = OutcomeList$outcome_id[j]) %>% get_data
           if(nrow(outcome)<5)
               {
               print(paste("Not enough outcomes for cohort:", target_id[i],"and outcome:",OutcomeList$outcome_id[j]))
               next
               }
        counter2=counter2+1
        print(paste( " Calculating Charlson-stratified incidence for target:", target_id[i],"and outcome:", OutcomeList$outcome_id[j]))
        ds <- data_prep_episodic(ae = outcome, target = target, charlson = charlson)
        Charlson_Group= charlson_cut(ds$charlson_score, breaks_charlson)
        modified_count <- ifelse(ds$count> 0, 1, 0)

        s2_charlson[[counter2]] <- survfit(Surv(ds$TAR, modified_count) ~ Charlson_Group,data = ds)
        s2_charlson[[counter2]]$target      =target_id[i]
        s2_charlson[[counter2]]$outcome     =OutcomeList$outcome_id[j]
        s2_charlson[[counter2]]$databaseId  =databaseId

        ds_ir     <- cbind(cut_and_aggregate(ds),databaseId)
        write_incidence_rates_csv(ds_ir,0)

        ds_ep=mcf_episodic_prep(ae = outcome, target = target, charlson = charlson)
        mcf1_charlson=reda::mcf(Recur(id=subject_id,time=time,event=event)~Charlson_Group,data=ds_ep)
        mcf_episodic_by_charlson[[counter2]]=mcf1_charlson@MCF
        mcf_episodic_by_charlson[[counter2]]$target=target_id[i]
        mcf_episodic_by_charlson[[counter2]]$outcome=OutcomeList$outcome_id[j]
        mcf_episodic_by_charlson[[counter2]]$databaseId=databaseId
        }
    }
   }

##########################################################################################
drop_charlson_tables()
saveRDS(s1_charlson,               paste0(outputFolderIR, "/chronic_tte_by_charlson.rds"))
saveRDS(s2_charlson,               paste0(outputFolderIR, "/episodic_tte_by_charlson.rds"))
saveRDS(mcf_episodic_by_charlson,  paste0(outputFolderIR, "/mcf_episodic_by_charlson.rds"))
DatabaseConnector::disconnect(conn)
}
