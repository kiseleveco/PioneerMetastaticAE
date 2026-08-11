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

#' Execute the cohort diagnostics
#'
#' @details
#' This function extracts the data needed to calculate incidence rates of outcomes across targets.
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
#'                                            will need to have write privileges in this schema. Note
#'                                            that for SQL Server, this should include both the
#'                                            database and schema name, for example 'cdm_data.dbo'.
#' @param cohortTable                         The name of the table that will be created in the work
#'                                            database schema. This table will hold the exposure and
#'                                            outcome cohorts used in this study.
#' @param tempEmulationSchema                 Some database platforms like Oracle and Impala do not
#'                                            truly support temp tables. To emulate temp tables,
#'                                            provide a schema with write privileges where temp tables
#'                                            can be created.
#' @param verifyDependencies                  Check whether correct package versions are installed?
#' @param outputFolderIR          	      Name of local folder to place results; make sure to use
#'                                            forward slashes (/). Do not use a folder on a network
#'                                            drive since this greatly impacts performance.
#' @param databaseId                          A short string for identifying the database (e.g.
#'                                            'Synpuf').
#' @param databaseName                        The full name of the database (e.g. 'Medicare Claims
#'                                            Synthetic Public Use Files (SynPUFs)').
execute_IR <- function(connectionDetails,
                       cdmDatabaseSchema,
                       cohortDatabaseSchema = cdmDatabaseSchema,
                       cohortTable = "cohort",
                       tempEmulationSchema = getOption("sqlRenderTempEmulationSchema"),
                       breaks_age = c(18, 55, 75, 120),
                       databaseId = databaseId,
					             outputFolderIR=outputFolderIR,
					             Grouping=TRUE) 
{
#check if output folder is present and if not create it 
if (!file.exists(outputFolderIR)) {
  dir.create(outputFolderIR, recursive = TRUE)
}
# the csv where the analysis settings are
IRsettings = read.csv("./inst/settings/IRsettings.csv")
 #Connection details
conn=DatabaseConnector::connect(connectionDetails)
##########################################################################################

# SQL scripts





###############################################
###############################################
###############################################
sql_all_outcomes='SELECT T.COHORT_DEFINITION_ID as target_id, O.COHORT_DEFINITION_ID as outcome_id, T.SUBJECT_ID ,T.COHORT_START_DATE,O.COHORT_START_DATE as OUTCOME_START_DATE,T.COHORT_END_DATE, p.YEAR_OF_BIRTH
FROM @cohortDatabaseSchema.@cohortTable T
JOIN @cohortDatabaseSchema.@cohortTable O ON T.subject_id = O.subject_id
JOIN @cdmDatabaseSchema.person p ON T.subject_id = p.person_id
WHERE T.cohort_definition_id = @target_id
AND YEAR(T.COHORT_START_DATE) - p.YEAR_OF_BIRTH > 18
  AND O.cohort_definition_id = @outcome_id
  AND O.cohort_start_date >= T.cohort_start_date
  AND O.cohort_end_date <= T.cohort_end_date'
###############################################
###############################################
###############################################
#first outcome during after target
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
###############################################
###############################################
###############################################
# All the target population 
sql_target ='
SELECT T.COHORT_DEFINITION_ID as target_id, T.SUBJECT_ID, T.COHORT_START_DATE, T.COHORT_END_DATE, p.YEAR_OF_BIRTH
FROM @cohortDatabaseSchema.@cohortTable T 
JOIN @cdmDatabaseSchema.person p ON T.subject_id = p.person_id 
WHERE T.COHORT_DEFINITION_ID=@target_id AND YEAR(T.COHORT_START_DATE) - p.YEAR_OF_BIRTH > 18'
###############################################
# For the chronic population patients having the event befor the cohort start date are excluded
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
###############################################
###############################################
###############################################

sql_grouping=' 
WITH new_cohort AS (
  SELECT 
    SUBJECT_ID,
    MIN(COHORT_START_DATE) as min_start_date,
    MAX(COHORT_END_DATE) as max_end_date
  FROM 
    @cohortDatabaseSchema.@cohortTable
  WHERE 
    COHORT_DEFINITION_ID IN (@cohort_list) -- replace with your cohort IDs
  GROUP BY 
    SUBJECT_ID
) 
INSERT INTO @cohortDatabaseSchema.@cohortTable (COHORT_DEFINITION_ID, SUBJECT_ID, COHORT_START_DATE, COHORT_END_DATE)
SELECT 
  @new_cohort_id, -- replace with your new cohort ID
  SUBJECT_ID,
  min_start_date,
  max_end_date
FROM 
  new_cohort;'
###############################################
###############################################
###############################################
sql_drop_cohort_id= "DELETE FROM @cohortDatabaseSchema.@cohortTable 
WHERE COHORT_DEFINITION_ID = @cohort_id;"



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
# Utility functions



#############################################################
#############################################################
#############################################################
# sql functions
# Translate the adverse events sql to the target dialect 
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

#############################################################
#############################################################
#############################################################
# Translate the target population sql to the target dialect
sql_translate_target=function(sql,target_id)
{
sql_rendered= SqlRender::render(sql,
				  cdmDatabaseSchema=cdmDatabaseSchema,
				  cohortDatabaseSchema=cohortDatabaseSchema,
				  cohortTable=cohortTable,
				  target_id=target_id)
return(SqlRender::translate(sql_rendered,targetDialect = DBMS))
}


#############################################################
#############################################################
#############################################################
# create treatment group sql 
sql_execute_grouping=function(sql,cohortDatabaseSchema,cohortTable,cohort_list,new_cohort_id)
{
sql_rendered= SqlRender::render(sql,
          cohortDatabaseSchema=cohortDatabaseSchema,
          cohortTable=cohortTable,
          cohort_list=cohort_list,
          new_cohort_id=new_cohort_id)
 DatabaseConnector::executeSql(conn, SqlRender::translate(sql_rendered,targetDialect = DBMS))
}
 #############################################################
# sql for dropping cohort id from the cohort table
sql_execute_drop_cohort_id=function(cohortDatabaseSchema,cohortTable,cohort_id)
{
sql_rendered=SqlRender::render(sql = sql_drop_cohort_id,
                              cohortDatabaseSchema=cohortDatabaseSchema,
                              cohortTable=cohortTable,
                              cohort_id=cohort_id)
  DatabaseConnector::executeSql(conn, SqlRender::translate(sql_rendered,targetDialect = DBMS))
}
#############################################################
#############################################################
#############################################################

# execute grouping if true
 
#############################################################
#############################################################
#############################################################
get_data=function(sql)
{
suppressMessages(as.data.frame(DatabaseConnector::querySql(conn,sql)))
}


#############################################################
#############################################################
#############################################################
# data preparation function for the chronic adverse events
data_prep_chronic <- function(ae, target){
    
# create a table for the adverse events
  df_ae <- mutate(ae,
                   target_id              = ae$TARGET_ID,
                   outcome_id             = ae$OUTCOME_ID,
                   subject_id             = ae$SUBJECT_ID,
                   age_at_index_T         = lubridate::year(ae$COHORT_START_DATE) - ae$YEAR_OF_BIRTH,
                   time_to_outcome_years  = as.numeric(difftime(ae$OUTCOME_START_DATE, ae$COHORT_START_DATE, units = "days"))/365.25,
                   index_year_target      = lubridate::year(ae$COHORT_START_DATE)) %>%
              group_by(target_id, outcome_id, subject_id, age_at_index_T, time_to_outcome_years, index_year_target) %>%
              summarise(count = n())
  
# Create the dataframe all patients at risk (i.e. all patients in the target population)
  df_target <-  mutate(target,
                         target_id         = target$TARGET_ID,
                         subject_id        = target$SUBJECT_ID,
                         age_at_index_T    = lubridate::year(target$COHORT_START_DATE) - target$YEAR_OF_BIRTH,
                         follow_up_years   = as.numeric(difftime( target$COHORT_END_DATE, target$COHORT_START_DATE, units = "days"))/365.25,
                         index_year_target = lubridate::year(target$COHORT_START_DATE)) %>% 
                      group_by(target_id, subject_id, age_at_index_T, index_year_target) %>%
                  summarise(follow_up_time = sum(follow_up_years))
  
  # Combine the dataframes with information on adverse events and all patients
  out_id  <- unique(df_ae$outcome_id)
  suppressMessages(
  df_comb <- (left_join(df_target, df_ae) %>% 
                mutate(TAR = if_else(!is.na(time_to_outcome_years), time_to_outcome_years, follow_up_time),
		                   outcome_id = out_id) %>%
                mutate(count = if_else(is.na(count), 0, count))))
  
  return(as.data.frame(df_comb))
}
#############################################################
#############################################################
# data preparation function for the episodic adverse events
# The only diffrence bettween the two functions is that the chronic function
# does not include the follow up time after the first event occurs
data_prep_episodic <- function(ae, target){
    
# create a table for the adverse events
  df_ae <- mutate(ae,
                   target_id              = ae$TARGET_ID,
                   outcome_id             = ae$OUTCOME_ID,
                   subject_id             = ae$SUBJECT_ID,
                   age_at_index_T         = lubridate::year(ae$COHORT_START_DATE) - ae$YEAR_OF_BIRTH,
                   time_to_outcome_years  = as.numeric(difftime(ae$OUTCOME_START_DATE, ae$COHORT_START_DATE, units = "days"))/365.25,
                   index_year_target      = lubridate::year(ae$COHORT_START_DATE)) %>%
              group_by(target_id, outcome_id, subject_id, age_at_index_T,  index_year_target) %>%
              summarise(count = n(),
                        time_to_outcome_years = min(time_to_outcome_years))
# Create the dataframe all patients at risk (i.e. all patients in the target population)
  df_target <-  mutate(target,
                         target_id         = target$TARGET_ID,
                         subject_id        = target$SUBJECT_ID,
                         age_at_index_T    = lubridate::year(target$COHORT_START_DATE) - target$YEAR_OF_BIRTH,
                         follow_up_years   = as.numeric(difftime( target$COHORT_END_DATE, target$COHORT_START_DATE, units = "days"))/365.25,
                         index_year_target = lubridate::year(target$COHORT_START_DATE)) %>% 
                      group_by(target_id, subject_id, age_at_index_T, index_year_target) %>%
                  summarise(follow_up_time = sum(follow_up_years))
  
  # Combine the dataframes with information on adverse events and all patients
  out_id  <- unique(df_ae$outcome_id)
  suppressMessages(
  df_comb <- (left_join(df_target, df_ae) %>% 
                mutate(TAR =  follow_up_time,
                                 outcome_id = out_id) %>%
                mutate(count = if_else(is.na(count), 0, count))))
  
  return(as.data.frame(df_comb))
}
#############################################################
#############################################################
#############################################################
# cut and aggregate the data at the desired age groups
cut_and_aggregate = function(df)
{
  break_points_age  <- c(18, breaks_age, Inf)
  df                <- mutate(df, age_cat = cut(age_at_index_T, breaks = break_points_age), include.lowest = TRUE)
  out_id            <- unique(df$outcome_id)
  suppressMessages(
  df <- (group_by(df, target_id, age_cat, index_year_target) %>%
           summarise(nr_events = sum(count, na.rm = TRUE), time_at_risk = sum(TAR), nr_patients = n())))
  
  # Include possibly missing groups
  res        <- expand.grid(unique(df$target_id), 
                            unique(df$age_cat), 
                            unique(df$index_year_target))
  names(res) <- c("target_id", "age_cat", "index_year_target")  
    suppressMessages(
  res        <- (left_join(res, df) %>% 
                   mutate(nr_events    = coalesce(nr_events, rep(0, nrow(res))), 
                          time_at_risk = coalesce(time_at_risk, rep(0, nrow(res))),
                          nr_patients  = coalesce(nr_patients, rep(0, nrow(res))))))
  res$outcome_id <- out_id
  return(as.data.frame(res))
}
#############################################################
# Episodic data preparation conforming to reda mcf function
mcf_episodic_prep=function(ae,target)
{
  df_ae <- mutate(ae,
                   target_id              = ae$TARGET_ID,
                   outcome_id             = ae$OUTCOME_ID,
                   subject_id             = ae$SUBJECT_ID,
                   age_at_index_T         = lubridate::year(ae$COHORT_START_DATE) - ae$YEAR_OF_BIRTH,
                   time_to_outcome_years  = as.numeric(difftime(ae$OUTCOME_START_DATE, ae$COHORT_START_DATE, units = "days"))/365.25,
                   index_year_target      = lubridate::year(ae$COHORT_START_DATE)) %>%
              group_by(target_id, outcome_id, subject_id, age_at_index_T, time_to_outcome_years, index_year_target) %>%
              summarise(count = n())
  
# Create the dataframe all patients at risk (i.e. all patients in the target population)
  df_target <-  mutate(target,
                         target_id         = target$TARGET_ID,
                         subject_id        = target$SUBJECT_ID,
                         age_at_index_T    = lubridate::year(target$COHORT_START_DATE) - target$YEAR_OF_BIRTH,
                         follow_up_years   = as.numeric(difftime( target$COHORT_END_DATE, target$COHORT_START_DATE, units = "days"))/365.25,
                         index_year_target = lubridate::year(target$COHORT_START_DATE)) %>% 
                    	 group_by(target_id, subject_id, age_at_index_T, index_year_target) %>%
                  summarise(follow_up_time = sum(follow_up_years))

 df_ae_1=data.frame(cbind(	                     df_ae$target_id,
							df_ae$outcome_id,
							df_ae$subject_id,
							df_ae$age_at_index_T,
							df_ae$index_year_target,
							df_ae$time_to_outcome_years,
                                                 event=1))		
names(df_ae_1)=c("target_id","outcome_id","subject_id","Age_Group","Year","time","event")
 
outcome_id=unique(df_ae_1$outcome_id)
 df_o_1=data.frame(cbind(                 df_target$target_id ,
						outcome_id,
						df_target$subject_id ,
						df_target$age_at_index_T,
						df_target$index_year_target,
						df_target$follow_up_time,
                                          event=0 ))
names(df_o_1)=c("target_id","outcome_id","subject_id","Age_Group","Year","time","event")

df=data.frame(rbind(df_o_1, df_ae_1))
df$Age_Group= cut(df$Age_Group, breaks = breaks_age)
return(df)
}
#############################################################
#############################################################
#############################################################
# Filter IR data 
get_IR_overall=function(df)
{
		 res2 <- df %>%
		  group_by(target_id,outcome_id  )%>%
		  summarize(
			total_events		=     sum(nr_events),
			total_time_at_risk 	= sum(time_at_risk),
			total_nr_patients 	= sum(nr_patients))%>%
	  	mutate(incidence_rate = (total_events / total_time_at_risk) * 1000)
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
#############################################################
get_IR_by_year_only=function(df)
{
  res2 <- df %>%
		group_by(target_id,outcome_id,index_year_target )%>%
		summarize(
		total_events		    = sum(nr_events),
		total_time_at_risk 	= sum(time_at_risk),
		total_nr_patients 	= sum(nr_patients)) %>%
		mutate(incidence_rate = (total_events / total_time_at_risk) * 1000)
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

#############################################################
get_IR_by_age_only=function(df)
{
res2 <- df %>%
		group_by(target_id,outcome_id,age_cat )%>%
		summarize(
		total_events		= sum(nr_events),
		total_time_at_risk 	= sum(time_at_risk),
		total_nr_patients 	= sum(nr_patients))%>%
		mutate(incidence_rate = (total_events / total_time_at_risk) * 1000)
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
#############################################################
get_IR_by_year_and_age=function(df)
{
res2 <- df %>%
		group_by(target_id, outcome_id, age_cat, index_year_target)%>%
		summarize(
		total_events		    = sum(nr_events),
		total_time_at_risk 	= sum(time_at_risk),
		total_nr_patients 	= sum(nr_patients))%>%
		mutate(incidence_rate = (total_events / total_time_at_risk) * 1000)
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

#############################################################
#############################################################
#############################################################
# function to write the data to a csv file
write_table=function(object,csv_file)
{
if (file.exists(csv_file)) {
		  write.table(object, csv_file, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
		} else 	{
		  write.table(object, csv_file, sep = ",", col.names = TRUE, row.names = FALSE)
				}     
}
#############################################################
#############################################################
#############################################################

# Censor any count less than 5 in the IR data
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

# Censor counts in the counts data frame
censor_counts_df=function(df)
{
	df$target_count[df$target_count<5 & df$target_count!=0]="<5"
	df$outcome_count[df$outcome_count<5 & df$outcome_count!=0]="<5"
	return(df)
}

	
 
#############################################################
#############################################################
#############################################################
# function to write the data to a csv file 
write_incidence_rates_csv=function(df,chronic_indicator)
{
IR_by_age_only=get_IR_by_age_only(df)
IR_by_age_only$databaseId=databaseId
IR_by_age_only=censor_counts_IR(IR_by_age_only)

IR_by_year_only=get_IR_by_year_only(df)
IR_by_year_only$databaseId=databaseId
IR_by_year_only=censor_counts_IR(IR_by_year_only)

IR_by_year_and_age=get_IR_by_year_and_age(df)
IR_by_year_and_age$databaseId=databaseId
IR_by_year_and_age=censor_counts_IR(IR_by_year_and_age)

IR_overall=get_IR_overall(df)
IR_overall$databaseId=databaseId
IR_overall=censor_counts_IR(IR_overall)

    if(chronic_indicator==1)
    {
    csv_file_1 <- paste0(outputFolderIR,"/chronic_IR_by_age_only.csv")
    write_table(IR_by_age_only,csv_file_1)
    csv_file_2 <- paste0(outputFolderIR,"/chronic_IR_by_year_only.csv")
    write_table(IR_by_year_only,csv_file_2)
    csv_file_3 <- paste0(outputFolderIR,"/chronic_IR_by_year_and_age.csv")
    write_table(IR_by_year_and_age,csv_file_3)
    csv_file_4 <- paste0(outputFolderIR,"/chronic_IR_overall.csv")
    write_table(IR_overall,csv_file_4)
    }
    if (chronic_indicator==0)
    {
    csv_file_1 <- paste0(outputFolderIR,"/episodic_IR_by_age_only.csv")
    write_table(IR_by_age_only,csv_file_1)
    csv_file_2 <- paste0(outputFolderIR,"/episodic_IR_by_year_only.csv")
    write_table(IR_by_year_only,csv_file_2)
    csv_file_3 <- paste0(outputFolderIR,"/episodic_IR_by_year_and_age.csv")
    write_table(IR_by_year_and_age,csv_file_3)
    csv_file_4 <- paste0(outputFolderIR,"/episodic_IR_overall.csv")
    write_table(IR_overall,csv_file_4)
    }
}

#############################################################
#############################################################
#############################################################

create_grouping_cohorts=function()
{
groups=read.csv(file="./inst/settings/groups.csv")
for(i in 1:nrow(groups))
{
groups$cohort_list[i]=gsub(";",",",groups$cohort_list[i]) # replace ; with , in the cohort list
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
#############################################################
#############################################################
#############################################################
if(Grouping==TRUE)
{
create_grouping_cohorts()
}
# execute the IR analysis
target_id  <- unique(IRsettings$target_id)
outcome_id <- unique(IRsettings$outcome_id)
counter1=0   #counter for chronic  events
counter2=0   #counter for episodic events
s1=NULL
s1_age=NULL
s2=NULL
s2_age=NULL
mcf1=NULL
mcf1_age=NULL
mcf_episodic=NULL
mcf_episodic_by_age=NULL
df_counts=  cbind(IRsettings,target_count=0,outcome_count=0)
# starting with the target cohorts
for(i in 1:length(target_id)){
    target <- sql_translate_target(sql = sql_target, target_id = target_id[i]) %>% get_data
	if(nrow(target)<5)
	{
	print(paste("No data for cohort:", target_id[i]))	 
	next
	}
# Select the outcomes for this target cohorts and start the loop
	OutcomeList=IRsettings[IRsettings$target_id==target_id[i],]
## Loop over the outcomes
	for(j in 1:nrow(OutcomeList)){
		if(OutcomeList$chronic[j]==1)
		{
		# target_chronic is needed to exclude those experiencing the adverse event before the 
			target_chronic= sql_translate_ae(sql = sql_target_chronic,target_id=target_id[i], outcome_id = OutcomeList$outcome_id[j]) %>% get_data
		     if(nrow(target_chronic)<5)
			{
				print(paste("Not enough data for target cohort:", target_id[i],"and outcome:",OutcomeList$outcome_id[j]))
 				next
			} 
		# add to the counts df
		indicator=which(		OutcomeList	$target_id		[j]	==	df_counts	$target_id		&
						        OutcomeList	$outcome_id		[j]	==	df_counts	$outcome_id		&
					        	OutcomeList	$episodic	  	[j]	==	df_counts	$episodic		&
					        	OutcomeList	$chronic	  	[j]	==	df_counts	$chronic)
		df_counts$target_count[indicator]=nrow(target_chronic)

        outcome <- sql_translate_ae(sql = sql_first_outcome,target_id=target_id[i], outcome_id = OutcomeList$outcome_id[j]) %>% get_data
       if(nrow(outcome)<5)
			{
				print(paste("Not enough outcomes for target cohort:", target_id[i],"and outcome:",OutcomeList$outcome_id[j]))
 				next
			}
	  df_counts$outcome_count[indicator]=nrow(outcome)
		counter1=counter1+1
		print(paste( " Calculating incidence for target:", target_id[i],"and outcome:", OutcomeList$outcome_id[j]))
		# Get the prepped data set
		ds <- data_prep_chronic(ae = outcome, target = target_chronic)
    Age_Group= cut(ds$age_at_index_T, breaks = breaks_age)

		# fit Keplan-Meier Curves (overall (s1) and by age s1_age) - chronic events
		 
		s1[[counter1]] <- survfit(Surv(ds$TAR, ds$count) ~  1,   data = ds)
		s1[[counter1]]$target		  =target_id[i]
		s1[[counter1]]$outcome		=OutcomeList$outcome_id[j]	
		s1[[counter1]]$databaseId	=databaseId

		s1_age[[counter1]] <- survfit(Surv(ds$TAR, ds$count) ~  Age_Group,   data = ds)
		s1_age[[counter1]]$target		  =target_id[i]
		s1_age[[counter1]]$outcome		=OutcomeList$outcome_id[j]	
		s1_age[[counter1]]$databaseId	=databaseId
		ds_ir     <- cbind(cut_and_aggregate(ds, breaks_age),databaseId)
    write_incidence_rates_csv(ds_ir,1)

		}
  
########################################################################
########################################################################
 # episodic events
		 if(OutcomeList$episodic[j]==1)
		{
           outcome <- sql_translate_ae(sql = sql_all_outcomes,target_id=target_id[i], outcome_id = OutcomeList$outcome_id[j]) %>% get_data
		     if(nrow(outcome)<5)
			   {
			   print(paste("Not enough outcomes for cohort:", target_id[i],"and outcome:",OutcomeList$outcome_id[j]))
			   next
			   }
			indicator=which(OutcomeList		$target_id		[j]	==	df_counts	$target_id		&
						        OutcomeList	$outcome_id		[j]	==	df_counts	$outcome_id		&
					        	OutcomeList	$episodic	  	[j]	==	df_counts	$episodic		&
					        	OutcomeList	$chronic	  	[j]	==	df_counts	$chronic)

		df_counts$outcome_count[indicator]=nrow(outcome)
		df_counts$target_count[indicator]=nrow(target)
 		counter2=counter2+1
	    print(paste( " Calculating incidence for target:", target_id[i],"and outcome:", OutcomeList$outcome_id[j]))
		ds <- data_prep_episodic(ae = outcome, target = target)
		Age_Group= cut(ds$age_at_index_T, breaks = breaks_age)
		modified_count <- ifelse(ds$count> 0, 1, 0)
	 
		s2[[counter2]] <- survfit(Surv(ds$TAR, modified_count) ~ 1,data = ds)
		s2[[counter2]]$target		=target_id[i]
		s2[[counter2]]$outcome		=OutcomeList$outcome_id[j]	
		s2[[counter2]]$databaseId	=databaseId
		
		s2_age[[counter2]] <- survfit(Surv(ds$TAR, modified_count) ~ Age_Group,data = ds)
		s2_age[[counter2]]$target		=target_id[i]
		s2_age[[counter2]]$outcome		=OutcomeList$outcome_id[j]	
		s2_age[[counter2]]$databaseId	=databaseId
  
 		ds_ir     <- cbind(cut_and_aggregate(ds, breaks_age),databaseId)
		write_incidence_rates_csv(ds_ir,0)
		
		ds_ep=mcf_episodic_prep(ae = outcome, target = target)
			
 		mcf1=reda::mcf(Recur(id=subject_id,time=time,event=event)~1,data=ds_ep)
		mcf1_age=reda::mcf(Recur(id=subject_id,time=time,event=event)~Age_Group,data=ds_ep)
		mcf_episodic[[counter2]]= mcf1@MCF
		mcf_episodic[[counter2]]$target=target_id[i]
		mcf_episodic[[counter2]]$outcome=OutcomeList$outcome_id[j]
		mcf_episodic[[counter2]]$databaseId=databaseId
		
		mcf_episodic_by_age[[counter2]]=mcf1_age@MCF
 		mcf_episodic_by_age[[counter2]]$target=target_id[i]
		mcf_episodic_by_age[[counter2]]$outcome=OutcomeList$outcome_id[j]
		mcf_episodic_by_age[[counter2]]$databaseId=databaseId
		}
	}

   }
##########################################################################################
   df_counts=cbind(df_counts,databaseId)
   df_counts=censor_counts_df(df_counts)
   saveRDS(s1,				        paste0(outputFolderIR,	 	"/chronic_tte.rds"))
   saveRDS(s1_age,			      	paste0(outputFolderIR,      "/chronic_tte_by_age.rds"))
   saveRDS(s2,				       	paste0(outputFolderIR, 		"/episodic_tte.rds"))
   saveRDS(s2_age,			       	paste0(outputFolderIR, 		"/episodic_tte_by_age.rds"))
   saveRDS(mcf_episodic,		    paste0(outputFolderIR,		"/mcf_episodic.rds"))
   saveRDS(mcf_episodic_by_age,		paste0(outputFolderIR, 	    "/mcf_episodic_by_age.rds"))
   write.csv(df_counts,			    paste0(outputFolderIR, 		"/counts.csv"))
}
	
  
