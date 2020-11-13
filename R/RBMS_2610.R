#adapted by Dafne Ram after SEBMS_indicator_example.R by Reto Schmucki

#setwd("/home/dafneram/Documents/R_butterfly_Oct2020")

#setwd("C:/DAFNE/SvenskDagfjarilsovervakning/RBMS/LUNARC")

#if(!requireNamespace("devtools")) install.packages("devtools")
#devtools::install_github("RetoSchmucki/rbms", ref = 'able_dev')

#### requires input arguments:
## [7] species number (1:28)
## [8] region (1:7)
## [9] Data source (1:4)
## [10] First year
## [11] Last year

if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("RetoSchmucki/rbms")

#load libraries
library(rbms)
library(data.table)
library(dplyr)

## load visit and count data
b_visit <- as.data.table(read.csv("data/visit_data.txt", sep = "|"))
b_count <- as.data.table(read.csv("data/count_data.txt", sep = "|"))

#list site ids with region
b_siteregion <- unique(b_count[,c('site','region')])
colnames(b_siteregion) <- c("SITE_ID", "CI_REGION") 


############## SET PARAMETERS ##############
args <- commandArgs()

#Select species number
p_Species <- as.numeric(1)  # args[7]
p_Speciess <- switch(p_Species, 8, 17, 19, 26, 28, 29, 30, 38, 40, 46, 50, 55, 67, 70, 71, 91, 92, 93, 95, 101, 105, 109, 110, 115, 117, 118, 119, 120)

## Species numbers:
#  1 = 8      8 = 38     16 = 91     23 = 110
#  2 = 17     9 = 40     17 = 92     24 = 115
#  3 = 19    10 = 46     18 = 93     25 = 117
#  4 = 26    11 = 50     19 = 95     26 = 118
#  5 = 28    12 = 55     20 = 101    27 = 119
#  6 = 29    13 = 67     21 = 105    28 = 120
#  7 = 30    14 = 70     22 = 109  

#Select region (for compiled indices, not flight curves)
p_Region  <- as.numeric(1) # args[8]
p_Regions <- switch(p_Region,
                    "SGot",                                                    # 1
                    "OGot",                                                    # 2
                    "VGotSve",                                                 # 3
                    "OSve",                                                    # 4
                    "NSveSNor",                                                # 5
                    "NNor",                                                    # 6
                    c("SGot", "OGot", "VGotSve", "OSve", "NSveSNor", "NNor"))  # 7 (whole Sweden)


#Select data sources
p_Data  <- as.numeric(1) # args[9]
p_Datas <- switch(p_Data,
                  c(54,56,63,64,67),                      #SeBMS    1
                  c(57),                                  #SLU      2
                  c(59,60,61,62,63,64,65,81),             #Lst      3
                  c(57,59,60,61,62,65,81,54,56,63,64,67)) #All data 4

b_visit <- b_visit[datasource %in% p_Datas,]
b_count <- b_count[datasource %in% p_Datas,]


#Select time period
p_FirstYear <- as.numeric(2010)   # args[10]
p_LastYear  <- as.numeric(2017)   # args[11]



############## FLIGHT CURVE & SITE INDEX ############

## generate an empty data.table to receive the site indices
sindex_bms <- data.table()

## compute flight curve for each species in Sweden

  ## subset visit and count data for that region and species
  r_visit <- b_visit
  r_count <- b_count[species == p_Speciess, ]

  ## set column names to comply with the format used by rbms functions
  setnames(r_count, c("site", "date", "species", "count"), c("SITE_ID", "DATE", "SPECIES", "COUNT"))
  setnames(r_visit, c("site", "date"), c("SITE_ID", "DATE"))
  
  ## initialize a data set that covers the time-series of interest
  ts_date <- rbms::ts_dwmy_table(InitYear = p_FirstYear, LastYear = p_LastYear, WeekDay1 = 'monday')
  
  ## define the specific monitoring season and add Anchor or "0" on both sides of the monitoring season to
  ## force the shape of the curve before and after the monitoring season and the time unit at which the data should be
  ## compiled ('w': week count; 'd': day count)
  ts_season <- rbms::ts_monit_season(ts_date, StartMonth = 4, EndMonth = 9, StartDay = 1, EndDay = NULL,
                                      CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 3, AnchorLag = 2, TimeUnit = 'w')
  
  ## add the visit data to the time-series
  ts_season_visit <- rbms::ts_monit_site(r_visit, ts_season)
  
  ## add the count data to the time-series
  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, r_count, sp = p_Speciess)
  
  ## compute the flight curve for each year, set minimum threshold to inform the GAM (minimum number of visit, minimum number of occurence,
  ## minimum number of site), and add the parameter for the functions (max of trial for convergence, error distribution family, GAM method [gam or bam])
  ## and select the time step at which the spline should be fitted ('w': week count; 'd': day count) - this must be consistent with the TimeUnit
  ## defined above
  
  ts_flight_curve <- rbms::flight_curve(ts_season_count, NbrSample = 1000, MinVisit = 3, MinOccur = 2, MinNbrSite = 3,
                                         MaxTrial = 4, GamFamily = 'nb', SpeedGam = FALSE, CompltSeason = TRUE, SelectYear = NULL,
                                         TimeUnit = 'w')
  
  ## retrieve the flight curves for the object returned above
  pheno <- ts_flight_curve$pheno
  
  ## use the relative flight curve to impute count where site have not been visited, here again, the TimeUnit must be consitent with the one
  ## define above
  impt_counts <- rbms::impute_count(ts_season_count=ts_season_count, ts_flight_curve=pheno, YearLimit= NULL, TimeUnit='w')
  
  ## estract a site index for each site that have been monitored at least 5% of the flight period observed in the region (note that this threshold
  ## is very liberal and could be set at higher level)
  sindex <- rbms::site_index(butterfly_count = impt_counts, MinFC = 0.05)
  sindex_sp <- sindex[ , FC_REGION := "Sweden"]
  
  ## append regional result in the full scheme data.table
  if(nrow(sindex_sp)>0){
    sindex_bms <- rbind(sindex_bms, sindex_sp[ , .(SPECIES, SITE_ID, M_YEAR, SINDEX, TOTAL_NM, FC_REGION)])
  }
  
  #add region to site ids for collated indices
  sindex_test <- merge(sindex_sp, b_siteregion, by = "SITE_ID")

############## BOOTSTRAP ##############

  ##select regions for compiled indices
  sindex_2 <- sindex_test[CI_REGION %in% p_Regions,]

## generate a object with n bootstrap samples (boot_n) and list of site (transect) to estimate the uncertainty
## of the collated index estimated from the statistical population of transect monitored in a given year
bootsample <- rbms::boot_sample(sindex_2, boot_n = 1000)

## compute the collated index for the original and the bootstrap samples, using the site indices in a glm of the form [sindex ~ factor(site) + factor(year)],
## corresponding to model 2 in TRIM, with a weight relative to the portion of the flight curve covered by the weekly visit recorded at each site. rm_zero is an
## argument to remove site where the species has never been observed from the model fitting, but these sites are later considered when computing the density.
## NOTE that the sindex can be standardized for a specific transect lenght (in the indicator, we do standardize to 2 km long transect which correspond to 1 ha
## monitoring area 2000 m long * 5 m wide = 10000 m2 [100x100m] )
co_index <- list()

#pb <- txtProgressBar(min = 0, max = dim(bootsample$boot_ind)[1], initial = 0, char = "*",  style = 3)


for(i in c(0, seq_len(nrow(bootsample$boot_ind)))){
  co_index[[i+1]] <- rbms::collated_index(data = sindex_2, 
                                          s_sp = p_Speciess,
                                          sindex_value = "SINDEX",
                                          bootID = i,
                                          boot_ind = bootsample,
                                          glm_weights = TRUE,
                                          rm_zero = TRUE)
  #  setTxtProgressBar(pb, i)
}

## collate and append all the result in a data.table format
co_index <- rbindlist(lapply(co_index, FUN = "[[","col_index"))
co_index[ , species := p_Speciess]

## exluded all indices below 0.0001 and above 100,000 individual per transect in a given year and compute the mean index and log transform for the observed set
## of transect monitored
co_index_b <- co_index[COL_INDEX > 0.0001 & COL_INDEX < 100000, ]
co_index_logInd <- co_index_b[BOOTi == 0, .(M_YEAR, COL_INDEX)][, log(COL_INDEX)/log(10), by = M_YEAR][, mean_logInd := mean(V1)]

## merge the mean log index with the full bootstrap dataset
setnames(co_index_logInd, "V1", "logInd"); setkey(co_index_logInd, M_YEAR); setkey(co_index_b, M_YEAR)
co_index_b <- merge(co_index_b, co_index_logInd, all.x = TRUE)

## compute the log collated index for each bootstap samples
setkey(co_index_b, BOOTi, M_YEAR)
co_index_b[ , boot_logInd := log(COL_INDEX)/log(10)]

## compute the metric used for the graph of the Collated Log-Index centered around 2 (observed, bootstap sampel, credible confidence interval, linear trend)
b1 <- data.table(M_YEAR = co_index_b$M_YEAR, LCI = 2 + co_index_b$boot_logInd - co_index_b$mean_logInd)
b2 <- data.table(M_YEAR = co_index_b[BOOTi == 0, M_YEAR], LCI= 2 + co_index_b[BOOTi == 0, logInd] - co_index_b[BOOTi == 0, mean_logInd])
b5 <- b1[co_index_b$BOOTi != 0, quantile(LCI, 0.025, na.rm = TRUE), by = M_YEAR]
b6 <- b1[co_index_b$BOOTi != 0, quantile(LCI, 0.975, na.rm = TRUE), by = M_YEAR]

## output data
LCI_1 <- merge(b2, b5,  by = "M_YEAR") 
LCI_out <- merge(LCI_1, b6, by = "M_YEAR")
LCI_out$SPECIES <- p_Speciess
colnames(LCI_out) <- c("YEAR", "LCI", "Q.025", "Q.975", "SPECIES")

write.csv2(LCI_out, paste0("S",p_Speciess,"_", "R", p_Region,"_", "D", p_Data, "_", p_FirstYear, "_", p_LastYear, ".csv"), row.names = FALSE)

