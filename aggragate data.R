# load packages 
{
library(sf)
library(ggplot2)
library(dplyr)
library(tigris)
library(viridisLite)
library(zctaCrosswalk)
library(zipcodeR)
library(spatstat)
library(spdep)

library(arrow)


library(spatstat)

library(usmap)
library(scales)
#library(spatstat)

vircol <-  viridis_pal()(100)

library(duckdb)
library(plotly)
library(lubridate)
library(tidyr)
}

# aggregated data for 9 counties in NY 
# Rockland, Westchester, Bronx, New York, Richmond, Kings, Queens, Nassau, Suffolk 



ny_data <- readRDS("/Volumes/HCUPdata-CC0941-MEDSPH/Rcode_Ke/spatiotemporal_RSV/subset_data_ny.rds")


ny_data_count_zip <- ny_data %>% 
  mutate(total = ifelse(c(rsv == 1 | rsv_prim == 1), 1, 0)) %>% 
  select(zip, AGE_GROUP, AGE_GROUP2, season, date, total) %>% 
  filter(!is.na(season)) %>% 
  filter(season %in% c("2004/05", "2005/06", "2006/07", "2007/08", "2008/09", "2009/10", "2010/11", "2011/12",   "2012/13", "2013/14")) %>% 
  group_by(date, zip, AGE_GROUP2) %>% 
  summarise(N = sum(total))  %>% 
  ungroup() %>% 
  filter(!is.na(zip))





{
options(tigris_use_cache = TRUE)
# Get ZIP code shapefile for New York
ny_zips <- zctas(state = "NY", year = 2010) # Geographic information for NY


# zcta_pop <- readRDS('/Volumes/HCUPdata-CC0941-MEDSPH/Rcode_Ke/NY_pop_2019.RDS') %>%
#   select(GEOID, population = estimate) %>%
#   rename(ZCTA5CE10 = GEOID) # give population for the ZIP code level



ZIP_COUNTY <- get_zcta_crosswalk() %>% 
  filter(state_fips == "36") %>% 
  select(-state_fips) %>% 
  group_by(zcta) %>%
  #slice_head(n = 1) %>% # remove multiple entries of ZIPs
  ungroup()

zip_multiple_counties <- ZIP_COUNTY %>%
  group_by(zcta) %>%
  summarise(n_counties = n_distinct(county_name)) %>%
  filter(n_counties > 1) %>%
  pull(zcta)

ZIP_COUNTY <- ZIP_COUNTY %>%
  filter(!zcta %in% zip_multiple_counties) %>% 
  select(-county_fips)

# Get county information
zip_info <- reverse_zipcode(zip_multiple_counties)

# Extract the county name
unique_zip_county <- data.frame(zcta = zip_info$zipcode,
                                county_name = zip_info$county)

ZIP_COUNTY <- rbind(ZIP_COUNTY, unique_zip_county)


merged_zip_county_data <- merge(ZIP_COUNTY, ny_data_count_zip, 
                                by.x = "zcta", by.y = "zip")

}

selected_counties <- paste0(c("New York", "Rockland", "Westchester", "Bronx", 
                              "Richmond", "Kings", "Queens", "Nassau", "Suffolk"), " County")


case_by_county_age <- merged_zip_county_data %>% 
  group_by(county_name, date, AGE_GROUP2) %>% 
  summarise(N = sum(N)) %>% 
  ungroup() %>% 
  filter(county_name %in% selected_counties)


case_by_age <- case_by_county_age %>% 
  mutate(AGE_GROUP = ifelse(AGE_GROUP2 %in% c("20-39", "40-59"), "20-59", 
                            ifelse(AGE_GROUP2 == "under 5", "under 5", 
                                   ifelse(AGE_GROUP2 == "5-19", "5-19", 
                                          ifelse(AGE_GROUP2 == "60 above", "60 above", NA))))) %>% 
  group_by(date, AGE_GROUP) %>% 
  summarise(N = sum(N)) %>% 
  ungroup() %>% 
  mutate(AGE_GROUP = factor(AGE_GROUP, levels = c("under 5", "5-19", "20-59", "60 above"))) %>% 
  arrange(AGE_GROUP)


### construct T x G matrix for rjags 

all_dates <- seq.Date(from = as.Date("2004-07-01"), to = as.Date("2014-06-30"), by = "month")

AGE_GROUP <- unique(case_by_age$AGE_GROUP)

# Create a data frame with all combinations of dates and regions
full_data <- expand.grid(date = all_dates, AGE_GROUP = AGE_GROUP)


observation <- merge(
  full_data,
  case_by_age,
  by = c("AGE_GROUP",  "date"),
  all.x = TRUE
)  %>% 
  arrange(AGE_GROUP)


observation$N[is.na(observation$N)] <- 0


# plot data ts by county name
# observation %>%
#   ggplot() +
#   geom_line(aes(x = date, y = N, group = AGE_GROUP)) +
#   facet_wrap(~AGE_GROUP, scales = "free_y")

observation2 <- observation %>% 
  tidyr::pivot_wider(names_from = AGE_GROUP, values_from = N, values_fill = list(N = 0)) %>% 
  dplyr::select(-date) 



observed_data <- as.matrix(observation2)

#saveRDS(observed_data, "/Volumes/HCUPdata-CC0941-MEDSPH/Rcode_Ke/rjags/rjag_rsv_spatiotemporal/NYwide.rds")




########## get population offset ###########

standard_popluation <- readxl::read_xlsx(path = "/Volumes/HCUPdata-CC0941-MEDSPH/Rcode_Ke/US_2020_population_standarded.xlsx", sheet = "Data")  %>% 
  mutate(AGE = unique(sort(ny_data$AGE_GROUP))) %>% 
  rename("AGE_GROUP" = "AGE",
         "Fraction" = "Percentage",
         "Population" = "Population") %>% 
  mutate(AGE_GROUP2 = ifelse(AGE_GROUP %in% c("under 5"), "under 5", 
                             ifelse(AGE_GROUP %in% c("5-9", "10-14", "15-19"), "5-19", 
                                    ifelse(AGE_GROUP %in% c("20-24", "25-29", "30-34", "35-39"), "20-39", 
                                           ifelse(AGE_GROUP %in% c("40-44", "45-49", "50-54", "55-59"), "40-59", "60 above")))))  %>% 
  group_by(AGE_GROUP2) %>% 
  summarise(Fraction = sum(Fraction), 
            Population = sum(Population))  %>% 
  rename("AGE_GROUP" = "AGE_GROUP2",
         "Fraction" = "Fraction",
         "Population" = "Population")  %>% 
  mutate(AGE_GROUP = factor(AGE_GROUP, 
                            levels = c("under 5", "5-19", "20-39", "40-59", "60 above"))) %>% 
  arrange(AGE_GROUP)


AGE_GROUP <- sort(unique(ny_data$AGE_GROUP2))
Fraction <- standard_popluation$Fraction



