source("R-code/00-preamble.R")

# load data ---------------------------------------------------------------
# this red list has been synthesized and harmonized by Laura Méndez Cuéllar
# explore the data here
drl <- read_excel("Data/rl_combined_REDLIST_updated.xlsx")

# clean and explore the limitations and caveats of the data
# some ISO3_codes are NA
drl %>% filter(is.na(ISO3_Code)) %>% select(Country) %>% distinct

# only for Ethiopia and Eritrea it would make sense to correct, because there
# is clear attribution to iso3 regions
# find the rows where the country is "Ethiopia and Eritrea"
ethiopia_eritrea_row <- drl[drl$Country == "Ethiopia and Eritrea", ]

# create two new rows for "Ethiopia" and "Eritrea"
ethiopia_row <- ethiopia_eritrea_row
ethiopia_row$Country <- "Ethiopia"
ethiopia_row$ISO3_Code <- "ETH"

eritrea_row <- ethiopia_eritrea_row
eritrea_row$Country <- "Eritrea"
eritrea_row$ISO3_Code <- "ERI"

# remove the old "Ethiopia and Eritrea" row
drl <- drl[drl$Country != "Ethiopia and Eritrea", ]

# bind the new rows back into the dataframe
drl <- rbind(drl, ethiopia_row, eritrea_row)

# remove other NA ISOs
drl <- drl %>% filter(!is.na(ISO3_Code))



# how many are near threatened, include? those are a lot of species in fact.
# probably good to include them as well, as they indicate decline (at least for
# an si analysis)
drl %>% filter(Redlist_cat == "NT") %>% nrow
drl %>% filter(Redlist_cat == "NT") %>% select(accepted_taxon_name) %>% distinct %>%  nrow
# drl %>% filter(redlist_source == "rl_europe") %>% View

# indicate threat status appropriately
drl <- drl %>%
  mutate(Threat = ifelse(Redlist_cat == "NT", "near_threatened", Threat)) %>%
  mutate(Threat = ifelse(Redlist_cat == "V", "near_threatened", Threat)) %>%
  mutate(Threat = ifelse(Threat == "extinct", "threatened", Threat)) 

# check, looks good
drl %>% group_by(Threat) %>% count



# inspect the data further
# how many species per list
drl %>% group_by(Country, ISO3_Code, region_name, ISO2_Code, Publ_year, redlist_source) %>% 
  summarise(nspec = n_distinct(accepted_taxon_name)) %>% View()

# how many subnational lists
drl %>% select(Country, ISO3_Code, region_name, ISO2_Code, Publ_year, redlist_source) %>% 
  distinct %>% 
  filter(!is.na(region_name)) %>% View()

# create a subnational data frame, this makes sense for the very large countries
# for india unfortunately, it is not really possible to assign these regions
# (many of which are duplicates in some way) to iso2 regions...
subnational_drl <- drl %>% 
  select(Country, ISO3_Code, region_name, ISO2_Code, Publ_year, redlist_source) %>% 
  distinct %>% 
  filter(!is.na(region_name)) %>% 
  filter(Country == "United States" |
           Country == "Russia"|
           Country == "Canada")

# how many subnational rls
subnational_drl %>% 
  select(Country, region_name) %>% distinct %>% 
  group_by(Country) %>% count()
# this is pretty complete, I checked with the respective country territories
# to find out how complete these are. USA has 51 lists because washington dc
# has an extra list.

# filter these countries then out for the national rls
subnational_drl <- left_join(subnational_drl %>% select(Country, ISO2_Code), drl)
drl <- anti_join(drl, subnational_drl %>% select(Country))

# correct ISO2 Codes for subnational
# add ISO3 Code in front of the ISO2 for USA and Canada, and change Newfoundland
# NF to NL
subnational_drl <- subnational_drl %>% 
  mutate(ISO2_Code = ifelse(ISO3_Code == "USA", 
                            paste0("USA-", ISO2_Code), ISO2_Code)) %>% 
  mutate(ISO2_Code = ifelse(ISO3_Code == "CAN", 
                            paste0("CAN-", ISO2_Code), ISO2_Code)) %>%
  mutate(ISO2_Code = ifelse(ISO2_Code == "CAN-NF", 
                            "CAN-NL", ISO2_Code))


# further inspection: 
# shippmann records only extinct (RE, EW, EX, likely extinct) species
# doesn't present complete RLs.
# drl %>% filter(redlist_source == "rl_shippmann") %>% View

# Argentinina only extinct species (does not really mean we have RL data for this
# country), but we keep it in nonetheless as these data are better than
# nothing

# check the regional lists again, there inclusion would be a bit
# overkill if a national assessment is present
drl %>% select(Country, ISO3_Code, region_name, ISO2_Code, Publ_year, redlist_source) %>% 
  distinct %>% 
  filter(!is.na(region_name)) %>% View()

# ok, are there normal ones for these countries actually, or are these the only ones?
drl %>% group_by(Country, ISO3_Code, region_name, ISO2_Code, Publ_year, redlist_source) %>% 
  summarise(nspec = n_distinct(accepted_taxon_name)) %>% 
  inner_join(drl %>% select(Country, ISO3_Code, region_name, ISO2_Code, Publ_year, redlist_source) %>% 
               distinct %>% 
               filter(!is.na(region_name)) %>% 
               select(Country) %>% distinct) %>% View

# for germany, malaysia, australia and spain throw out (for UK/India there is no entire national one),
# this makes the rest of these red lists national RLs.
drl <- drl %>% 
  filter(is.na(region_name) | 
           !region_name %in% c("Victoria", "Baden-Wuerttemberg", "Berlin", 
                               "Nordrhein-Westfalen", "Andalucia", "Galicia", 
                               "Valencia", "Peninsular Malaysia", "Sabah", 
                               "Sabah and Sarawak", "Sarawak"))

# figure for publication years of RLs
drl %>% group_by(Country, ISO3_Code, region_name, ISO2_Code) %>% 
  filter(!is.na(Publ_year)) %>% 
  filter(Publ_year != 0) %>% 
  summarise(ntimes = n_distinct(Publ_year), years = Publ_year) %>% 
  ungroup() %>% 
  select(Country, ntimes, years) %>% distinct %>% 
  ggplot(aes(x = years, y = reorder(Country, years), group = Country)) +
  geom_point(size = 2, color = "#FAC55F") +
  geom_line(size = 0.7, color = "#FAC55F") +
  labs(x = "Year", y = "") +
  theme_minimal(base_family = "Arial Narrow", base_size = 13) +
  theme(panel.grid = element_blank())


showtext_opts(dpi=600)
ggsave("Figures/pub_years.png", dpi = 600, bg = "white", height = 13, width = 12)
showtext_opts(dpi=96)


# how many countries have lists from more than one year
drl %>% group_by(Country, ISO3_Code, region_name, ISO2_Code) %>% 
  filter(!is.na(Publ_year)) %>% 
  filter(Publ_year != 0) %>%
  summarise(ntimes = n_distinct(Publ_year), years = Publ_year) %>% 
  filter(ntimes > 1) %>%
  ungroup() %>% 
  select(Country) %>% distinct

# check how often it happens that publication years are in a 5 year window,
# this would make it unlikely that these are entire updated assessments
drl %>% group_by(Country, ISO3_Code, region_name, ISO2_Code) %>% 
  filter(!is.na(Publ_year)) %>% 
  mutate(ntimes = n_distinct(Publ_year, na.rm = TRUE)) %>%
  filter(ntimes > 1) %>%
  ungroup() %>% 
  select(Country, ntimes, Publ_year) %>% distinct %>%
  group_by(Country) %>%
  filter(!is.na(Publ_year)) %>%
  mutate(Publ_year = as.numeric(Publ_year)) %>% 
  filter(within_5_years(Publ_year)) %>%
  ungroup() %>% 
  select(Country) %>% distinct %>% View()


# okay, now it would be important to exclude downlisted species.
# that means if I have a species that was listed as threatened in 2002, and
# not anymore in 2009 this species should be dropped
downlisted_species <- drl %>% 
  group_by(Country, ISO3_Code, region_name, ISO2_Code, accepted_taxon_name) %>% 
  arrange(Publ_year) %>% # Arrange by publication year within each group
  mutate(
    # Get the previous publication year
    previous_year = lag(Publ_year),
    # Check if the previous status was "threatened", the current is not, and the years are different
    previous_threat_status = lag(Threat),
    downlisted = if_else(previous_threat_status == "threatened" & 
                           Threat == "near_threatened" & Publ_year != previous_year, 
                         TRUE, FALSE, missing = FALSE)
  ) %>%
  ungroup() %>% 
  filter(downlisted == TRUE)
downlisted_species %>% count(downlisted)
downlisted_species %>% View()
# okay let's get rid of these instances, but actually, not of the downlisted
# ones but of that species generally for that country/region

# check also how many rows will be removed
nrow(drl)
drl <- drl %>% anti_join(downlisted_species %>% 
                           select(accepted_taxon_name, ISO3_Code) %>% 
                           distinct)
nrow(drl)


# check for how many countries we have lists from multiple sources
drl %>% select(Country, redlist_source) %>% distinct %>% 
  group_by(Country) %>% 
  count() %>% 
  filter(n > 1) %>% 
  nrow()


# sometimes a country can have for the same RL year a conflicting threat status
# based on multiple sources, be conservative and take the lower status,
# wherever this is the case
# reduce the dimensionality of this data frame first
drl <- drl %>% select(accepted_taxon_name,
                      accepted_taxon_authors,
                      status = Threat,
                      Redlist_cat,
                      ISO3_Code,
                      ISO2_Code, 
                      Country,
                      region_name,
                      Publ_year,
                      redlist_source)


View(drl)
unique(drl$Redlist_cat)

# create a vector of the Red List categories ordered from least to most severe
priority_order <- c("NE", "LC", "NT", "R", "DD", "V", "WL", "threatened", "VU", "EN",  
                    "CR*", "CR", "CR (PEW)", "CR (PE)", "CR(PE)", "RE?", 
                    "REW", "RE", "EW", "PEW", "PRE", "PE", "EX")

# convert this to a factor in the data frame, applying the priority order
drl_sorted <- drl %>%
  mutate(Redlist_cat = factor(Redlist_cat, levels = priority_order))

# arrange by Country, Year, and Redlist_cat (the factor will keep this in priority order)
drl_sorted <- drl %>%
  arrange(ISO3_Code, Publ_year, accepted_taxon_name, Redlist_cat) 

# select the least severe threat category for each Country-Year group
drl_sorted <- drl_sorted %>%
  group_by(ISO3_Code, Publ_year, accepted_taxon_name, ) %>%
  slice(1) %>%
  ungroup()



# final look
drl_sorted %>% group_by(Country, ISO3_Code, region_name, ISO2_Code, Publ_year, redlist_source) %>% 
  summarise(nspec = n_distinct(accepted_taxon_name)) %>% View()
# Thailand with just 7 species assessed can also not count as a country
# with a RL, not ideal, but better than nothing
drl_sorted %>% filter(Country == "Thailand") %>% View()

# bind national and subnational lists together
drl_final <- bind_rows(drl_sorted %>% mutate(ISO2_Code = NA), 
                       subnational_drl %>% 
                         select(accepted_taxon_name,
                                accepted_taxon_authors,
                                status = Threat,
                                Redlist_cat,
                                ISO3_Code,
                                ISO2_Code, 
                                Country,
                                region_name,
                                Publ_year))


# OK now reduce the dimensionality of this data frame
drl_final <- drl_final %>% select(accepted_taxon_name,
                      accepted_taxon_authors,
               status,
               Redlist_cat,
               ISO3_Code,
               ISO2_Code, 
               Country,
               region_name,
               Publ_year) %>% 
  filter(status != "not_threatened") %>% 
  filter(status != "not_evaluated_dd") %>% 
  distinct()


# rls for how many countries?
drl_final %>% select(Country) %>% distinct %>% nrow()
View(drl_final)

# how many species are threatened in this final list?
# but wait for taxonomic harmonization
# and exclusion of subs and vars and hybrids for final number reported in text
drl_final %>% filter(status == "threatened") %>% 
  select(accepted_taxon_name) %>% distinct %>% nrow()
unique(drl_final$Redlist_cat)
drl_final %>% count(Redlist_cat) %>% View()

# R (rare) species should be excluded
drl_final <- drl_final %>% filter(Redlist_cat != "R") 
drl_final %>% View()


# whats the temporal span of final rl data
drl_final %>% select(Country, region_name, Publ_year) %>% 
  distinct %>% 
  mutate(Publ_year = as.numeric(Publ_year)) %>% 
  summary()

drl_final[which.max(drl_final$Publ_year),] %>% View()
drl_final[which.min(drl_final$Publ_year),] %>% View()
View(drl_final)

# how many species with threat status
drl_final %>% filter(status != "near_threatened") %>% count(Redlist_cat)

# define category mappings
iucn_mapping <- list(
  "CR (Critically Endangered)" = c("CR", "CR (PE)", "CR (PEW)", "CR(PE)", "CR*"),
  "EN (Endangered)" = c("EN"),
  "VU (Vulnerable)" = c("VU", "G"),
  "EW (Extinct in the Wild)" = c("EW", "REW", "PEW"),
  "EX (Extinct)" = c("EX", "PEX",  "PRE", "RE", "RE?", "PE"),
  "Threatened" = c("threatened")
)

consolidated_data <- drl_final %>% 
  filter(status != "near_threatened") %>% 
  count(Redlist_cat) %>%
  mutate(IUCN_Category = case_when(
    Redlist_cat %in% iucn_mapping[["CR (Critically Endangered)"]] ~ "CR (Critically Endangered)",
    Redlist_cat %in% iucn_mapping[["EN (Endangered)"]] ~ "EN (Endangered)",
    Redlist_cat %in% iucn_mapping[["VU (Vulnerable)"]] ~ "VU (Vulnerable)",
    Redlist_cat %in% iucn_mapping[["EW (Extinct in the Wild)"]] ~ "EW (Extinct in the Wild)",
    Redlist_cat %in% iucn_mapping[["EX (Extinct)"]] ~ "EX (Extinct)",
    Redlist_cat %in% iucn_mapping[["Threatened"]] ~ "Threatened",
    TRUE ~ NA_character_  # For any categories not included in the mapping
  )) %>%
  filter(!is.na(IUCN_Category)) %>%
  group_by(IUCN_Category) %>%
  summarise(Count = sum(n)) %>%
  arrange(desc(Count))

# ok 
sum(consolidated_data$Count)
drl_final %>% filter(status != "near_threatened") %>% count(Redlist_cat) %>% summarise(sum(n))
nrow(drl_final %>% filter(status == "near_threatened"))

write.csv(drl_final, "Data/redlist_checked.csv", row.names = F)
