source("R-code/00-preamble.R")

# load data ---------------------------------------------------------------
# red list data
drl <- read_csv("Data/redlist_checked.csv")
nrow(drl)
# View(drl)

# glonaf data
dgl <- read_delim("Data/GLONAF/Taxon_x_List_GloNAF_vanKleunenetal2018Ecology.csv", 
                  delim = "\t", escape_double = FALSE, 
                  locale = locale(encoding = "utf-16"), 
                  trim_ws = TRUE)
names(dgl)


# taxonomic harmonization and further data cleaning -----------------------


#### red list data
# select only rl species and authors, make unique
drl_species <- drl %>% select(accepted_taxon_name, accepted_taxon_authors) %>% 
  distinct

# match with wcvp
drl_wcvp <- wcvp_match_names(drl_species, 
                             name_col = "accepted_taxon_name",
                             author_col = "accepted_taxon_authors",
                             fuzzy = F)

# multiple matches
drl_wcvp %>% count(accepted_taxon_name) %>% arrange(desc(n)) %>% print(n=20)

# if multiple matches keep accepted
drl_wcvp_acc <- drl_wcvp %>% select(accepted_taxon_name, 
                                    wcvp_status, 
                                    wcvp_accepted_id) %>% 
  arrange(desc(wcvp_status %in% c("Accepted", "Synonym")), 
          desc(wcvp_status == "Accepted")) %>%
  distinct(accepted_taxon_name, .keep_all = TRUE)

# check
drl_wcvp_acc %>% count(accepted_taxon_name) %>% arrange(desc(n)) %>% print(n=20)

# join accepted name
names <- wcvp_names %>% select(wcvp_accepted_id = plant_name_id, ipni_id, taxon_name)
drl_wcvp_acc <- left_join(drl_wcvp_acc, names) 

# join back to original data frame
drl <- left_join(drl, drl_wcvp_acc)

# how many species?
drl %>% select(taxon_name) %>% distinct %>% nrow()

# exclude hybrids, subspecies and varieties
drl <- drl %>% 
  left_join(wcvp_names %>% 
              select(plant_name_id, species_hybrid),
            by = c("wcvp_accepted_id" = "plant_name_id")) %>% 
  filter(is.na(species_hybrid)) %>% select(-species_hybrid)

drl %>% select(taxon_name) %>% distinct %>% nrow()
36710 - 36549

# remove unrecognized taxa
drl %>% filter(is.na(taxon_name)) # none

# how many subsp. and varieties in data, remove
drl %>% filter(grepl("subsp\\.|var\\.", taxon_name)) %>% nrow
drl %>% 
  filter(grepl("subsp\\.|var\\.", taxon_name)) %>% 
  select(taxon_name) %>% 
  distinct

# remove as this is impossible to match globally with RLs (see Humphrey et al NEE)
drl <- drl %>% filter(!grepl("subsp\\.|var\\.", taxon_name))

# also flag apomictic taxa
# download latest version of the Apomixis Database from: 
# https://uni-goettingen.de/en/433689.html
# (i) click on "Search" in query box. 
# (ii) scroll down in species list. 
# (iii) click "Export Report as CSV" at bottom of species list. 

# read text file as downloaded from URL
apomixis <- read.table("Data/search_result.txt", header = T, sep = "")  
# View(apomixis)
# filter rows and manipulate Genus column
apomixis <- apomixis %>%
  filter(Apomixis.Yes.Uncertain. == "Y") %>%
  mutate(
    Genus = as.character(Genus),
    Genus1 = sapply(strsplit(Genus, split = ' [()]'), `[`, 1),
    Genus2 = sapply(strsplit(Genus, split = '[()]'), `[`, 2)
  )

# extract unique genera
apomixis_genus <- unique(c(
  apomixis$Genus[!grepl("\\s*\\([^\\)]+\\)", apomixis$Genus)],
  apomixis$Genus1[!is.na(apomixis$Genus1)],
  apomixis$Genus2[!is.na(apomixis$Genus2)]
))

# extract genus from species and flag apomictic genera
drl <- drl %>%
  mutate(
    work_genus = sapply(taxon_name, function(x) strsplit(x, " ")[[1]][1]),
    apomictic = ifelse(work_genus %in% apomixis_genus, 1, 0)
  ) %>% select(-work_genus)


# how many apomictic species?
drl %>% filter(apomictic == 1) %>% select(taxon_name) %>% distinct %>% nrow()

drl <- drl %>% filter(apomictic == 0)

# how many unique species
drl %>% select(taxon_name) %>% distinct %>% nrow
# threatened
drl %>% select(taxon_name, status, ISO3_Code, ISO2_Code) %>% distinct %>% 
  filter(status == "threatened") %>% 
  select(taxon_name) %>% distinct %>% nrow

# exclusively near-threatened
drl %>%
  select(taxon_name, status, ISO3_Code, ISO2_Code) %>% 
  distinct() %>% 
  group_by(taxon_name) %>% 
  # Filter to include only species with all statuses as "near_threatened"
  filter(all(status == "near_threatened")) %>%
  ungroup() %>% 
  select(taxon_name) %>% 
  distinct() %>%
  nrow()

write.csv(drl, "Data/drl-harmonized.csv", row.names = FALSE)


#### glonaf data
# select only glonaf species and authors, make unique
dgl_species <- dgl %>% select(standardized_name, author) %>% 
  distinct

# match with wcvp
dgl_wcvp <- wcvp_match_names(dgl_species, 
                             name_col = "standardized_name",
                             author_col = "author",
                             fuzzy = F)

# multiple matches
dgl_wcvp %>% count(standardized_name) %>% arrange(desc(n)) %>% print(n=20)

# if multiple matches keep accepted
dgl_wcvp_acc <- dgl_wcvp %>% select(standardized_name, 
                                    wcvp_status, 
                                    wcvp_accepted_id) %>% 
  arrange(desc(wcvp_status %in% c("Accepted", "Synonym")), 
          desc(wcvp_status == "Accepted")) %>%
  distinct(standardized_name, .keep_all = TRUE)

# check
dgl_wcvp_acc %>% count(standardized_name) %>% arrange(desc(n)) %>% print(n=20)

# join accepted name
names <- wcvp_names %>% select(wcvp_accepted_id = plant_name_id, ipni_id, taxon_name)
dgl_wcvp_acc <- left_join(dgl_wcvp_acc, names) 

# join back to original data frame
dgl <- left_join(dgl, dgl_wcvp_acc)

write.csv(dgl, "Data/dgl-harmonized.csv", row.names = FALSE)


# join glonaf data and red list data --------------------------------------

#### prepare glonaf 
# only naturalized species
dgl <- read_csv("Data/dgl-harmonized.csv")

# how many species
dgl %>% select(taxon_name) %>% distinct %>% nrow

# only include naturalized species
dgl <- dgl %>% 
  filter(status == "naturalized") 

dgl %>% select(taxon_name) %>% distinct %>% nrow
13081 - 12315

# remove hybrids, as these are not really distinct species
dgl <- dgl %>% 
  left_join(wcvp_names %>% 
              select(plant_name_id, species_hybrid),
            by = c("wcvp_accepted_id" = "plant_name_id")) %>% 
  filter(is.na(species_hybrid)) %>% select(-species_hybrid)
dgl %>% select(taxon_name) %>% distinct %>% nrow
12315 - 11954

# remove unrecognized taxa
dgl %>% filter(is.na(taxon_name))
dgl %>% filter(is.na(taxon_name)) %>% select(standardized_name) %>% distinct %>% nrow
dgl <- dgl %>% filter(!is.na(taxon_name))
dgl %>% select(taxon_name) %>% distinct %>% nrow

# how many subsp. and varieties in data, remove
dgl %>% filter(grepl("subsp\\.|var\\.", taxon_name)) %>% select(taxon_name) %>% 
  distinct %>% nrow()

# remove as this is impossible to match globally with RLs (see Humphrey et al NEE)
dgl <- dgl %>% filter(!grepl("subsp\\.|var\\.", taxon_name))

# remove and flag apomictic genera
dgl <- dgl %>%
  mutate(
    work_genus = sapply(taxon_name, function(x) strsplit(x, " ")[[1]][1]),
    apomictic = ifelse(work_genus %in% apomixis_genus, 1, 0)
  ) %>% select(-work_genus)

# how many apomictic species?
dgl %>% filter(apomictic == 1) %>% select(taxon_name) %>% distinct %>% nrow()
dgl <- dgl %>% filter(apomictic == 0)

# select most important cols
dgl <- dgl %>% 
  select(taxon_name, wcvp_accepted_id, status, region_id)

# how many unique species at last in Glonaf that we will work with?
dgl %>% select(taxon_name) %>% distinct %>% nrow


# instead of region, join iso3 code to glonaf, and for Russia, USA and Canada
# we could have more detailed info, but for Russia really difficult to match 
# Glonaf and RLs, thus simply take ISO3 here
glonaf_region <- read_csv("Data/GLONAF/Region_GloNAF_vanKleunenetal2018Ecology.csv")
dgl_iso3 <- left_join(dgl, glonaf_region %>% select(region_id, country_ISO))              
dgl_iso3 <- dgl_iso3 %>% select(-region_id) %>% rename(ISO3_Code = country_ISO)

# but also create a seperate dataset where we have the tdwg3 level abbreviations,
# we will use this resolution when we calculate species non-native AOOs,
# which is more precise
dgl_tdwg3 <- left_join(dgl, glonaf_region %>% select(region_id, tdwg3))              
dgl_tdwg3 <- dgl_tdwg3 %>% select(-region_id)
# write because we will not continue wokring with these data for now...
write.csv(dgl_tdwg3, "Data/glonaf_tdwg3.csv", row.names = F)


# go back to iso3 glonaf
# restructure status column so that you can match it with the rl cat data below
dgl <- dgl_iso3 %>% rename(status_detail = status) %>% 
  mutate(status = "non-native") %>% 
  select(-status_detail) %>% 
  # because sometimes the region links to mutiple ISOs of the same when the
  # assessments are subnational, for example Germany has statelevel non-native
  # data but not a national list (in Glonaf), so make distinct
  distinct

# how many countries data are available
dgl %>% select(ISO3_Code) %>% distinct %>% nrow


#### prepare red list data
drl <- read_csv("Data/drl-harmonized.csv")

# now simply use ISO-3 aggregation
drl <- drl %>% select(taxon_name, 
                            wcvp_accepted_id, 
                            status, 
                            ISO3_Code) %>% 
  distinct


# bind rows
d <- bind_rows(drl, dgl) %>% arrange(taxon_name)
head(d)
nrow(d)
View(d)


# make plot for spatial data availability ---------------------------------

# in which countries do we have data for threatened and non-native
threatened_data <- d %>% select(status, ISO3_Code) %>% 
  distinct  %>% 
  filter(status == "threatened") %>% 
  mutate(status = "Red List data")
non_native_data <- d %>% select(status, ISO3_Code) %>% 
  distinct  %>% 
  filter(status == "non-native")%>% 
  mutate(status = "Naturalization data")

# join with world data (loaded in preamble)
world_threatened <- left_join(world, threatened_data, by = c("iso_a3" = "ISO3_Code"))
world_threatened$status[is.na(world_threatened$status)] <- "Missing"
world_non_native <- left_join(world, non_native_data, by = c("iso_a3" = "ISO3_Code"))
world_non_native$status[is.na(world_non_native$status)] <- "Missing"

custom_theme <- theme_minimal(base_family = "Arial Narrow") +
  theme(
    legend.position = "top",
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank(), 
    legend.margin = margin(0, 0, 0, 0),
    axis.text = element_blank(),
    legend.box.margin = margin(0, 0, 0, 0)
  )

# plot threatened species
cc_threatened <- c("Red List data" = "#FF6F61")
plot_threatened <- ggplot() +
  geom_sf(data = world_threatened, aes(fill = status), size = 0.1, col = NA) +
  scale_fill_manual(values = cc_threatened, name = "", na.value = "grey80") +
  coord_sf(crs = "+proj=robin") +
  custom_theme

# plot non-native species
cc_non_native <- c("Naturalization data" = "#a3ffda")
plot_non_native <- ggplot() +
  geom_sf(data = world_non_native, aes(fill = status), size = 0.1, col = NA) +
  scale_fill_manual(values = cc_non_native, name = "", na.value = "grey80") +
  coord_sf(crs = "+proj=robin") +
  custom_theme

# Combine plots using patchwork
plot_threatened / plot_non_native +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Save combined plot
# save plot
showtext_opts(dpi=600)
ggsave("Figures/data_availabiltiy.png", dpi = 600, bg = "white", height = 5, width = 5)
showtext_opts(dpi=96)


# double classification ---------------------------------------------------

# now there is following issues, sometimes red lists also mention
# non-native species. that is these sp have an assessment
# for example Datura stramonium is non-native everywhere in Europe
# but switzerland assessed it once as threatened in 2002

# for non-native
d %>%
  group_by(taxon_name, ISO3_Code) %>%
  filter(any(status == "threatened") & any(status == "non-native")) %>%
  ungroup() %>% 
  View

# how many species?
d %>%
  group_by(taxon_name, ISO3_Code) %>%
  filter(any(status == "threatened") & any(status == "non-native")) %>%
  ungroup() %>% 
  select(taxon_name) %>% distinct %>% nrow
# that's a lot

# well this shows another problem. Abies alba, for example in Germany, is
# in some part native (and threatened  there) and in others (north) non-native
# so it is in fact possible to have both statuses within one country

# we have to find a way to disentangle these two cases: RLs listing non-natives,
# and countries comprising both genuinely a native threatened species, which is in
# another part of the country non-native

# ok simplify, which species and where is mismatch
mismatch <- d %>%
  group_by(taxon_name, ISO3_Code) %>%
  filter(any(status == "threatened") & any(status == "non-native")) %>%
  ungroup() %>% 
  select(taxon_name, ISO3_Code) %>% distinct

# can exclude reunion because as an island its simply too small for that to
# occur + rnaturalearth does not list REU as an iso3, weirdly
species_exclude <- mismatch %>% filter(ISO3_Code == "REU")
mismatch <- mismatch %>% filter(ISO3_Code != "REU")

# use powo wcvp_distribution shapefiles to intersect for each species
# the mismatch country shapefile if there is an intersection, keep, if not
# we need to remove this species
mismatch$overlap <- NA 
for (i in 1:nrow(mismatch)) {
  
  # species iso3 polygons as of our data
  mismatch_proj <- world %>% 
    select("iso_a3") %>% 
    left_join(mismatch[i,], by = c("iso_a3" = "ISO3_Code")) %>% 
    filter(taxon_name == mismatch$taxon_name[i])%>% 
    st_transform(3857) %>%
    st_make_valid()
  
  # species polygons as of wcvp data
  powo_proj <- wcvp_distribution(mismatch$taxon_name[i], taxon_rank = "species") %>% 
    filter(occurrence_type == "native") %>% 
    select(LEVEL3_COD) %>% 
    st_transform(3857) %>%
    st_make_valid()
  
  # check if there is any overlap (intersection)
  intersection <- st_intersection(powo_proj, mismatch_proj)
  
  # check if there's an overlap and update the 'overlap' column
  if (nrow(intersection) > 0) {
    mismatch$overlap[i] <- TRUE  # There is overlap
  } else {
    mismatch$overlap[i] <- FALSE  # No overlap
  }
  
  print(i)
}

print(mismatch)

# how many have genuine overlap?
mismatch %>% filter(overlap == "TRUE") %>% select(taxon_name) %>% distinct %>% nrow

# species that need exclusion
species_exclude <- mismatch %>% filter(overlap == "FALSE") %>% 
  bind_rows(species_exclude)

# save files because it takes a while to generate
# write.csv(mismatch, "Data/mismatch.csv", row.names = F)
# write.csv(species_exclude, "Data/species_exclude.csv", row.names = F)

mismatch <- read_csv("Data/mismatch.csv")
species_exclude <- read_csv("Data/species_exclude.csv")

# remove species country combos where there is no overlap
d <- d %>%
  anti_join(species_exclude %>% select(-overlap) %>% mutate(status = "threatened"))
nrow(d)
View(d)


# now there might be one issue in particular, the old swiss RL list
# some species and these species are in fact also non-native there
# but not listed by glonaf. This is the only list, to my knowledge,
# where this issue emerges
# let's clean manually.
d <- d %>% 
  filter(!(taxon_name == "Scilla luciliae" & status == "threatened" & ISO3_Code == "CHE")) %>% 
  filter(!(taxon_name == "Scilla forbesii" & status == "threatened" & ISO3_Code == "CHE"))



# create a separate dataset where countries are the same ------------------

# right now the glonaf dataset holds more countries than the new rl dataset
# check
d %>% group_by(ISO3_Code) %>% count(status) %>% spread(status, n) %>% 
  filter(is.na(`non-native`)|is.na(threatened)) %>% View()

# for how many countries do we have naturalization data
d %>% group_by(ISO3_Code) %>% count(status) %>% spread(status, n) %>% 
  filter(!is.na(`non-native`)) %>% filter(!is.na(ISO3_Code)) %>% nrow()
# for how many countries do we have red list data
d %>% group_by(ISO3_Code) %>% count(status) %>% spread(status, n) %>% 
  filter(!is.na(threatened)) %>% filter(!is.na(ISO3_Code)) %>% nrow()

# verify with original data
drl %>% select(ISO3_Code) %>% filter(!is.na(ISO3_Code)) %>% distinct %>% nrow
dgl %>% select(ISO3_Code) %>% filter(!is.na(ISO3_Code)) %>% distinct %>% nrow
# matches!

# ok remove countries where just one data type exists
d_countrysame <- d %>% group_by(ISO3_Code) %>% count(status) %>% spread(status, n) %>% 
  filter(!(is.na(`non-native`)|is.na(threatened))) %>% select(ISO3_Code) %>% 
  left_join(d)


# now only species that are non-native and threatened --------------------

# perform this for both datasets separately

# all data
# exclude species that are ONLY threatened, these have not naturalized anywhere

# for this purpose create again a dataset that excludes near-threatned species
# and one that assumes that they also count as somewhat threatened species

d_all <- d %>% mutate(status = ifelse(status=="near_threatened", 
                                      "threatened",
                                      status))
nrow(d_all)

d_all %>%
  group_by(taxon_name) %>%
  summarize(status_unique = list(unique(status))) %>% 
  select(status_unique) %>% distinct %>% View()

species_to_remove <- d_all %>%
  group_by(taxon_name) %>%
  summarize(status_unique = list(unique(status))) %>% 
  filter(status_unique == "threatened") %>% 
  pull(taxon_name)

#  filter out the rows where species are only threatened
d_final_all <- d_all %>%
  filter(!taxon_name %in% species_to_remove)

nrow(d_final_all)
# for how many non-native, invasive and threatened species?
d_final_all %>% filter(status == "threatened") %>% select(taxon_name) %>% distinct %>% nrow()
d_final_all %>% filter(status == "non-native") %>% select(taxon_name) %>% distinct %>% nrow()

write.csv(d_final_all, "Data/non-native-red-list-species-NTincluded.csv", row.names = F)

d_onlythreatened <- d %>% filter(status != "near_threatened")
nrow(d_onlythreatened)

d_onlythreatened %>%
  group_by(taxon_name) %>%
  summarize(status_unique = list(unique(status))) %>% 
  select(status_unique) %>% distinct %>% View()

species_to_remove <- d_onlythreatened %>%
  group_by(taxon_name) %>%
  summarize(status_unique = list(unique(status))) %>% 
  filter(status_unique == "threatened") %>% 
  pull(taxon_name)

#  filter out the rows where species are only threatened
d_final_onlythreatened <- d_onlythreatened %>%
  filter(!taxon_name %in% species_to_remove)
nrow(d_final_onlythreatened)

# for how many non-native, invasive and threatened species?
d_final_onlythreatened %>% filter(status == "threatened") %>% select(taxon_name) %>% distinct %>% nrow()
d_final_onlythreatened %>% filter(status == "non-native") %>% select(taxon_name) %>% distinct %>% nrow()

write.csv(d_final_onlythreatened, "Data/non-native-red-list-species-NTexcluded.csv", row.names = F)


# country same data
# exclude species that are ONLY threatened, these have not naturalized anywhere
# so far
d_countrysame <- d_countrysame %>% filter(status != "near_threatened")

species_to_remove <- d_countrysame %>%
  group_by(taxon_name) %>%
  summarize(status_unique = list(unique(status))) %>% 
  filter(status_unique == "threatened") %>% 
  pull(taxon_name)

#  filter out the rows where species are only threatened
d_countrysame <- d_countrysame %>%
  filter(!taxon_name %in% species_to_remove)

nrow(d_countrysame)
head(d_countrysame)

write.csv(d_countrysame, "Data/non-native-red-list-species-countrysame.csv", row.names = F)

