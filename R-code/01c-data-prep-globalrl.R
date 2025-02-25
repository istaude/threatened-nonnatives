source("R-code/00-preamble.R")


# load data ---------------------------------------------------------------

# iucn global rl and only vascular plants
diucn <- read_csv("Data/simple_summary.csv") %>% filter(phylumName == "TRACHEOPHYTA")
names(diucn)
str(diucn)
unique(diucn$phylumName)

# taxonomic harmonization -------------------------------------------------

diucn_species <- diucn %>% select(scientificName, authority) %>% 
  distinct

# match with wcvp
diucn_wcvp <- wcvp_match_names(diucn_species, 
                             name_col = "scientificName",
                             author_col = "authority",
                             fuzzy = T)


# multiple matches
diucn_wcvp %>% count(scientificName) %>% arrange(desc(n)) %>% print(n=20)

# if multiple matches keep accepted
diucn_wcvp_acc <- diucn_wcvp %>% select(scientificName, 
                                    wcvp_status, 
                                    wcvp_accepted_id) %>% 
  arrange(desc(wcvp_status %in% c("Accepted", "Synonym")), 
          desc(wcvp_status == "Accepted")) %>%
  distinct(scientificName, .keep_all = TRUE)

# check
diucn_wcvp_acc %>% count(scientificName) %>% arrange(desc(n)) %>% print(n=20)

# join accepted name
names <- wcvp_names %>% select(wcvp_accepted_id = plant_name_id, ipni_id, taxon_name)
diucn_wcvp_acc <- left_join(diucn_wcvp_acc, names) 
View(diucn_wcvp_acc)

# join back to original data frame
diucn <- left_join(diucn, diucn_wcvp_acc)
View(diucn)
# how many species?
diucn %>% select(taxon_name) %>% distinct %>% nrow()

# exclude hybrids, subspecies and varieties
diucn <- diucn %>% 
  left_join(wcvp_names %>% 
              select(plant_name_id, species_hybrid),
            by = c("wcvp_accepted_id" = "plant_name_id")) %>% 
  filter(is.na(species_hybrid)) %>% select(-species_hybrid)

diucn %>% select(taxon_name) %>% distinct %>% nrow()
71996 - 71930

# remove unrecognized taxa, needs fuzzy matching apparently
diucn %>% filter(is.na(taxon_name)) %>% View()
# if this doesn't even work with fuzzy matching we must remove
diucn <- diucn %>% filter(!is.na(taxon_name))

# how many subsp. and varieties in data, remove
diucn %>% filter(grepl("subsp\\.|var\\.", taxon_name)) %>% nrow
diucn %>% 
  filter(grepl("subsp\\.|var\\.", taxon_name)) %>% 
  select(taxon_name) %>% 
  distinct

# remove as this is impossible to match globally with RLs (see Humphrey et al NEE)
diucn <- diucn %>% filter(!grepl("subsp\\.|var\\.", taxon_name))

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
diucn <- diucn %>%
  mutate(
    work_genus = sapply(taxon_name, function(x) strsplit(x, " ")[[1]][1]),
    apomictic = ifelse(work_genus %in% apomixis_genus, 1, 0)
  ) %>% select(-work_genus)


# how many apomictic species?
diucn %>% filter(apomictic == 1) %>% select(taxon_name) %>% distinct %>% nrow()
diucn <- diucn %>% filter(apomictic == 0)

# how many unique species
diucn %>% select(taxon_name) %>% distinct %>% nrow
# threatened
diucn %>% select(taxon_name, redlistCategory) %>% distinct %>% 
  count(redlistCategory) %>% 
  filter(redlistCategory == "Extinct" |
           redlistCategory == "Extinct in the Wild" |
           redlistCategory == "Critically Endangered" |
           redlistCategory == "Endangered" |
           redlistCategory == "Vulnerable"
           ) %>% 
  mutate(perc = n/24523)

# near-threatened
diucn %>% select(taxon_name, redlistCategory) %>% distinct %>% 
  count(redlistCategory) %>% 
  filter(redlistCategory == "Near Threatened")



# integrate with glonaf data ----------------------------------------------

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
  select(taxon_name, wcvp_accepted_id, status) %>% distinct

# how many unique species at last in Glonaf that we will work with?
dgl %>% select(taxon_name) %>% distinct %>% nrow

# join with global IUCN RL
diucn_gl <- left_join(dgl, diucn %>% select(taxon_name, redlistCategory) %>% distinct)

# save file
write.csv(diucn_gl, "Data/globalrl-non-natives.csv", row.names = FALSE)


# inspect data ------------------------------------------------------------

diucn_gl %>% head
diucn_gl %>% count(redlistCategory)  %>% 
  filter(redlistCategory == "Extinct" |
           redlistCategory == "Extinct in the Wild" |
           redlistCategory == "Critically Endangered" |
           redlistCategory == "Endangered" |
           redlistCategory == "Vulnerable"
  ) %>% 
  summarize(sum_threatened = sum(n)) %>% 
  mutate(perc = sum_threatened / 9195)

# if near-threatened are included
diucn_gl %>% count(redlistCategory)  %>% 
  filter(redlistCategory == "Extinct" |
           redlistCategory == "Extinct in the Wild" |
           redlistCategory == "Critically Endangered" |
           redlistCategory == "Endangered" |
           redlistCategory == "Vulnerable" |
           redlistCategory == "Near Threatened"
  ) %>% 
  summarize(sum_threatened = sum(n)) %>% 
  mutate(perc = sum_threatened / 9195)
