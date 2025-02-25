source("R-code/00-preamble.R")

# load data ---------------------------------------------------------------
# only threatened non-natives
d <- read_csv("Data/non-native-red-list-species-NTexcluded.csv")

df <- d %>%
  group_by(taxon_name, wcvp_accepted_id) %>%
  summarize(status_unique = list(unique(status))) %>% 
  mutate(status_unique = ifelse(status_unique == 'c("threatened", "non-native")',
                                "threatenednonnative", "nonnative")
  )

threatened_and_nonnative <- df %>% 
  filter(status_unique == "threatenednonnative")

# because we want to classify ranges at the smallest resolution possible
# load RL data where subnational assessments are still preserved for the
# USA, Russia and Canada
drl_harmonized <- read_csv("Data/drl-harmonized.csv")
rl_species <- left_join(d %>% 
                          filter(status == "threatened") %>% 
                          select(taxon_name, ISO3_Code) %>% 
                          distinct,
                        drl_harmonized %>% 
                          select(taxon_name, ISO3_Code, ISO2_Code))

# also load a subnational iso2 shapefile that matches with our ISO2 Codes
subnational_shp <- ne_states(returnclass = "sf")
# only select subnational shapes for USA, RUS and Can
subnational_shp <- subnational_shp %>% 
  filter(adm0_a3 == "USA"|adm0_a3 == "RUS"|adm0_a3 == "CAN") %>% 
  rename(ISO2_Code=iso_3166_2)

# make shape valid
subnational_shp <- subnational_shp %>% 
  st_transform(8857) %>% 
  st_make_valid()


# correct rl iso2s
rl_species$ISO2_Code <- gsub("USA", "US", rl_species$ISO2_Code)
rl_species$ISO2_Code <- gsub("RUS", "RU", rl_species$ISO2_Code)
rl_species$ISO2_Code <- gsub("CAN", "CA", rl_species$ISO2_Code)

# check for USA if fits
left_join(
  rl_species %>% filter(ISO3_Code == "USA") %>% select(ISO2_Code) %>% distinct,
  subnational_shp %>% filter(adm0_a3 == "USA") 
)
# check for RUS if fits
left_join(
  rl_species %>% filter(ISO3_Code == "RUS") %>% select(ISO2_Code) %>% distinct,
  subnational_shp %>% filter(adm0_a3 == "RUS")
)
# check for CAN if fits
left_join(
  rl_species %>% filter(ISO3_Code == "CAN") %>% select(ISO2_Code) %>% distinct,
  subnational_shp %>% filter(adm0_a3 == "CAN")
)
# perfect, the rl ISO2 is now harmonized with our subnational_shp
# also reproject world shapefile, better here than in the loop
world <- world %>% st_transform(8857)


# also load the Glonaf data (shapefiles) so that this doesn't happen at the
# ISO3 level either
glonaf_shp <- st_read("Data/GLONAF/shapefile/regions2.shp")
# make shape valid
glonaf_shp <- glonaf_shp %>% 
  st_transform(8857) %>% 
  st_make_valid()

glonaf_region <- read_csv("Data/GLONAF/Region_GloNAF_vanKleunenetal2018Ecology.csv")
# select region id, which links to glonaf_sp, and OBJDsic which links to
# the shapefile
glonaf_region <- glonaf_region %>% select(region_id, OBJIDsic) %>% distinct

# load harmonized glonaf data
dgl_harmonized <- read_csv("Data/dgl-harmonized.csv")
# select only taxon_name and region_id
glonaf_sp <- dgl_harmonized %>% select(taxon_name, region_id) %>% 
  # crop to the selected species in the data, df
  right_join(df %>% select(taxon_name) %>% distinct)
# join OBJDsic via region connector
glonaf_sp <- glonaf_sp %>% left_join(glonaf_region)


# loop for all species ----------------------------------------------------

# from here, run the script best over night, takes a bit to loop over all species
# here we retrieve gbif occurrence records and classify them into
# native, threatened and non-native

# Set up the species vector
species <- unique(threatened_and_nonnative$taxon_name)


# Initialize an empty data frame to store results
species_keys_df <- data.frame(species = character(), usageKey = integer(), stringsAsFactors = FALSE)

# Loop over each species to find GBIF usage keys
for (i in 1:length(species)) {
  # Retrieve the usageKey from GBIF
  taxon_result <- name_backbone(species[i], rank = "species")
  
  # Extract the usageKey
  usage_key <- ifelse("usageKey" %in% names(taxon_result), taxon_result$usageKey, NA)
  
  # Add species name and usageKey to the dataframe
  species_keys_df <- rbind(species_keys_df, data.frame(species = species[i], 
                                                       usageKey = usage_key))
  # Print progress
  print(paste("Processed species", i, "of", length(species)))
}

# Check if all species have keys
species_keys_df %>% filter(is.na(usageKey))

# Manually add the usageKeys for the two species
species_keys_df[species_keys_df$species == "Acmella pusilla", "usageKey"] <- 5386762

# Write and read
# write.csv(species_keys_df, "Data/species_keys_df.csv", row.names = FALSE)
species_keys_df <- read_csv("Data/species_keys_df.csv")
nrow(species_keys_df)

#which(species_keys_df$species == "x")
# Initialize an empty list to store results
all_sf_combined <- list()
# i = 1 # trial species
# Loop over each species
for (i in 1:nrow(species_keys_df)) {
  
  species_name <- species_keys_df$species[i]
  # Fetch species polygons from GBIF
  gbif_proj <- mvt_fetch(taxonKey = species_keys_df$usageKey[i], 
                                      bin = "square", 
                                      squareSize = 8,
                                      hasCoordinateIssue = FALSE,
                                      basisOfRecord = c("HUMAN_OBSERVATION", 
                                                        "MACHINE_OBSERVATION", 
                                                        "OBSERVATION"))
  
  # Check if the fetch result is empty or NULL
  if (is.null(gbif_proj) || nrow(gbif_proj) == 0) {
    print(paste("No GBIF data for species", species_keys_df$species[i], "- skipping to next."))
    next  # Skip to the next iteration of the loop
  }
  
  gbif_proj <- gbif_proj %>% 
    mutate(species = species_name) %>%
    st_transform(8857) %>%
    st_make_valid()
  
  # Crop with POWO native ranges
  powo_proj <- wcvp_distribution(species_name, taxon_rank = "species") %>% 
    filter(occurrence_type == "native") %>% 
    select(LEVEL3_COD) %>% 
    st_transform(8857) %>% 
    st_make_valid()
  
  gbif_cropped_native <- gbif_proj %>% 
    st_intersection(powo_proj)  %>% 
    mutate(status = "native")
  
  # Crop with threatened countries, but consider that some countries
  # have ISO2 (subnational) threat data
  threatened_countries <- rl_species %>% 
    filter(taxon_name == species_name)
  
  # Process subnational and national data separately
  threatened_countries_iso2 <- threatened_countries %>% filter(!is.na(ISO2_Code))
  threatened_countries_iso3 <- threatened_countries %>% filter(is.na(ISO2_Code))
  
  # Handle ISO2 level data (subnational regions)
  if (nrow(threatened_countries_iso2) > 0) {
    threatened_subnational <- subnational_shp %>%
      select(ISO2_Code) %>%
      left_join(threatened_countries_iso2) 
    
    gbif_cropped_threatened_iso2 <- gbif_proj %>% 
      st_intersection(threatened_subnational %>% filter(taxon_name == species_name)) %>% 
      mutate(status = "threatened")
    
  } else {
    gbif_cropped_threatened_iso2 <- NULL
  }
  
  # Handle ISO3 level data (national regions)
  if (nrow(threatened_countries_iso3) > 0) {
    threatened_countries_proj <- world %>% 
      select(ISO3_Code = iso_a3) %>% 
      left_join(threatened_countries_iso3)
    
    gbif_cropped_threatened_iso3 <- gbif_proj %>% 
      st_intersection(threatened_countries_proj %>% filter(taxon_name == species_name)) %>% 
      mutate(status = "threatened")
  } else {
    gbif_cropped_threatened_iso3 <- NULL
  }
  
  # Combine ISO2 and ISO3 results (subnational takes precedence if both exist)
  gbif_cropped_threatened <- rbind(gbif_cropped_threatened_iso2, gbif_cropped_threatened_iso3)
  
  
  # Crop with non-native GLONAF regions
  non_native_countries <- glonaf_sp %>% 
    filter(taxon_name == species_name)
  non_native_countries_proj <- glonaf_shp %>% 
    select("OBJIDsic") %>% 
    left_join(non_native_countries, by = "OBJIDsic")

  gbif_cropped_nonnative <- gbif_proj %>% 
    st_intersection(non_native_countries_proj %>% filter(taxon_name == species_name)) %>% 
    mutate(status = "non-native")
  

  # Combine POWO range with RL range = full native range
  sf_combined <- rbind(gbif_cropped_native %>% select(species, status), 
                       gbif_cropped_threatened %>% select(species, status),
                       gbif_cropped_nonnative %>% select(species, status))
  
  # Add to the overall list
  all_sf_combined[[i]] <- sf_combined
  
  # Print progress
  print(paste("Processed species", i, "of", length(species), species_name))
}


# Save the list to a file
saveRDS(all_sf_combined, "Data/species-shps/all_sf_combined.rds")

# Load the list back into R
all_sf <- readRDS("Data/species-shps/all_sf_combined.rds")
length(all_sf)


# how many species didnt have any gbif data?
species_keys_df[sapply(all_sf, is.null),]
# some of these like Arthraxon castratus do have occ records
# on GBIF, but these are preserved specimen records, and thus
# were not included here

# filter these species out
all_sf <- Filter(Negate(is.null), all_sf)
length(all_sf)


# which species
species_list <- lapply(all_sf, function(x) unique(x$species))
length(species_list)
null_or_empty <- sapply(species_list, function(x) is.null(x) || length(x) == 0)
if (any(null_or_empty)) {
  print("There are NULL or empty elements in species_list.")
  print(which(null_or_empty))  # Locations of NULL or empty elements
}

(all_species <- unlist(species_list))
length(all_species)

# ok if empty and null elments are removed 2486 species have some gbif data

# visualize a species as example
# k = 256
species_list[k]
which(species_list == "Barbarea vulgaris")
ggplot() +
  #geom_sf(data = world) +
  geom_sf(data = all_sf[[k]] %>% filter(status != "non-native"), col = "green") +
  geom_sf(data = all_sf[[k]] %>% filter(status == "threatened"), col= "red") +
  geom_sf(data = all_sf[[k]] %>% filter(status == "non-native"), col = "purple")


# calculate mapping resolution for species with latitudinal spread
x <- mvt_fetch(taxonKey = species_keys_df$usageKey[k], 
          bin = "square", 
          squareSize = 8,
          hasCoordinateIssue = FALSE,
          basisOfRecord = c("HUMAN_OBSERVATION", 
                            "MACHINE_OBSERVATION", 
                            "OBSERVATION"))
ggplot() +
  geom_sf(data = world) +
  geom_sf(data = x,  col = "green") 

b <- NULL
for(i in 1:nrow(x)){
  b[i] <-  x[i,] %>% st_area() / 1e6
}
range(b)


# calculate the fraction of species range threatened ----------------------

# Initialize an empty data frame to store the results
results_df <- data.frame(
  species = character(),
  native_area_km2 = numeric(),
  threatened_area_km2 = numeric(),
  ratio = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each element in the list
for (i in 1:length(all_sf)) {
  
  # Get the species name (assuming there's a species column)
  species_name <- unique(all_sf[[i]]$species)
  
  # Check if species_name is empty (character(0))
  if (length(species_name) == 0) {
    print(paste("Skipping species", i, "because species_name is empty."))
    next  # Skip to the next iteration of the loop
  }
  
  # Exclude non-native occurrences -> native range
  native_range <- all_sf[[i]] %>% filter(status != "non-native") %>% 
    st_union() 
  native_area <- sum(st_area(native_range)) # Area in m²
  
  # Just threatened -> threatened range
  threatened_range <- all_sf[[i]] %>% filter(status == "threatened") %>% 
    st_union() 
  threatened_area <- sum(st_area(threatened_range)) # Area in m²
  
  # Calculate the ratio safely by removing units for the comparison
  native_area_km2 <- as.numeric(native_area) / 1e6  # Convert to km²
  threatened_area_km2 <- as.numeric(threatened_area) / 1e6  # Convert to km²
  
  ratio <- ifelse(native_area_km2 > 0, threatened_area_km2 / native_area_km2, NA)
  
  # Store the results in the data frame
  results_df <- rbind(
    results_df,
    data.frame(
      species = species_name,
      native_area_km2 = native_area_km2,
      threatened_area_km2 = threatened_area_km2,
      ratio = ratio,
      stringsAsFactors = FALSE
    )
  )
  
  print(paste("Processed species", i, "of", length(all_sf), species_name))
  
}

write.csv(results_df, "Data/native_threatened_areas.csv", row.names = FALSE)

results_df <- read.csv("Data/native_threatened_areas.csv")
nrow(results_df)
length(all_species)

# Calculate the non-native area -------------------------------------------

# Initialize the results data frame
results_dfnn <- data.frame(species = character(), 
                           nonnative_area_km2 = numeric(), 
                           stringsAsFactors = FALSE)

# Loop over each element in the list
for (i in 1:length(all_sf)) {
  
  # Get the species name (assuming there's a species column)
  species_name <- unique(all_sf[[i]]$species)
  
  # Check if species_name is empty (character(0))
  if (length(species_name) == 0) {
    print(paste("Skipping species", i, "because species_name is empty."))
    next  # Skip to the next iteration of the loop
  }
  
  # Calculate the area for non-native occurrences
  nonnative_range <- all_sf[[i]] %>% filter(status == "non-native") %>% 
    st_union() 
  nonnative_area <- sum(st_area(nonnative_range)) # Area in m²
  
  # Convert the non-native area to km²
  nonnative_area_km2 <- as.numeric(nonnative_area) / 1e6  # Convert to km²
  
  # Store the results in the data frame
  results_dfnn <- rbind(
    results_dfnn,
    data.frame(
      species = species_name,
      nonnative_area_km2 = nonnative_area_km2,
      stringsAsFactors = FALSE
    )
  )
  
  print(paste("Processed species", i, "of", length(all_sf), species_name))
  
}

# Save the results
write.csv(results_dfnn, "Data/nonnative_areas.csv", row.names = FALSE)



# Join with native and threatened df
results_df <- read.csv("Data/native_threatened_areas.csv")
results_dfnn <- read.csv("Data/nonnative_areas.csv")
results_df <- left_join(results_df, results_dfnn)


# Check the data ----------------------------------------------------------

# Now make it a requirement that a species had at least one observation
# in its native, threatened or non-native range part
nrow(results_df)

# how many 
results_df %>% 
  filter(native_area_km2 == 0) %>% nrow()
results_df %>% 
  filter(nonnative_area_km2 == 0) %>% nrow()
results_df %>% 
  filter(threatened_area_km2 == 0) %>% nrow()

filtered_results <- results_df %>% 
  filter(!is.na(nonnative_area_km2)) %>% 
  filter(nonnative_area_km2 != 0) %>% 
  filter(threatened_area_km2 != 0) %>% 
  filter(native_area_km2 != 0)

filtered_results %>% nrow

# Ensure that all_sf and results df match in terms of species
# Extract species names from the dataframe
species_names_in_df <- filtered_results$species

# Filter the list to keep only elements with species in species_df
# all_sf <- readRDS("species-shps/all_sf_combined.rds")

all_sf_red <- lapply(all_sf, function(element) {
  if (any(element$species %in% species_names_in_df)) {
    return(element)  # Keep the element if species matches
  } else {
    return(NULL)  # Otherwise, return NULL
  }
})

all_sf_red <- all_sf_red[!sapply(all_sf_red , is.null)]
# test
all_sf_red[[1421]]
species_names_in_df[1421]

# Save the list to a file
saveRDS(all_sf_red, "Data/species-shps/all_valid_sf.rds")

# Write filtered_resutls
write.csv(filtered_results, "Data/species_areas.csv", row.names = FALSE)

# data visualization ------------------------------------------------------
d <- read.csv("Data/species_areas.csv")
d %>% arrange(desc(ratio))
d %>% arrange((ratio)) %>% View
quantile(d$ratio, na.rm = T)
mean(d$ratio, na.rm = T)
median(d$ratio, na.rm = T)
median(d$nonnative_area_km2, na.rm = T)
median(d$threatened_area_km2, na.rm = T)
#View(d)


# for scaling
histogram_data <- ggplot_build(
  ggplot(d, aes(x = ratio)) +
    geom_histogram(binwidth = 0.05)
)$data[[1]]
max_count <- max(histogram_data$count)
density_scaling_factor <- nrow(d) * 0.05 

# histogram
d %>% 
  ggplot(aes(x = ratio)) +
  geom_histogram(aes(y = after_stat(count)), fill = "#FAC55F", binwidth = 0.05) + 
  geom_density(aes(y = ..density.. * density_scaling_factor), color = "#95846d", size =0.3) + 
  #stat_ecdf(aes(y = ..y.. * max_count), geom = "step", color = "grey30") + 
  scale_y_continuous(
    position = "right"
  ) +
  scale_x_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1), 
                     expand = c(0, 0),
                     labels = c(">0%", "25%", "50%", "75%", "100%")) +
  labs(y ="No. of\nspecies", x = "Native range threatened (%)") +
  theme_minimal(base_family = "Arial Narrow", base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = 13),
        axis.title.y.right = element_text(angle = 0,
                                          lineheight = 0.8,
                                          vjust = 1,
                                          hjust = 0
        )) +
  annotate(
    "segment", 
    x = 0.136, xend = 0.136, 
    y = 110, yend = 450, 
    color = "#95846d", 
    size = 0.3, 
    linetype = "dotted", 
    arrow = arrow(length = unit(0.13, "cm"), type = "closed", angle = 35)  # Solid arrowhead
  ) +
  annotate("text", x = 0.136, y = 475, label = "Average range at risk (13.6%)", 
           color = "#95846d", size = 4.5, family = "Arial Narrow", hjust = 0.1)

ggsave("Figures/histogram.svg", height = 5, width = 5)


# how many have ranges that are 70% threatened
d %>% filter(ratio >= 0.7) %>% nrow /
  d %>% nrow

d %>% filter(ratio >= 0.7) %>% View()

# what fraction is in the first bin?
d %>% filter(ratio <= 0.05) %>% nrow /
  d %>% nrow

# how many have lost more than 90%
d %>% filter(ratio >= 0.9) %>% nrow /
  d %>% nrow


# one:one line with gains and losses in range
# have same range for both x and y axes
range_limit <- range(c(d$nonnative_area_km2, d$threatened_area_km2))

# fit the major axis regression model
ma_mod <- lmodel2(log10(threatened_area_km2) ~ log10(nonnative_area_km2), data = d)

# extract slope and intercept of the major axis regression
ma_slope <- ma_mod$regression.results[ma_mod$regression.results$Method == "MA", "Slope"]
ma_slope_ci <- ma_mod$confidence.intervals[2, c("2.5%-Slope", "97.5%-Slope")]
ma_intercept <- ma_mod$regression.results[ma_mod$regression.results$Method == "MA", "Intercept"]

color_gradient <- scale_fill_gradientn(
  colours = c('#fce0a9', '#fbce78', '#f9bc46', '#f8aa15'),
  breaks = c(0.0, 0.25, 0.5, 0.75, 1),  
  na.value = "grey90"
)

(p0 <- ggplot(d, aes(x = nonnative_area_km2, 
                                  y = threatened_area_km2,
                                  size = ratio, fill = ratio)) +
    geom_point(alpha = 0.6, shape = 21, col = "black") +
    geom_point(alpha = 0.4, shape = 1, col = "white") +
    color_gradient +
    geom_abline(slope = 1, intercept = 1, col = "grey30", lty = "dashed", size = 0.3) +  # 1:1 line
    geom_abline(slope = ma_slope, intercept = ma_intercept, color = '#95846d') +  
    scale_x_log10(limits = range_limit) +  
    scale_y_log10(limits = range_limit) +  
    scale_size_continuous(
      name = "Native range threatened (%)",
      labels = scales::percent_format(accuracy = 1) 
    ) +
    labs(x = 'Range naturalized (km²)',
         y = 'Native range threatened (km²)', tag = "a") + 
    theme_minimal(base_family = "Arial Narrow", base_size = 13) +
    guides(
      size = guide_legend(
        override.aes = list(fill = c('#fce0a9', '#fbce78', '#f9bc46', '#f8aa15'), 
                            alpha = 0.7, shape = 21, color = "white"),
        direction = "horizontal"
      ), 
      fill = "none"
    ) +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 13),
      legend.position = "top",  
      legend.direction = "horizontal", 
      legend.box = "horizontal",  
      legend.title = element_text(size = 12),  
      legend.text = element_text(size = 11)  
    )
)


# proportion net gain vs. partial compensation plot
p1 <- d %>%
  mutate(diff = nonnative_area_km2 - threatened_area_km2) %>% 
  arrange(desc(diff)) %>%
  mutate(rank = 1:nrow(.),
         diff_sign = ifelse(diff >= 0, "positive", "negative"),
         new_range = native_area_km2 + nonnative_area_km2 - threatened_area_km2) %>% 
  ggplot(aes(x = rank, y = diff)) +
  geom_point(alpha = 0.4, shape = 1, col = "white") +
  geom_line(aes(color = diff_sign), size = 0.8, alpha = 1) + 
  geom_ribbon(aes(ymin = 0, ymax = diff, fill = diff_sign), alpha = 0.3) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(labels = scales::label_scientific()) +
  scale_color_manual(values = c("positive" = "#C3C73C", "negative" = '#7C7AA9')) + 
  scale_fill_manual(values = c("positive" = "#C3C73C", "negative" = '#7C7AA9')) +
  labs(y = "Range naturalized - threatened (km²)", tag = "b") +
  theme_minimal(base_family = "Arial Narrow", base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  annotate("text", x = -1, y = Inf * 1.2, label = "Net increase", vjust = 1,
           hjust = 0, size = 5, family = "Arial Narrow",
           col = "#C3C73C") + 
  annotate("text", x = Inf, y = Inf * 1.2, label = "Partial compensation", vjust = 1,
           hjust = 1, size = 5, family = "Arial Narrow", 
           col = '#7C7AA9') 

p1


# how many species do have net gains vs partial compensation
d %>% mutate(type = ifelse(threatened_area_km2 > nonnative_area_km2, 
                           "partial compensation",
                           "net gain")) %>% 
  count(type)

data <- data.frame(
  category = c("Net increase", "Partial compensation"),
  percentage = c(59.6, 40.4),
  group = c("Net increase", "Partial compensation")
)

data$group <- factor(data$group, levels = c("Partial compensation", "Net increase"))
# Create a stacked horizontal bar plot
p2 <- ggplot(data, aes(x = 1, y = percentage, fill = group)) +
  geom_bar(stat = "identity", col = "white", width = 1.5) +
  coord_flip() + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_fill_manual(values = c("Net increase" = "#C3C73C", 
                               "Partial compensation" = "#7C7AA9")) + 
  geom_text(aes(label = paste0(percentage, "%")), family = "Arial Narrow",
            position = position_stack(vjust = 0.5), size = 5, color = "white",
            fontface = "bold") +
  labs(x = "", y = "") +
  theme_void() +
  theme(
    axis.text = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")
  )+
  annotate("text", x = 2.2, y = 0, label = "N = 1,716", 
           hjust = 0, size = 4, family = "Arial Narrow", col = "#95846d")


p0/ p1 / p2 + plot_layout(heights = c(6, 6, 1))

showtext_opts(dpi=600)
ggsave("Figures/species_areas.png", dpi = 600, bg = "white", height = 9, width = 8)
showtext_opts(dpi=96)



# SI
# paired wilcox test
wilcox.test(d$nonnative_area_km2, d$threatened_area_km2, paired = TRUE)

# bland altman plot for supplement
d$mean_AOO <- (d$nonnative_area_km2 + d$threatened_area_km2) / 2
d$difference_AOO <- d$nonnative_area_km2 - d$threatened_area_km2

ggplot(d, aes(x = mean_AOO, y = difference_AOO)) +
  geom_point(aes(color = difference_AOO > 0), alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "#95846d") +
  scale_color_manual(values = c("TRUE" = "#1D6774", "FALSE" = "#FFD700")) + 
  labs(x = "Mean of naturalized and threatened range size (km²)", 
       y = "Difference (naturalized - threatened) range size (km²)") +
  theme_minimal(base_family = "Arial Narrow", base_size = 13) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "none"  # Remove legend if unnecessary
  )

  
showtext_opts(dpi=600)
ggsave("Figures/bland-altman.png", dpi = 600, bg = "white", height = 6, width = 6)
showtext_opts(dpi=96)

