source("R-code/00-preamble.R")

# load data ---------------------------------------------------------------
d_final_all <- read_csv("Data/non-native-red-list-species-NTincluded.csv")
d_final_onlythreatened <- read_csv("Data/non-native-red-list-species-NTexcluded.csv")
d_countrysame <- read_csv("Data/non-native-red-list-species-countrysame.csv")


# inspect -----------------------------------------------------------------

# how many non-natives are also red-listed
d_final_onlythreatened %>%
  group_by(taxon_name) %>%
  summarize(status_combined = paste(sort(unique(status)), collapse = " and ")) %>% 
  count(status_combined)

# which fraction is that?
2513 / (6682 + 2513)
# 11306 non-native species (excluding hybrids, and accounted for newest WCVP taxonomy)
d %>% select(taxon_name) %>% distinct %>% nrow


# how many non-natives are also red-listed
d_final_all %>%
  group_by(taxon_name) %>%
  summarize(status_combined = paste(sort(unique(status)), collapse = " and ")) %>% 
  count(status_combined)

# which fraction is that?
2862 / (6333 + 2862)
# 11306 non-native species (excluding hybrids, and accounted for newest WCVP taxonomy)
d %>% select(taxon_name) %>% distinct %>% nrow



# pie/doughnut chart for viz
data <- data.frame(
  category = c("no", "yes"),
  value = c(72.7, 27.3)
)
data <- data %>%
  mutate(cumulative = c(65, 15),  # Position for labels
         label = paste0(value, "%"))  # Add percentage labels


# Create the chart
ggplot(data, aes(x = 2, y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +  
  coord_polar(theta = "y") +
  theme_void() +  
  xlim(0.5, 2.5) + 
  labs(fill = "") +  
  scale_fill_manual(values = c("#69acac", "#ff7062")) + 
  theme(legend.position = "none") + 
  geom_text(aes(x = 2, y = c(68, 12), label = c("No", "Yes")), 
            color = "white", size = 6, fontface = "bold", 
            family = "Arial Narrow", hjust = 0) +
  geom_text(aes(x = 2, y = c(65, 15), label = paste0(value, "%")), 
            color = "white", size = 5,
            family = "Arial Narrow") +
  annotate("text", x = 1.5, y = 0, 
           label = "Naturalizing species\nthreatened in parts of\ntheir native ranges", 
           color = "#06576D", size = 4.5,lineheight = 0.8,
           hjust = 0.5, vjust = 2,
           family = "Arial Narrow") +
  annotate("text", x = 1.5, y = 0, 
           label = "N = 9,195\nspecies assessed", 
           color = "grey40", size = 4, lineheight = 0.8,
           hjust = 0.5, vjust = 6,
           family = "Arial Narrow")

ggsave("Figures/pie_chart.svg", height = 5, width = 5)




# inspect now when countries are same -------------------------------------

# how many non-natives are also red-listed
d_countrysame %>%
  group_by(taxon_name) %>%
  summarize(status_combined = paste(sort(unique(status)), collapse = " and ")) %>% 
  count(status_combined)

# which fraction is that?
2365 / (6360 + 2365)
d_countrysame %>% select(taxon_name) %>% distinct %>% nrow



# overestimation? ---------------------------------------------------------

# test whether the range size of all species in wcvp predicts, listing
# on a national RL
drl <- read_csv("Data/drl-harmonized.csv")

# select species that are threatened somewhere
rl_specs <- drl %>% filter(status == "threatened") %>% 
  select(wcvp_accepted_id, status) %>% distinct

all_specs <- wcvp_names %>% 
  filter(taxon_status == "Accepted") %>% 
  filter(is.na(species_hybrid)) %>% 
  filter(!grepl("subsp\\.|var\\.", taxon_name)) %>% 
  mutate(
    work_genus = sapply(taxon_name, function(x) strsplit(x, " ")[[1]][1]),
    apomictic = ifelse(work_genus %in% apomixis_genus, 1, 0)
  ) %>% 
  filter(apomictic == 0) %>% 
  select(plant_name_id) %>% 
  distinct %>% 
  # anti_join the non-native taxa in our data 
  full_join(rl_specs,
            by = c("plant_name_id" = "wcvp_accepted_id") 
            
  )

nrow(all_specs)
nrow(rl_specs)

all_range <- all_specs %>% 
  left_join(
    wcvp_distributions %>% filter(introduced == 0)
  )  %>% 
  group_by(plant_name_id) %>% 
  summarise(range = n_distinct(area_code_l3))

summary(all_range)

all <- left_join(all_specs, all_range) %>% 
  mutate(status = ifelse(is.na(status)== T, 0, 1))

# raw data plot
ggplot(all, aes(x = range, y = status)) +
  geom_jitter(alpha = 0.05, height = 0.05, size = 0.1, col = "grey30") + 
  geom_smooth(col = "#FF6F61") +
  labs(
    x = "Range size (# botanical countries in which spp. is native)",
    y = "Probability of being threatened"
  ) +
  theme_minimal(base_family = "Arial Narrow", base_size = 13) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) -> figs2a



# bin data into range intervals and calculate the proportion of threatened species in each bin
binned_data <- all %>%
  mutate(range_bin = cut(range, breaks = 51)) %>% 
  group_by(range_bin) %>%
  summarize(
    median_range = median(range, na.rm = TRUE),
    threatened_prob = mean(status)
  )
nrow(binned_data)

# fit a GAM model
gam_model <- gam(threatened_prob ~ s(median_range, k = 10), data = binned_data)
summary(gam_model)
ggplot(binned_data, aes(x = median_range, y = threatened_prob)) +
  geom_point(col = "grey30") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10), col = "#FF6F61", se = T, size = 1) +
  labs(
    x = "Range size (binned)",
    y = "Probability of being threatened"
  ) +
  theme_minimal(base_family = "Arial Narrow", base_size = 13) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) -> figs2b

combined_plot <- 
  figs2a /
  figs2b +
  plot_annotation(
    tag_levels = 'a')

combined_plot

showtext_opts(dpi=600)
ggsave("Figures/probthreat_vs_range.png", dpi = 600, bg = "white", height = 8, width = 8)
showtext_opts(dpi=96)

