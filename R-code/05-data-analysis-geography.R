source("R-code/00-preamble.R")

# load data ---------------------------------------------------------------
d <- read_csv("Data/non-native-red-list-species-NTexcluded.csv")

# make map that shows fraction of threatened non-natives in country -------

# we will use the iso 3 code, of course there are finer scales available
# as well, but not systematically. So the assumption here is if a plant
# was threatened at least once in time or in one part of a political country
# then this species is somewhere threatened in this political country (not 
# necessarily meaning the species is threatened everywhere in that country)

# first of all calculate the number of plant species that are non-native
# per country
nn <- d %>% filter(status == "non-native") %>% group_by(ISO3_Code) %>% count(status) %>% 
  arrange(desc(n))

# then find out how many of these species are threatened in their native range
df <- d %>%
  group_by(taxon_name, wcvp_accepted_id) %>%
  summarize(status_unique = list(unique(status))) %>% 
  mutate(status_unique = ifelse(status_unique == 'c("threatened", "non-native")',
                                "threatenednonnative", "nonnative")
  )


# now count only species that are non-native in a given country AND have the status_unique
# of threatenednonnative
nt <- df %>% 
  right_join(
    d
  ) %>% filter(status_unique == "threatenednonnative" & status == "non-native") %>% 
  group_by(ISO3_Code) %>% 
  count(status_unique) %>% arrange(desc(n)) %>% rename(threatened = n)


map_data <- full_join(nn %>% select(-status), nt %>% select(-status_unique)) %>% 
  mutate(threatened = ifelse(is.na(threatened)==T, 0, threatened)) %>% 
  mutate(ratio = threatened/n) %>% 
  arrange(desc(ratio))

# which iso3s have 1?
map_data %>% filter(ratio == 1)
# very low count of non-natives, makes it easy to have a high fraction

# what are the quantiles?
quantile(map_data$ratio, na.rm = T)

# how many more than half?
map_data %>% arrange(desc(ratio)) %>% View()
map_data %>% filter(ratio >= 0.5)
mean(map_data$ratio)
# for SI
# show histogram and cdf for the number of countries with a certain % threatened
# non-natives

# for scaling
histogram_data <- ggplot_build(
  ggplot(map_data %>% filter(n  > 20), aes(x = ratio)) +
    geom_histogram(binwidth = 0.05)
)$data[[1]]
max_count <- max(histogram_data$count)

ggplot(map_data %>% filter(n > 20), aes(x = ratio)) +
  geom_histogram(aes(y = ..count..), fill = "grey80", col = "white", binwidth = 0.05) + 
  stat_ecdf(aes(y = ..y.. * max_count), geom = "step", color = "grey30") + 
  scale_y_continuous(
    name = "Count", 
    sec.axis = sec_axis(~ . / max_count, name = "Cumulative density") 
  ) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "% threatened non-natives in country") +
  theme_minimal(base_family = "Arial Narrow", base_size = 13) +
  theme(panel.grid = element_blank())


showtext_opts(dpi=600)
ggsave("Figures/supp-histogram.png", dpi = 600, bg = "white", height = 6, width = 6)
showtext_opts(dpi=96)



# create accomp bar plot where you show the top countries, filtering countries
# that have at least 100 non-natives
map_data_pruned <- map_data %>% filter(n>=20) 
map_data_pruned$Country_Name <- countrycode(map_data_pruned$ISO3_Code, "iso3c", "country.name")
map_data_pruned$Country_Name[map_data_pruned$ISO3_Code == "ANT"] <- "Netherlands Antilles"
map_data_pruned$Continent <- countrycode(map_data_pruned$ISO3_Code, "iso3c", destination = "continent")
map_data_pruned$Continent[map_data_pruned$ISO3_Code == "ANT"] <- "Americas"       # Netherlands Antilles (historically classified here)
map_data_pruned$Continent[map_data_pruned$ISO3_Code == "ATF"] <- "Antarctica"          # French Southern Territories
map_data_pruned$Continent[map_data_pruned$ISO3_Code == "CCK"] <- "Asia"                # Cocos (Keeling) Islands
map_data_pruned$Continent[map_data_pruned$ISO3_Code == "IOT"] <- "Asia"                # British Indian Ocean Territory
map_data_pruned$Continent[map_data_pruned$ISO3_Code == "UMI"] <- "Oceania"             # United States Minor Outlying Islands
#View(map_data_pruned)

ggplot(map_data_pruned, aes(x = reorder(Country_Name, -ratio))) +
  facet_wrap(.~Continent, scale = "free_y")+
    geom_col(aes(y = n), width = 0.7, fill = "#78c1b8") +
    geom_point(aes(y = ratio * max(n), color = "Threatened Ratio")) +
    geom_hline(yintercept = 0.29* max(map_data_pruned$n), lty = "dashed", col = "grey30") +
    scale_y_continuous(
      name = "# non-natives",
      sec.axis = sec_axis(~ . / max(map_data_pruned$n),
                          labels = scales::percent, name = "% threatened")
    ) +
    scale_color_manual(name = "Threatened Ratio", values = c("Threatened Ratio" = "#FF6F61")) +
    labs(x = "") +
    theme_minimal(base_family = "Arial Narrow", base_size = 13) +
    theme(
      plot.margin = unit(c(0, 1, 0, 0), "cm"),
      axis.text.y = element_text(size = 7),
      legend.position = "none",
      panel.grid = element_blank()
    ) +
  coord_flip()



showtext_opts(dpi=600)
ggsave("Figures/supp-countrydata.png", dpi = 600, bg = "white", height = 8, width = 10)
showtext_opts(dpi=96)

# mapping -----------------------------------------------------------------

# world is loaded in the preamble
# merge the counts with the world map data
world_data <- world %>%
  left_join(map_data, by = c("iso_a3" = "ISO3_Code"))


# define the number of quantiles
n <- 5

# create quantiles and pretty labels for threatened species
world_data <- world_data %>%
  mutate(quantile_threatened = cut(ratio, 
                                   breaks = c(quantile(ratio, probs = seq(0, 1, by = 1/5), 
                                                       na.rm = T)),
                                   include.lowest = T, dig.lab=10)) %>% 
  mutate(pretty_labs_threatened = pretty_labels_percent(quantile_threatened))


# define color palette
cc <- c("#caffea", "#a3ffda", "#78c1b8", "#6aadad", "#3B7984")

# world map for threatened species
(fig_map <- ggplot() +
  geom_sf(data = world_data, aes(fill = quantile_threatened), size = 0.1, col = NA) +
  scale_fill_manual(values = cc, na.value = "grey90", 
                    name = "% threatened non-natives",
                    labels = levels(world_data$pretty_labs_threatened)) +
  coord_sf(crs = "+proj=robin") +
  theme_minimal(base_family = "Arial Narrow", base_size = 13) +  # Ensure the base size is the same
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 10),  # Set consistent legend text size
    legend.margin = margin(t = 5, r = 5, b = 5, l = 5),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.5, 'cm'),
    axis.text = element_blank()
  )
)


# create multipanel plot with phylogeny

layout <- "
A
B
"
combined_plot <- fig_phylo + fig_map + 
  plot_layout(design = layout) +
  plot_layout(heights = c(2, 1.4)) +
  plot_annotation(
    tag_levels = 'a',
    theme = theme(
      plot.tag = element_text(family = "Arial Narrow", face = "bold", size = 14),
      legend.title = element_text(family = "Arial Narrow", size = 14),
      legend.text = element_text(family = "Arial Narrow", size = 10)  # Ensure consistency
    )
  )

# display the combined plot
combined_plot

showtext_opts(dpi=600)
ggsave("Figures/map-phylo.png", dpi = 600, bg = "white", height = 8, width = 8)
showtext_opts(dpi=96)


# where are the threatened non-natives coming from ------------------------
d <- read_csv("Data/non-native-red-list-species-NTexcluded.csv")
threatened_where <- d %>% filter(status == "threatened") %>% 
  group_by(ISO3_Code) %>% 
  summarize(threatened_where = n_distinct(taxon_name))
nrow(threatened_where)
# smaller than RL number because in some countries none of the threatened species
# is non-native elswhere

# merge the counts with the world map data
world_data <- world %>%
  left_join(threatened_where, by = c("iso_a3" = "ISO3_Code"))


# define the number of quantiles
n <- 5

# create quantiles and pretty labels for threatened species
world_data <- world_data %>%
  mutate(quantile_threatened = cut(threatened_where, 
                                   breaks = c(quantile(threatened_where, 
                                                       probs = seq(0, 1, by = 1/5), 
                                                       na.rm = T)),
                                   include.lowest = T, dig.lab=10)) %>% 
  mutate(pretty_labs_threatened = pretty_labels_raw(quantile_threatened))


# define color palette
cc <- c("#caffea", "#a3ffda", "#78c1b8", "#6aadad", "#3B7984")

# world map for threatened species
(fig_map <- ggplot() +
    geom_sf(data = world_data, aes(fill = quantile_threatened), size = 0.1, col = NA) +
    scale_fill_manual(values = cc, na.value = "grey90", 
                      name = "% threatened non-natives",
                      labels = levels(world_data$pretty_labs_threatened)) +
    coord_sf(crs = "+proj=robin") +
    theme_minimal(base_family = "Arial Narrow", base_size = 13) +  # Ensure the base size is the same
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      panel.grid = element_blank(),
      legend.text = element_text(size = 10),  # Set consistent legend text size
      legend.margin = margin(t = 5, r = 5, b = 5, l = 5),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.spacing.y = unit(0.5, 'cm'),
      axis.text = element_blank()
    )
)

