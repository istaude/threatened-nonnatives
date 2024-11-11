source("R-code/00-preamble.R")

# load data ---------------------------------------------------------------
d <- read_csv("Data/non-native-red-list-species-NTexcluded.csv")

# phylogeny of non-natives that are also red-listed -----------------------

# what species are that? 
df <- d %>%
  group_by(taxon_name, wcvp_accepted_id) %>%
  summarize(status_unique = list(unique(status))) %>% 
  mutate(status_unique = ifelse(status_unique == 'c("threatened", "non-native")',
                                "threatenednonnative", "nonnative")
  ) %>% 
  # join back again higher level taxa
  left_join(wcvp_names %>% 
              select(plant_name_id, 
                     family, 
                     genus,
                     lifeform_description,
                     climate_description),
            by = c("wcvp_accepted_id" = "plant_name_id"))

# calculate for each family the fraction of threatened vs non-native species
family_threatened <- df %>% 
  group_by(family) %>% 
  count(status_unique) %>% 
  spread(status_unique, n) %>% 
  mutate_all(~replace_na(., 0)) %>% 
  mutate(sum = nonnative + threatenednonnative,
         perc = threatenednonnative/sum) %>% 
  arrange(desc(perc))

# what fraction of non-native families has threatened non-natives
1 - (
  family_threatened %>% filter(perc == 0) %>% nrow /
    family_threatened %>% nrow
)
# 74% of seed plant families containing non-native species had non-natives that
# were threatened in their native range.



# create family level phylogenetic tree -----------------------------------

head(family_threatened, n = 30)
data <- family_threatened %>% dplyr::select(family, perc) %>% 
  mutate(perc = ifelse(perc == 0, NA, perc))

# load phylogenetic tree from https://doi.org/10.1186/s12915-021-01166-2
file_path <- "Data/tree/4782_ML_433_12_Family_bootstrap.tre"  
phylo_tree <- read.nexus(file_path)
phylo_tree$edge.length <- NULL

# tip labels
tip_data_tree <- data.frame(label = phylo_tree$tip.label, 
                            family = str_replace(phylo_tree$tip.label, "^[^\\s]+\\s", ""))


# find the tips that correspond to the families in your data
tip_data_match <- inner_join(data, tip_data_tree)

# some families apparently wrong taxonomy try to find these manually
# ok these are all non-seed plants (mostly fern families)
# left_join(data, tree_data) %>% filter(is.na(label)) %>% select(fam = family)


# prune the tree to include only these tips
pruned_tree <- drop.tip(phylo_tree, setdiff(tip_data_tree$label, tip_data_match$label))
tree_data <- merge(tip_data_match, data, all.x = TRUE)

# join the tree with data
pruned_tree_data <- full_join(pruned_tree, tree_data, by = "label")

# plot the tree
# custom color gradient for the tree plot
color_gradient <- scale_color_gradientn(
  colours = c("#a3ffda", "#78c1b8", "#6aadad", "#3B7984", "#1E303A"),
  values = c(0.05, 0.25, 0.5, 0.75, 1),
  na.value = "grey90",
  labels = scales::label_percent()
)

# phylogenetic tree plot
(fig_phylo <- ggtree(pruned_tree_data, layout = 'circular', ladderize = FALSE) +
  geom_tiplab(aes(label = family, color = perc), size = 1.5, family = "Arial Narrow") +
  color_gradient +
  theme(legend.position = "right") +
  labs(color = "% threatened non-natives") +
  theme_minimal(base_family = "Arial Narrow", base_size = 13) +  # Ensure the base size is the same
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 10),  # Set consistent legend text size
    legend.margin = margin(t = 5, r = 5, b = 5, l = 5),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.5, 'cm'),
    axis.text = element_blank()
  ))


# calculate phylogenetic signal -------------------------------------------

# load tree again but do not set edge lengths to zero
phylo_tree <- read.nexus(file_path)
which(is.na(phylo_tree$edge.length))
# replace NA values with the mean of the existing edge lengths
mean_edge_length <- mean(phylo_tree$edge.length, na.rm = TRUE)
phylo_tree$edge.length[is.na(phylo_tree$edge.length)] <- mean_edge_length

# prune tree again to our data
pruned_tree <- drop.tip(phylo_tree, setdiff(tip_data_tree$label, tip_data_match$label))

# match percentage column with order in tree
row.names(tree_data) <- tree_data$label
tree_data <- tree_data[match(pruned_tree$tip.label, row.names(tree_data)), ]
trait_vector <- tree_data$perc
names(trait_vector) <- row.names(tree_data)

# calculate Pagel's lambda
pagel_lambda <- phylosig(pruned_tree, trait_vector, method = "lambda", test = TRUE)
print(pagel_lambda)

# calculate Blomberg’s K
blomberg_k <- phylosig(pruned_tree, trait_vector, method = "K", test = TRUE)
print(blomberg_k)

# prepare the data for the D-statistic
# there are some nodes that are not unique
node_labels <- pruned_tree$node.label
duplicated_nodes <- node_labels[duplicated(node_labels)]
print(duplicated_nodes)

# make unqiue
if (!is.null(pruned_tree$node.label)) {
  pruned_tree$node.label <- make.unique(pruned_tree$node.label)
}

# check
node_labels <- pruned_tree$node.label
print(any(duplicated(node_labels)))

# make percentage var binary
tree_data$perc_binary <- ifelse(is.na(tree_data$perc) == T, 0, 1)
tree_data <- tree_data %>% dplyr::select(label, perc_binary)

data_caper <- comparative.data(pruned_tree, tree_data, 
                               names.col = "label", vcv = TRUE, na.omit = T)

d_stat <- phylo.d(data_caper, binvar = perc_binary)
print(d_stat)

