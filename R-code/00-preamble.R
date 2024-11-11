library(tidyverse)
library(rWCVP)
library(rWCVPdata)
library(readxl)
library(patchwork)
library(paletteer)
library(ggExtra)
library(data.table)
library(conflicted)
library(countrycode)
library(rgbif)
library(data.table)
library(lwgeom)
library(mgcv)
library(lmodel2)

# map
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(MazamaSpatialUtils)

# phylogenetic analyses
library(ggtree) # BiocManager::install("ggtree")
library(ape)
library(phytools)
library(motmot)
library(caper)

# load fonts
library(hrbrthemes)
library(sysfonts)
library(showtext)
font_paths("C:/Windows/Fonts")
font_add("Arial Narrow", regular = "arialn.ttf", italic = "arialni.ttf", 
         bold = "arialnb.ttf")
showtext_auto()
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")

# own functions
# pretty labels
pretty_labels_percent <- function(x, digs = 2){
  
  a <- sapply(x, 
              
              FUN = function(x, digits = digs){
                
                a <- gsub("[^0-9.-]", " ", x)
                a <- trimws(a)
                a1 <- round(as.numeric(str_split(a, " ")[[1]][1]), digits)*100
                a2 <- round(as.numeric(str_split(a, " ")[[1]][2]), digits)*100
                y <- paste0(prettyNum(a1, big.mark=","), " to ", prettyNum(a2, big.mark=","), "%")
                
                y <- ifelse(is.na(x) ==  TRUE, "Missing", y)
                return(y)
              }
              
  )
  
  factor(a,  levels = unique(a[order(x)]))
  
}


pretty_labels_raw <- function(x, digs = 2){
  
  a <- sapply(x, 
              
              FUN = function(x, digits = digs){
                
                a <- gsub("[^0-9.-]", " ", x)
                a <- trimws(a)
                a1 <- round(as.numeric(str_split(a, " ")[[1]][1]), digits)
                a2 <- round(as.numeric(str_split(a, " ")[[1]][2]), digits)
                y <- paste0(prettyNum(a1, big.mark=","), " to ", prettyNum(a2, big.mark=","))
                
                y <- ifelse(is.na(x) ==  TRUE, "Missing", y)
                return(y)
              }
              
  )
  
  factor(a,  levels = unique(a[order(x)]))
  
}

# function to check if there are any entries within a 5-year window
within_5_years <- function(years) {
  n <- length(years)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (abs(years[i] - years[j]) <= 5) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

# load world map data
world <- ne_countries(scale = "medium", returnclass = "sf")
# print countries with missing iso_a3 codes
print(world %>% filter(iso_a3 == '-99') %>% select(name))
# fix for bug https://github.com/geopandas/geopandas/issues/1041
world <- world %>%
  mutate(
    iso_a3 = case_when(
      name == 'France' ~ 'FRA',
      name == 'Norway' ~ 'NOR',
      name == 'N. Cyprus' ~ 'CYP',
      name == 'Somaliland' ~ 'SOM',
      name == 'Kosovo' ~ 'RKS',
      TRUE ~ iso_a3
    )
  )

