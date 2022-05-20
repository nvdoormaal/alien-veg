######################
## ALIEN VEGETATION ##
######################

#Empty environment
rm(list=ls())

#Load packages
pacman::p_load(tidyverse, ggplot2, sf, paletteer, lubridate, ggspatial)
               #readxl, raster, adehabitatHR)

#Load OWNR shapfile
OWNR.WGS <- st_read(dsn = "./Data/OWNRBufferzone.shp")  ## Shapefile including buffer zone
OWNR.UTM <- st_transform(OWNR.WGS, crs = 32736)         ## Convert shp to UTM36

#Load Ekuthuleni shapefile
Eku.UTM <- st_read(dsn = "./Data/Ekuthuleni.shp")

#Buffer
OWNR.buf <- st_buffer(OWNR.UTM, dist = 500)

#Load data
AllData <- read_csv("./Data/AlienVeg_2020_2021.csv")

CleaningData <- AllData %>% 
  mutate(
    Species = str_to_lower(Species),
    Biocontrol = str_to_lower(Biocontrol),
    `Chemical c` = str_to_lower(`Chemical c`),
    Notes = str_to_lower(Notes),
    timestamp = as_datetime(timestamp),
    NewQuarter = quarter(timestamp, with_year = TRUE),
    NewQuarter = str_replace(NewQuarter, pattern = "[.]", replacement = " - Quarter "),
    Quarter = coalesce(Quarter, NewQuarter),
    NewControl = case_when(
      str_detect(Notes, pattern = "dead") | str_detect(Notes, pattern = "dying") ~ "Dead",
      str_detect(Notes, pattern = "nurse") ~ "Removed",
      str_detect(Biocontrol, pattern = "bug") & str_detect(`Chemical c`, pattern = "msma") ~ "Both",
      str_detect(Biocontrol, pattern = "bug") ~ "Bio",
      str_detect(`Chemical c`, pattern = "msma") ~ "Chemical",
      TRUE ~ "Untreated"
    )
  ) %>% 
  unique() %>% 
  subset(!(Species == "uknown")) %>% 
  select(c(Species, Quarter, NewControl, Latitude, Longitude, Northing, Easting))

##Boxing glove cactus dataset
BoxGlove <- read_csv("./Data/BGC_2021.csv")

BoxGlove <- select(BoxGlove, c("Species", "Quarter", "NewControl", "Latitude", "Longitude", "Northing", "Easting"))


##Final dataset
AlienVeg <- rbind(
  CleaningData[str_which(CleaningData$Quarter, pattern = "2019", negate = TRUE),],
  BoxGlove
)

AlienVeg <- AlienVeg %>% 
  mutate(
    NewSpecies = case_when(
      c(Species == "boxing glove cactus" | Species == "jointed cactus" | Species == "moon cactus") ~ "other",
      TRUE ~ Species
      )
  )

##Create maps
# Categorised by species and year-quarter
QuarterMap <- ggplot() +
  geom_sf(data = OWNR.UTM, fill = "#90ee90", alpha = 0.2, color = "black") +
  geom_sf(data = Eku.UTM, fill = "#90ee90", alpha = 0.2, color = "black") + 
  geom_point(data = AlienVeg, aes(x = Easting, y = Northing, shape = NewSpecies, color = NewControl), size = 3) +
  scale_shape("Species") + 
  scale_color_paletteer_d(name = "Treatment", "ggthemes::colorblind") + 
  facet_wrap(~Quarter) + 
  ggtitle("Locations and treatment of alien vegetation species on Olifants West per quarter",
          subtitle = "Other species include boxing glove, jointed, and moon cactus") +
  annotation_north_arrow(which_north = "grid", location = "tr", height = unit(1, "cm"), width = unit(1, "cm")) +
  annotation_scale(location = "bl") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

# ggsave(filename = paste0("AlienVegQuarters_", today(), ".svg"), plot = QuarterMap, device = grDevices::svg, path = "./Figures",
#        height = unit("6", "cm"), width = unit(10, "cm"))

ggsave(filename = paste0("AlienVegQuarters_", today(), ".png"), 
       plot = QuarterMap, height = unit("6", "cm"), width = unit(10, "cm"), 
       device = "png", path = "./Figures", dpi = 500)
                         

##################33



#Select only prickly pears and queen of the night
DataPP <- AllData %>% 
  subset(c(Species == "Prickly pear" | Species == "Queen of the night")) %>%
  drop_na(c(X, Y)) %>%
  st_as_sf(coords = c(x= "X", y = "Y"), crs = 4326) %>%
  st_transform(crs = 32736) %>%
  group_by(Species) %>% 
  group_split()

sf_to_KDE <- function(data, rangevalues){

  data.sp <- data %>% 
    st_geometry() %>%
    as_Spatial()
  
  data.kde <- kernelUD(data.sp, h=500, grid=1000, extent=0.5)
  data.kde <- getvolumeUD(data.kde)

  data.r <- raster(as(data.kde,"SpatialPixelsDataFrame"))
  data.r[data.r>90] <- NA
  data.r <- trim(data.r)
  data.rp <- rasterToPoints(data.r)
  data.rp <- data.frame(data.rp)
  
  colnames(data.rp) <- c("X", "Y", "Kernel")
  data.rp <- add_column(data.rp, Species = data$Species[1])
}

B <- purrr::map(DataPP, sf_to_KDE)

C <- bind_rows(B)

ggplot(data = C, aes(x=X, y=Y, fill=Kernel)) +
  geom_raster() +
    geom_sf(data = OWNR.UTM, inherit.aes = FALSE) +
  facet_wrap(~Species) +
  geom_point(alpha =0.5)

#######################################################################################
### ALIEN VEG ######
### 31 AUG 2021 ####
####################

## Load packages
pacman::p_load(tidyverse, ggplot2, sf, paletteer, lubridate, ggspatial, zoo, grid, gridExtra)
  
#Load OWNR shapfile
OWNR.WGS <- st_read(dsn = "./Data/OWNRBufferzone.shp")  ## Shapefile including buffer zone
OWNR.UTM <- st_transform(OWNR.WGS, crs = 32736)         ## Convert shp to UTM36

#Load Ekuthuleni shapefile
Eku.UTM <- st_read(dsn = "./Data/Ekuthuleni.shp")

#Buffer
OWNR.buf <- st_buffer(OWNR.UTM, dist = 500)

#Load data
Andy_raw <- read_csv("./Data/Central region Andy.csv")
Izel_raw <- read_csv("./Data/Central region Izel.csv")

##Cleaning attempt - Andy
Andy_clean <- Andy_raw %>% 
  filter(!is.na(Title)) %>% 
  filter(!is.na(`Date Created`)) %>% 
  filter(!is.na(Longitude) & !is.na(Latitude)) %>% 
  select(1:9)

##Cleaning attempt - Izel
Izel_clean <- Izel_raw %>% 
  filter(!is.na(Title)) %>% 
  filter(!is.na(`Date Created`)) %>% 
  filter(!is.na(Longitude) & !is.na(Latitude)) %>% 
  select(1:8)

##Combine data
AllData <- bind_rows(Andy_clean, Izel_clean)

##Keep Jan - Jun
AllData <- subset(AllData, as.Date(`Date Created`) < as.Date("2021-07-01"))

##Final Tree wrapping dataset
AllData_wrap <- AllData %>% 
  mutate(
    Title = str_to_lower(Title)) %>% 
  filter(
    c(str_detect(Title, pattern = "mesh") | str_detect(Title, pattern = "wrap"))
  )

Wrap_final <- AllData_wrap %>% 
  mutate(
    Species = case_when(
      str_detect(AllData_wrap$Title, pattern = "arula") | str_detect(AllData_wrap$Title, pattern = "mar") ~ "Marula",
      str_detect(AllData_wrap$Title, pattern = "thorn") ~ "Knobthorn",
      TRUE ~ "Other"
    ),
    YearQtr = as.yearqtr(as_date(`Date Created`))
  )

Wrap.sf <- Wrap_final %>% 
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(crs = 32736)

##Final Alien Veg dataset
AllData_Alien <- AllData %>% 
  unite('Merged_X', 2:6, remove = TRUE) %>% 
  mutate(
    Title = str_to_lower(Title),
    Merged_X = str_to_lower(Merged_X)) %>% 
  filter(
    c(str_detect(Title, pattern = "^pp") | str_detect(Title, pattern = "pp car") | 
        str_detect(Title, pattern = "qon") | str_detect(Title, pattern = "queen") | 
        str_detect(Title, pattern = "coc")) | str_detect(Title, pattern = "box")
  )

Alien_final <- AllData_Alien %>% 
  mutate(
    Species = case_when(
      str_detect(Title, pattern = "^pp") | str_detect(Title, pattern = "pp car") ~ "Prickly pear",
      str_detect(Title, pattern = "qon") | str_detect(Title, pattern = "queen") ~ "Queen of the night",
      str_detect(Title, pattern = "coc") ~ "Cocklebur",
      str_detect(Title, pattern = "box") ~ "Boxing glove",
      TRUE ~ "Other"
    ),
    Treatment = case_when(
      str_detect(Title, pattern = "bug") | str_detect(Merged_X, pattern = "bc") | str_detect(Merged_X, pattern = "bug") ~ "Treated",
      str_detect(Merged_X, pattern = "no bugs") ~ "Untreated",
      str_detect(Merged_X, pattern = "dying") | str_detect(Merged_X, pattern = "nursery") ~ "Dead/Removed",
      TRUE ~ "Untreated"
    ),
    YearQtr = as.yearqtr(as_date(`Date Created`))
  )

Alien.sf <- Alien_final %>% 
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(crs = 32736)

##Create maps
# Tree wrapping maps by quarter and species
WrapMap <- ggplot() +
  geom_sf(data = OWNR.UTM, fill = "#90ee90", alpha = 0.2, color = "black") +
  geom_sf(data = Eku.UTM, fill = "#90ee90", alpha = 0.2, color = "black") + 
  geom_sf(data = Wrap.sf, aes(color = Species), size = 2) +
  scale_color_paletteer_d(name = "Species", "ggthemes::wsj_red_green") + 
  facet_grid(as.factor(YearQtr)~Species) + 
  ggtitle("Locations of tree wrapping on Olifants West in 2021",
          subtitle = "Categorised by species and quarter") +
  annotation_north_arrow(which_north = "grid", location = "tr", height = unit(1, "cm"), width = unit(1, "cm")) +
  annotation_scale(location = "bl") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
  )

# ggsave(filename = paste0("WrapMapQtr2021_", today(), ".svg"), plot = WrapMap,
#        device = grDevices::svg, path = "./Figures",
#        height = unit("6", "cm"), width = unit(10, "cm"))

##Summary graph
Wrap_sum <- Wrap_final %>% 
  count(YearQtr, Species) %>% 
  mutate(
    YearQtr = as.factor(YearQtr)
  )

 Wrap_plot <- ggplot(data = Wrap_sum) +
    geom_bar(aes(x=YearQtr, y = n, fill = Species), stat = "identity", position = position_dodge2(width = 0.9, preserve = "single")) +
    scale_fill_paletteer_d(name = "Species", "ggthemes::wsj_red_green") +
    ylab("Number of trees wrapped") +
    xlab("Time (in quarters)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    text = element_text(size = 15)
  )

 ## Create table
 Wrap_sum <- Wrap_sum %>% 
   arrange(YearQtr, Species) %>% 
   rename("Quarter" = YearQtr, "Number" = n)
 
 Wrap_table <- tableGrob(Wrap_sum, rows = NULL, theme = ttheme_minimal())
 
 ## Combine table and graph
Wrap_plottable <- Wrap_plot + 
  annotation_custom(Wrap_table, xmin = -1.2, ymax=170)
  

# ggsave(filename = paste0("WrapPlotQtr2021_", today(), ".svg"), plot = Wrap_plottable,
#        device = grDevices::svg, path = "./Figures",
#        height = unit("6", "cm"), width = unit(10, "cm"))
 
 
 ##Create maps
 # Alien Veg maps by quarter and species
 AlienMap <- ggplot(data = Alien.sf) +
   geom_sf(data = OWNR.UTM, fill = "#90ee90", alpha = 0.2, color = "black") +
   geom_sf(data = Eku.UTM, fill = "#90ee90", alpha = 0.2, color = "black") + 
   geom_sf(aes(color = Species, shape = Treatment), size = 2) +
   scale_shape_discrete(solid = TRUE) +
   ggthemes::scale_color_colorblind() + 
   facet_grid(Species~as.factor(YearQtr)) + 
   ggtitle("Locations of alien plants on Olifants West in 2021",
           subtitle = "Categorised by species and quarter") +
   annotation_north_arrow(which_north = "grid", location = "tr", height = unit(0.5, "cm"), width = unit(0.5, "cm")) +
   annotation_scale(location = "bl") +
   theme_bw() +
   theme(
     axis.title = element_blank(),
     axis.text = element_blank(),
     axis.ticks = element_blank(),
     panel.grid = element_blank(),
     legend.position = "bottom",
     legend.title = element_text(face = "bold")
   )
 
 # ggsave(filename = paste0("AlienMapQtr2021_", today(), ".svg"), plot = AlienMap,
 #        device = grDevices::svg, path = "./Figures",
 #        height = unit("6", "cm"), width = unit(10, "cm"))
 
 ##Summary graph
Alien_sum <- Alien_final %>% 
   count(YearQtr, Species, Treatment) %>% 
   mutate(
     YearQtr = as.factor(YearQtr)
   )
 
 Alien_plot <- ggplot(data = Alien_sum) +
   geom_bar(aes(x=YearQtr, y = n, fill = Species), color="black", stat = "identity", position = position_dodge2(width = 0.9, preserve = "single")) +
   ggthemes::scale_fill_colorblind() +
   facet_wrap(~Treatment) +
   ylab("Number of alien plants") +
   xlab("Time (in quarters)") +
   theme_bw() +
   theme(
     legend.position = "bottom",
     legend.title = element_text(face = "bold"),
     panel.grid.major.x = element_blank()
   )

 ## Create table
Alien_sum <- Alien_sum %>% 
   arrange(YearQtr, Species) %>% 
   rename("Quarter" = YearQtr, "Number" = n)
 
Alien_table <- tableGrob(Alien_sum, rows = NULL, theme = ttheme_minimal())
 
 ## Combine table and graph
Alien_plottable <- grid.arrange(Alien_plot, Alien_table, ncol=2)

# ggsave(filename = paste0("AlienPlotQtr2021_", today(), ".svg"), plot = Alien_plottable,
#        device = grDevices::svg, path = "./Figures",
#        height = unit("6", "cm"), width = unit(10, "cm"))
