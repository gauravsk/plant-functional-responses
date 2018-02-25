# A first pass at implementing some of the models I describe in "motivation.md"
# I won't be able to work with community data yet- so I focus on just the traitXenvironment effects on vital rates.
# I will work with just a small number of species for now- these are species that
# have fairly high performance (seed output in the absence of competitors) across the reserve.
# The species I will focus on for now are:
# PLER, HECO, SACO, LACA, CEME, HOMU, AGHE (==URLI), LOWR (== ACWR), AMME, and XXX. 

# I will work with just a small number of traits for now- namely, 
# I will focus on Seed Size, SLA, SRL, and flowering phenology.

# I will work with just a small number of environmental axes, as well- namely, 
# I will focus on soil Ca:Mg ratio, sand content, and soil moisture.

# Data import, munge, and setup ----------


library(tidyverse); theme_set(theme_bw())
library(lme4)
library(vegan)
library(glmmTMB)
library(bbmle)
thinkpad = TRUE

if(thinkpad){
  dat <- read_delim("data/performance/seed_production_processed.csv", 
                    col_names = T, delim = ",")
} else {
  dat <- read_delim("~/Dropbox/spatial_tapioca/data/performance/seed_production_processed.csv", 
                    col_names = T, delim = ",")
}
focal_species <- c("PLER", "LACA", "HECO", "HOMU", "CEME", "AGHE", "LOWR", "AMME", "SACO", "CHGL")
focal_species <- unique(dat$sp_code)

dat <- dat %>% filter(plot_type == "L") %>% 
  mutate(plot_num = paste0("plot_", plot_num),
         replicate = as.factor(replicate),
         num_seeds_produced = ifelse(is.na(num_seeds_produced), 0, num_seeds_produced)) %>%
  rename(species = sp_code, seed_production = num_seeds_produced, site = plot_num) %>% 
  select(-plot_type) %>%
  filter(species %in% focal_species)

dat
ggdat <- ggplot(dat)
gg_raw_seed <- ggplot(dat) + geom_boxplot(aes(site,seed_production+1)) + 
  facet_wrap(~species, ncol = 5) + scale_y_log10()
gg_raw_seed

gg_rawseed_hist <- ggplot(dat) + geom_histogram(aes(seed_production + 1))
gridExtra::grid.arrange(gg_rawseed_hist + ggtitle("Unlogged"), 
                        gg_rawseed_hist + scale_x_log10() + ggtitle("Logged"), ncol = 2) 

# And now we bring in the environmental data for each site
if(thinkpad){
  env_dat <- read_delim("data/environmental/all_environmental_data.csv", delim = ",")
} else {
  env_dat <- read_delim("~/Dropbox/spatial_tapioca/data/environmental/all_environmental_data.csv",
                        delim = ",") 
}
env_dat <- env_dat %>% rename(site = plot) %>% mutate(site = paste0("plot_", site))
# skimr::skim(env_dat) %>% select(-missing, -complete, -n)

plot(env_dat %>% select(-site, -lat, -lon, -type, - microsite, -ele), 
     pch = 21, bg = alpha("black", .25))

# There's a ton of variables, and so we can do an NMDS to wrap our heads around the dimensionality

env_dat_c <- env_dat %>% select(-type, -microsite, -lat, -lon, -Tmin, -ele) %>% 
  as.data.frame %>% tibble::column_to_rownames("site") 

soilmds <- metaMDS(env_dat_c)
par(mfrow = c(1,3))
plot(soilmds, main = "NMDS", xlim = c(-.25, .25))
orditorp(soilmds, display="species", cex = 1.2)
pc <-prcomp(env_dat_c %>% mutate_all(scale))
plot(pc)
biplot(pc, main = "PCA on scaled", cex = 1.3)

# OK, there's a lot of dimensionality, but let's focus on these:
# Ca:Mg ratio (data tells us it matters, and also we biologically know that it does)
# Soil texture (let's go with Sand for now)
# Nitrate (ppm), which goes strongly with organic_matter_ENR

env_vars_for_analysis <- env_dat_c %>% data.frame %>%
  tibble::rownames_to_column("site") %>%
  mutate(ca_mg = Ca_ppm/Mg_ppm) %>%
  select(site, ca_mg, nitrate_ppm = Nitrate_ppm, sand) %>% 
  mutate_if(is.numeric, scale) %>% 
  rename(ca_mg_scaled = ca_mg, nitrate_ppm_scaled = nitrate_ppm, sand_scaled = sand)

# Now we merge this env_vars_for_analysis df with the performance data frame `dat`

dat <- left_join(dat, env_vars_for_analysis)
tbl_df(dat)
ggdat <- ggplot(dat)
gg_camg <- ggdat + geom_point(aes(x = ca_mg_scaled, y = seed_production+1), alpha = 1/5) + facet_wrap(~species) + scale_y_log10()
gg_nitrate <- ggdat + geom_point(aes(x = nitrate_ppm_scaled, y = seed_production+1), alpha = 1/5) + facet_wrap(~species) + scale_y_log10()
gg_sand <- ggdat + geom_point(aes(x = sand_scaled, y = seed_production+1), alpha = 1/5) + facet_wrap(~species) + scale_y_log10()

gridExtra::grid.arrange(gg_camg,gg_nitrate, gg_sand, ncol = 3)

# Now we merge in the species traits
traits <- read_delim("data/trait/merged_sp_avg.csv", delim = ",") 
traits <- traits %>% mutate(species = ifelse(species == "SACA", "SACO", species))
traits <- traits %>% mutate(species = ifelse(species == "VUMA", "VUMI", species))

#traits <- traits %>% filter(species %in% unique(merged_df$species)) %>% arrange(species)
traits %>% tbl_df
traits$species
traits_sel <- traits %>% select(species, sla_cm2_g, max_height_cm, leaf_size_cm2) %>% 
  filter(species %in% unique(dat$species))

dat <- left_join(dat, traits_sel)
dat$seed_production_l <- log10(dat$seed_production+1)
a <- glmmTMB(seed_production_l ~ sla_cm2_g + ca_mg_scaled + sla_cm2_g:ca_mg_scaled + (1|species), data = dat)
b <- glmmTMB(seed_production_l ~ sla_cm2_g + sand_scaled + sla_cm2_g:sand_scaled + (1|species), data = dat)
c <- glmmTMB(seed_production_l ~ sla_cm2_g + nitrate_ppm_scaled + sla_cm2_g:nitrate_ppm_scaled + (1|species), data = dat)
summary(a)
summary(b)
summary(c)
anova(a)
dat2 <- dat %>% filter(seed_production > 0)

# Let's scale the seed production per species
dat_s <- dat %>% select(site, replicate, species, seed_production) %>% tidyr::spread(species, seed_production) %>% 
  mutate_if(is.numeric, scale_this) %>% 
  unite(id, c("site", "replicate"), sep = "XX") %>%
  tidyr::gather("species", "seed_production",2:18) %>% 
  separate(id, c("site", "replicate"), sep = "XX") %>% 
  left_join(., dat, by = c('site', 'species', 'replicate')) %>% select(-seed_production.y) %>% rename(seed_production = seed_production.x)

summary(a1)
