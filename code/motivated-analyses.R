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
thinkpad = F

# Bring in seed production data ------------
if(thinkpad){
  dat <- read_delim("data/performance/seed_production_processed.csv", 
                    col_names = T, delim = ",")
} else {
  dat <- read_delim("~/Dropbox/spatial_tapioca/data/performance/seed_production_processed.csv", 
                    col_names = T, delim = ",")
}
# focal_species <- c("PLER", "LACA", "HECO", "HOMU", "CEME", "AGHE", "LOWR", "AMME", "SACO", "CHGL")
focal_species <- unique(dat$sp_code)
# focal_plots <- 740:755 # Toggle to select Candy Valley only
focal_plots <- 740:763
min_seed <- 0
dat <- dat %>% # filter(plot_type == "L") %>% 
  filter(plot_num %in% focal_plots) %>%
  mutate(plot_num = paste0("plot_", plot_num),
         replicate = as.factor(replicate),
         num_seeds_produced = ifelse(is.na(num_seeds_produced), 0, num_seeds_produced)) %>%
  filter(num_seeds_produced >= min_seed) %>%
  rename(species = sp_code, seed_production = num_seeds_produced, site = plot_num) %>% 
  filter(species %in% focal_species)


# Read in the ML estimates of lambda and "alpha" (i.e. sensitivity to competitors) ----------
l_and_e_estimates <- read_delim("model_outputs/mod5pars.csv", delim = ",") %>% rename(species = X1)
l_and_e_estimates <- l_and_e_estimates %>% select(-sigma,-convergence_code) %>% 
  gather("type", "value", 2:49) %>% 
  separate(type, sep = "_", into = c("delme","plot","type")) %>% unite(site, delme, plot, remove = T) %>%
  spread("type", "value") 
l_and_e_estimates
# Merge in with `dat`
dat <- left_join(dat, l_and_e_estimates)

# Read in germination data ----------
germ <- read_delim("~/Dropbox/spatial_tapioca/data/performance/calculated/plot_level_performance_summary.csv", delim = ",")
germ <- germ %>% select(site = plot_num, plot_type, species = sp_code, mean_germ_percentage = germ_percentage) %>%
  mutate(site = paste0("plot_", site))
# Merge in with `dat`
dat <- left_join(dat, germ)

# Dat now has the following columns:
# Site, replicate, species, seed_production, alpha, lambda, mean_germ_percentage

dat

ggdat <- ggplot(dat)
gg_raw_seed <- ggplot(dat) + geom_boxplot(aes(site,seed_production+1)) + 
  facet_wrap(~species, ncol = 5) + scale_y_log10()
X11(); gg_raw_seed

gg_rawseed_hist <- ggplot(dat) + geom_histogram(aes(seed_production + 1))
X11(); gridExtra::grid.arrange(gg_rawseed_hist + ggtitle("Unlogged seed production"), 
                        gg_rawseed_hist + scale_x_log10() + ggtitle("Logged seed production"), ncol = 2) 

dat_l <- dat %>% filter(plot_type == "L")
# Print out the mean and the sd
dat_l %>% select(site, replicate, species, seed_production) %>% 
  tidyr::spread(species, seed_production) %>% summarize_if(is.numeric, funs(sd(.,na.rm = T)))
dat_l %>% select(site, replicate, species, seed_production) %>% 
  tidyr::spread(species, seed_production) %>% summarize_if(is.numeric, funs(sd(.,na.rm = T)))




# And now we bring in the environmental data for each site -----------
if(thinkpad){
  env_dat <- read_delim("data/environmental/all_environmental_data.csv", delim = ",")
} else {
  env_dat <- read_delim("~/Dropbox/spatial_tapioca/data/environmental/all_environmental_data.csv",
                        delim = ",") 
}
env_dat <- env_dat %>% rename(site = plot) %>% mutate(site = paste0("plot_", site))
# skimr::skim(env_dat) %>% select(-missing, -complete, -n)

X11(); plot(env_dat %>% filter(!(site %in% c("plot_758", "plot_759"))) %>% select(-site, -lat, -lon, -type, - microsite, -ele), 
     pch = 21, bg = alpha("black", .25))

# There's a ton of variables, and so we can do an NMDS to wrap our heads around the dimensionality

env_dat_c <- env_dat %>% select(-type, -microsite, -lat, -lon, -Tmin, -ele) %>% 
  as.data.frame %>% tibble::column_to_rownames("site") 

soilmds <- metaMDS(env_dat_c)
X11()
par(mfrow = c(1,3))
plot(soilmds, main = "NMDS", xlim = c(-.25, .25))
orditorp(soilmds, display="species", cex = 1.2)
pc <-prcomp(env_dat_c %>% mutate_all(scale))
plot(pc)
biplot(pc, main = "PCA on scaled", cex = 1.3)

# Merge in NMDS scores into env_dat
env_dat <- left_join(env_dat, soilmds$points %>% data.frame %>% tibble::rownames_to_column("site"))

# Now we merge this env_vars_for_analysis df with the performance data frame `dat`
dat <- left_join(dat, env_dat)

dat_plot_yesno <- dat %>% group_by(site, plot_type,species) %>% mutate(seed_production = if_else(seed_production))

ggdat <- ggplot(dat)

gg_mds1 <- ggdat + geom_point(aes(x = MDS1, y = seed_production+1), alpha = 1/5) + facet_wrap(~species) + scale_y_log10()
gg_mds2 <- ggdat + geom_point(aes(x = MDS2, y = seed_production+1), alpha = 1/5) + facet_wrap(~species) + scale_y_log10()

gridExtra::grid.arrange(gg_camg,gg_nitrate, gg_sand, ncol = 3)

# Now we merge in the species traits ----------
if(thinkpad){
  traits <- read_delim("data/trait/merged_sp_avg.csv", delim = ",")
} else {
  traits <- read_delim("~/Dropbox/spatial_tapioca/data/trait/merged_sp_avg.csv",
                        delim = ",") 
}
traits <- traits %>% mutate(species = ifelse(species == "SACA", "SACO", species))
traits <- traits %>% mutate(species = ifelse(species == "VUMA", "VUMI", species))

#traits <- traits %>% filter(species %in% unique(merged_df$species)) %>% arrange(species)
traits %>% tbl_df
traits$species
traits_sel <- traits %>% select(species, sla_cm2_g, max_height_cm, leaf_size_cm2) %>% 
  filter(species %in% unique(dat$species))

dat <- left_join(dat, traits_sel)


# Do some more data prep ------------
# Make a new column seed_production_l that has log-transformed seed production
dat$seed_production_l <- log10(dat$seed_production+1)

# Make a column seed_production_scaled that includes seed production, scaled per species
scale_this <- function(x) as.vector(scale(x))
dat <- dat %>% select(site, replicate, species, seed_production) %>% tidyr::spread(species, seed_production) %>% 
  mutate_if(is.numeric, scale_this) %>% 
  unite(id, c("site", "replicate"), sep = "XX") %>%
  tidyr::gather("species", "seed_production",2:18) %>% 
  separate(id, c("site", "replicate"), sep = "XX") %>% 
  left_join(., dat, by = c('site', 'species', 'replicate')) %>% rename(seed_production = seed_production.y,
                                                                       seed_production_scaled  = seed_production.x)







# Run some models --------------
fit_zipoisson <- glmmTMB(seed_production~ca_mg_scaled + nitrate_ppm_scaled + sand_scaled + sla_cm2_g + 
                           sla_cm2_g*ca_mg_scaled + sla_cm2_g*nitrate_ppm_scaled + sla_cm2_g*sand_scaled + 
                           (seed_production|species) + (1|site),
                         ziformula = ~1, family = poisson, data = dat)

fit_zipoisson2 <- update(fit_zipoisson, ziformula = ~ species)
fit_nbinom <- update(fit_zipoisson2, family = nbinom2)
# fit_nbinom2 <- update(fit_nbinom, family = nbinom1) # This model has a non-positive-definite Hessian matrix and so is excluded

AICtab(fit_zipoisson, fit_zipoisson2, fit_nbinom)

# We try a hurdle-model, using fit_nbinom as the starting point since it had the lowest AIC

# fit_hnbinom1 <-  update(fit_nbinom,
#                         ziformula=~.,
#                         family=list(family="truncated_nbinom1",link="log"))
# But we exclude it as well, since it has a non-positive-definite Hessian matrix

# AICtab still suggests fit_nbinom as the optimal, but I want to leave out the terms

AICtab(fit_zipoisson, fit_zipoisson2, fit_nbinom, fit_nbinom_norandoms, fit_nbinom_notrait, fit_nbinom_nocamg)

# What if we go back and add in terms one by one
fit_nbinom_camgOnly <- glmmTMB(seed_production ~ ca_mg_scaled, family = nbinom2, data = dat)
fit_nbinom_nitrateOnly <- glmmTMB(seed_production ~ nitrate_ppm_scaled, family = nbinom2, data = dat)
fit_nbinom_sandOnly <- glmmTMB(seed_production ~ sand_scaled, family = nbinom2, data = dat)

AICtab(fit_zipoisson, fit_zipoisson2, fit_nbinom, fit_nbinom_norandoms, fit_nbinom_notrait, fit_nbinom_nocamg,
       fit_nbinom_sandOnly, fit_nbinom_nitrateOnly, fit_nbinom_camgOnly)
 

# ---------------
# Bring in estimates of L and r/alpha


l_and_e_estimates <- read_delim("model_outputs/mod5pars.csv", delim = ",") %>% rename(species = X1)
l_and_e_estimates <- l_and_e_estimates %>% select(-sigma,-convergence_code) %>% 
  gather("type", "value", 2:49) %>% 
  separate(type, sep = "_", into = c("delme","plot","type")) %>% unite(site, delme, plot, remove = T) %>%
  spread("type", "value") 
l_and_e_estimates

ggplot(l_and_e_estimates) + 
  geom_point(aes(x = lambda, y = alpha)) + 
  facet_wrap(~species, scales = "free")
l_and_e_estimates %>% filter(species == "HECO") %>% View

env_dat_c <- left_join(env_dat_c %>% tibble::rownames_to_column("site"), 
                       soilmds$points %>% data.frame %>% tibble::rownames_to_column("site"))
# OK, there's a lot of dimensionality, but let's focus on these:
# Ca:Mg ratio (data tells us it matters, and also we biologically know that it does)
# Soil texture (let's go with Sand for now)
# Nitrate (ppm), which goes strongly with organic_matter_ENR

env_vars_for_analysis <- env_dat_c %>% data.frame %>%
  mutate(ca_mg = Ca_ppm/Mg_ppm) %>%
  select(site, ca_mg, nitrate_ppm = Nitrate_ppm, sand, MDS1, MDS2, depth) %>% 
  mutate_if(is.numeric, scale) %>% 
  rename(ca_mg_scaled = ca_mg, nitrate_ppm_scaled = nitrate_ppm, sand_scaled = sand)

joined <- left_join(l_and_e_estimates, env_vars_for_analysis)

# Add a "scaled_lambda" column to joined that has the lambda value, scaled by species
joined <- left_join(joined, candi_v %>% select(species, site, lambda) %>% 
                      tidyr::spread(., species, lambda) %>% mutate_if(is.numeric, scale) %>% 
                      tidyr::gather(.,species, scaled_lambda, 2:18))

ggplot(joined, aes(x = MDS1, y = scaled_lambda)) + geom_boxplot() + 
  geom_point(color = alpha("black", .5), size = 2.5)

plots_to_keep <- paste0("plot_", 740:763)

candi_v <- joined %>% filter(site %in% plots_to_keep)
ggplot(candi_v) + geom_point(aes(x = ca_mg_scaled, y = lambda)) + facet_wrap(~species, scales = "free_y")
ggplot(candi_v) + geom_point(aes(x = nitrate_ppm_scaled, y = lambda)) + facet_wrap(~species, scales = "free_y")
ggplot(candi_v) + geom_point(aes(x = sand_scaled, y = lambda)) + facet_wrap(~species, scales = "free_y")


X11();ggplot(joined) + geom_point(aes(x = depth, y = scaled_lambda)) + facet_wrap(~species, scales = "free_y") + 
  geom_smooth(aes(x = depth, y = scaled_lambda), method = "lm")
X11();ggplot(joined) + geom_point(aes(x = depth, y = lambda)) + facet_wrap(~species, scales = "free_y") + scale_y_log10() +
  geom_smooth(aes(x = depth, y = lambda), method = "lm")

X11();ggplot(joined) + geom_point(aes(x = MDS1, y = scaled_lambda)) + facet_wrap(~species, scales = "free_y") + 
  geom_smooth(aes(x = MDS1, y = scaled_lambda), method = "lm")
X11();ggplot(joined) + geom_point(aes(x = MDS1, y = lambda)) + facet_wrap(~species, scales = "free_y") + scale_y_log10() +
  geom_smooth(aes(x = MDS1, y = lambda), method = "lm")
