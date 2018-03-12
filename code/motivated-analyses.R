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
thinkpad = T

# Bring in seed production data ------------
if(thinkpad){
  dat <- read_delim("data/performance/seed_production_processed.csv", 
                    col_names = T, delim = ",")
} else {
  dat <- read_delim("~/Dropbox/spatial_tapioca/data/performance/seed_production_processed.csv", 
                    col_names = T, delim = ",")
}
# focal_species <- c("PLER", "LACA", "HECO", "HOMU", "CEME", "AGHE", "LOWR", "AMME", "SACO", "CHGL")
# focal_species <- unique(dat$sp_code)
# focal_plots <- 740:755 # Toggle to select Candy Valley only
focal_plots <- 740:763
min_seed <- 0
dat <- dat %>% # filter(plot_type == "L") %>% 
  filter(plot_num %in% focal_plots) %>%
  mutate(plot_num = paste0("plot_", plot_num),
         replicate = as.factor(replicate),
         num_seeds_produced = ifelse(is.na(num_seeds_produced), 0, num_seeds_produced),
         sp_code = ifelse(sp_code == "VUMI", "VUMA", sp_code)) %>%
  filter(num_seeds_produced >= min_seed) %>%
  rename(species = sp_code, seed_production = num_seeds_produced, site = plot_num) 


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
if(thinkpad){
  germ <- read_delim("data/performance/calculated/plot_level_performance_summary.csv", 
                    col_names = T, delim = ",")
} else {
  germ <- read_delim("~/Dropbox/spatial_tapioca/data/performance/calculated/plot_level_performance_summary.csv", 
                    col_names = T, delim = ",")
}
germ <- germ %>% select(site = plot_num, plot_type, species = sp_code, 
                        mean_germ_percentage = germ_percentage) %>%
  mutate(site = paste0("plot_", site),
         species = ifelse(species == "VUMI", "VUMA", species))
# Merge in with `dat`
dat <- left_join(dat, germ)

# Dat now has the following columns:
# Site, replicate, species, seed_production, alpha, lambda, mean_germ_percentage
dat

# Let's make a dat_sum that has these numbers summarized at the Sp*site*plot-type level
dat_sum <- dat %>% group_by(site, species) %>% summarize_if(is.numeric, funs(mean(., na.rm = T))) %>%
  arrange(site,species) %>% select(-seed_production)


ggdat <- ggplot(dat_sum)
gg_raw_seed <- ggdat + geom_boxplot(aes(site,lambda))
X11(); gg_raw_seed

gg_lambda_hist <- ggdat + geom_histogram(aes(lambda))
X11(); gridExtra::grid.arrange(gg_lambda_hist + ggtitle("Unlogged estimates of Lambda (sp x site)"), 
                               gg_lambda_hist + scale_x_log10() + ggtitle("Logged estimates of Lambda (sp x site)"), ncol = 2) 
dat_sum$log_lambda <- log10(dat_sum$lambda)

# And now we bring in the environmental data for each site -----------
if(thinkpad){
  env_dat <- read_delim("data/environmental/all_environmental_data.csv", delim = ",")
} else {
  env_dat <- read_delim("~/Dropbox/spatial_tapioca/data/environmental/all_environmental_data.csv",
                        delim = ",") 
}
env_dat <- env_dat %>% rename(site = plot) %>% mutate(site = paste0("plot_", site))

# skimr::skim(env_dat) %>% select(-missing, -complete, -n)

X11(); plot(env_dat %>% select(-site, -lat, -lon, -type, - microsite, -ele), 
     pch = 21, bg = alpha("black", .25))

# There's a ton of variables, and so we can do an NMDS to wrap our heads around the dimensionality

env_dat_c <- env_dat %>% select(-type, -microsite, -lat, -lon, -Tmin, -ele) %>% 
  as.data.frame %>% tibble::column_to_rownames("site") 

soilmds <- metaMDS(env_dat_c)
# X11()
# par(mfrow = c(1,3))
# plot(soilmds, main = "NMDS", xlim = c(-.25, .25))
# orditorp(soilmds, display="species", cex = 1.2, col = "red")
# orditorp(soilmds, display="sites", cex = 1, labels = as.character(740:763))
# 
# pc <-prcomp(env_dat_c %>% mutate_all(scale))
# plot(pc)
# biplot(pc, main = "PCA on scaled", cex = 1.3)
soilmds_envVars <- data.frame(soilmds$species)
soilmds_sites <- data.frame(soilmds$points)
ggplot(soilmds_envVars) + geom_point(aes(x = MDS1, y = MDS2), col = "darkred", size = 2) + xlim(c(-.25,.25)) +
  geom_text_repel(aes(MDS1, MDS2, label = rownames(soilmds_envVars)), col = "darkred", size = 5.2) +
  geom_point(data = soilmds_sites, aes(MDS1, MDS2), size = .75, col = "grey25") +
  geom_text_repel(data = soilmds_sites, aes(MDS1, MDS2, label = as.character(740:763)), 
                  col = "grey25", force = 4, nudge_x = .005, nudge_y = .005)

# Merge in NMDS scores into env_dat
env_dat <- left_join(env_dat, soilmds$points %>% data.frame %>% tibble::rownames_to_column("site"))
env_dat$ca_mg <- env_dat$Ca_ppm/env_dat$Mg_ppm
# Now we merge this env_vars_for_analysis df with the performance data frame `dat`
dat <- left_join(dat, env_dat)
dat_sum <- left_join(dat_sum, env_dat)
# dat_plot_yesno <- dat %>% group_by(site, plot_type,species) %>% mutate(seed_production = if_else(seed_production))

# PLOTS -----------
# Explore patterns in germination and lambda
library(gridExtra)
gg_plot_l <- ggplot(dat_sum, aes(x = plot_num, y = log_lambda)) + 
  geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, scales = "free_y", ncol = 4) + geom_line() + ggtitle("Lambda  across plots")

gg_plot_g <- ggplot(dat_sum, aes(x = plot_num, y = mean_germ_percentage)) + 
  geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, scales = "free_y", ncol = 4) + geom_line() + 
  ggtitle("Germination rates across plots ")

# grid.arrange(gg_plot_l, gg_plot_g, ncol = 2)

gg_camg_l <- ggplot(dat_sum, aes(x = ca_mg, y = log_lambda)) + geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, scales = "free_y", ncol = 4) + geom_smooth(method = "lm") + scale_x_log10()
gg_camg_g <- ggplot(dat_sum, aes(x = ca_mg, y = mean_germ_percentage)) + geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, ncol = 4) + geom_smooth(method = "lm") + scale_x_log10()
grid.arrange(gg_camg_l, gg_camg_g, ncol = 2)

gg_depth_l <- ggplot(dat_sum, aes(x = depth, y = log_lambda)) + geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, scales = "free_y", ncol = 4) + geom_smooth(method = "lm")
gg_depth_g <- ggplot(dat_sum, aes(x = depth, y = mean_germ_percentage)) + geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species,  ncol = 4) + geom_smooth(method = "lm")

gg_sand_l <- ggplot(dat_sum, aes(x = sand, y = log_lambda)) + geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, scales = "free_y", ncol = 4) + geom_smooth(method = "lm")
gg_sand_g <- ggplot(dat_sum, aes(x = sand, y = mean_germ_percentage)) + geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species,  ncol = 4) + geom_smooth(method = "lm")

gg_nitrate_l <- ggplot(dat_sum, aes(x = Nitrate_ppm, y = log_lambda)) + geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, scales = "free_y", ncol = 4) + geom_smooth(method = "lm",formula = y ~ x + I(x^2)) 
gg_nitrate_g <- ggplot(dat_sum, aes(x = Nitrate_ppm, y = mean_germ_percentage)) + geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, ncol = 4) + geom_smooth(method = "lm") 

gg_nitrate_l <- ggplot(dat_sum, aes(x = K_ppm, y = log_lambda)) + geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, scales = "free_y", ncol = 4) + geom_smooth(method = "lm",formula = y ~ x + I(x^2)) 
gg_nitrate_g <- ggplot(dat_sum, aes(x = K_ppm, y = mean_germ_percentage)) + geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, ncol = 4) + geom_smooth(method = "lm") 


gg_mds1_l <- ggplot(dat_sum, aes(x = MDS1, y = log_lambda)) + 
  geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, scales = "free_y", ncol = 4) + geom_smooth(method = "lm")
gg_mds1_g <- ggplot(dat_sum, aes(x = MDS1, y = mean_germ_percentage)) + 
  geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, ncol = 4) + geom_smooth(method = "lm")

gg_mds2_l <- ggplot(dat_sum, aes(x = MDS2, y = log_lambda)) + geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, scales = "free_y", ncol = 4) + geom_smooth(method = "lm")
gg_mds2_g <- ggplot(dat_sum, aes(x = MDS2, y = mean_germ_percentage)) + geom_point(alpha = 1/2, size = 3) + 
  facet_wrap(~species, ncol = 4) + geom_smooth(method = "lm")


# Now we merge in the species traits ----------
if(thinkpad){
  traits <- read_delim("data/trait/merged_sp_avg.csv", delim = ",")
} else {
  traits <- read_delim("~/Dropbox/spatial_tapioca/data/trait/merged_sp_avg.csv",
                        delim = ",") 
}
traits <- traits %>% mutate(species = ifelse(species == "SACA", "SACO", species))
traits <- traits %>% mutate(species = ifelse(species == "VUMA", "VUMI", species))
traits <- traits %>% select(-year_measured)
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}
#traits <- traits %>% filter(species %in% unique(merged_df$species)) %>% arrange(species)
traits %>% tbl_df
traits_logged <- traits %>% mutate_if(is.numeric, funs(log(abs(.))))
colnames(traits_logged)[2:18] <- paste0("log_",colnames(traits_logged)[2:18])
traits_logged_scaled <- traits %>% data.frame %>% mutate_if(is.numeric, scale_this)
colnames(traits_logged_scaled)[2:18] <- paste0("scaled_",colnames(traits_logged)[2:18])
dat_sum <- left_join(dat_sum, traits)
dat_sum <- left_join(dat_sum, traits_logged)
dat_sum <- left_join(dat_sum, traits_logged_scaled)

scaled_ls <- dat_sum %>% select(site, species, lambda) %>% tidyr::spread(species, lambda) %>% 
  data.frame %>% mutate_if(is.numeric, scale_this) %>%  #unite(id, c("site", "replicate"), sep = "XX") %>%
  tidyr::gather("species", "scaled_lambda",2:18) 

dat_sum <- left_join(dat_sum, scaled_ls) 
dat_sum_sub <- dat_sum %>% filter(species %in% c("PLER", "LACA", "SACO", "AGHE", 
                                                 "HOMU", "CEME", "HECO", "EUPE", 
                                                 "LOWR", "VUMA", "CHGL")) %>%
  filter(site %in% paste0("plot_", 740:755))

colnames(dat_sum)

sla_fit <- lmer(log_lambda ~ scaled_log_sla_cm2_g*MDS1 + scaled_log_sla_cm2_g*MDS2 + 
              (1+species) + (1|site), data = dat_sum)
sla_fit_sub <- lmer(log_lambda ~ scaled_log_sla_cm2_g*MDS1 + scaled_log_sla_cm2_g*MDS2 + 
                   (1+species) + (1|site), data = dat_sum_sub)
summary(sla_fit)
summary(sla_fit_sub)
car::Anova(sla_fit)
car::Anova(sla_fit_sub)

seed_fit <- lmer(log_lambda ~ scaled_log_seed_mass_g*MDS1 + scaled_log_seed_mass_g*MDS2 + 
                   (1+species) + (1|site), data = dat_sum)
seed_fit_sub <- lmer(log_lambda ~ scaled_log_seed_mass_g*MDS1 + scaled_log_seed_mass_g*MDS2 + 
                    (1+species) + (1|site), data = dat_sum_sub)
summary(seed_fit)
summary(seed_fit_sub)
car::Anova(seed_fit)
car::Anova(seed_fit_sub)

seed_fit_g <- lmer(mean_germ_percentage ~ scaled_log_seed_mass_g*MDS1 + scaled_log_seed_mass_g*MDS2 + 
                    (1+species) + (1|site), data = dat_sum)
seed_fit_sub_g <- lmer(mean_germ_percentage ~ scaled_log_seed_mass_g*MDS1 + scaled_log_seed_mass_g*MDS2 + 
                      (1+species) + (1|site), data = dat_sum_sub)
summary(seed_fit_g)
summary(seed_fit_sub_g)
car::Anova(seed_fit_g)
car::Anova(seed_fit_sub_g)


sla_fit_camg <- lmer(log_lambda ~ scaled_log_sla_cm2_g*scale(ca_mg) + 
                   (1+species) + (1|site), data = dat_sum)
sla_fit_camg_sub <- lmer(log_lambda ~ scaled_log_sla_cm2_g*scale(ca_mg) +
                       (1+species) + (1|site), data = dat_sum_sub)

summary(sla_fit_camg)
summary(sla_fit_camg_sub)
car::Anova(sla_fit_camg)
car::Anova(sla_fit_camg_sub)


seed_fit_camg <- lmer(log_lambda ~ scaled_log_seed_mass_g*scale(ca_mg) + 
                       (1+species) + (1|site), data = dat_sum)
seed_fit_camg_sub <- lmer(log_lambda ~ scaled_log_seed_mass_g*scale(ca_mg) +
                           (1+species) + (1|site), data = dat_sum_sub)
summary(seed_fit_camg)
summary(seed_fit_camg_sub)
car::Anova(seed_fit_camg)
car::Anova(seed_fit_camg_sub)


seed_fit_camg_g <- lmer(mean_germ_percentage ~ scaled_log_seed_mass_g*scale(ca_mg) + 
                     (1+species) + (1|site), data = dat_sum)
seed_fit_sub_camg_g <- lmer(mean_germ_percentage ~ scaled_log_seed_mass_g*scale(ca_mg) +
                         (1+species) + (1|site), data = dat_sum_sub)

summary(seed_fit_camg_g)
summary(seed_fit_sub_camg_g)
car::Anova(seed_fit_camg_g)
car::Anova(seed_fit_sub_camg_g)
