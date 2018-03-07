# make_site_matrix
# takes in a dataframe with a column named `plot_num`
# returns a dataframe with 24 columns and as many rows as in the original df
# with a `1` in the appropriate column indicating which plot that row 
# is coming from.
make_site_matrix <- function(df) {
  to_return <- df %>% mutate(site1 = ifelse(plot_num==740, 1, 0), 
                             site2 = ifelse(plot_num==741, 1, 0),
                             site3 = ifelse(plot_num==742, 1, 0),
                             site4 = ifelse(plot_num==743, 1, 0),
                             site5 = ifelse(plot_num==744, 1, 0),
                             site6 = ifelse(plot_num==745, 1, 0),
                             site7 = ifelse(plot_num==746, 1, 0),
                             site8 = ifelse(plot_num==747, 1, 0),
                             site9 = ifelse(plot_num==748, 1, 0),
                             site10 = ifelse(plot_num==749, 1, 0),
                             site11 = ifelse(plot_num==750, 1, 0),
                             site12 = ifelse(plot_num==751, 1, 0),
                             site13 = ifelse(plot_num==752, 1, 0),
                             site14 = ifelse(plot_num==753, 1, 0),
                             site15 = ifelse(plot_num==754, 1, 0),
                             site16 = ifelse(plot_num==755, 1, 0),
                             site17 = ifelse(plot_num==756, 1, 0),
                             site18 = ifelse(plot_num==757, 1, 0),
                             site19 = ifelse(plot_num==758, 1, 0),
                             site20 = ifelse(plot_num==759, 1, 0),
                             site21 = ifelse(plot_num==760, 1, 0),
                             site22 = ifelse(plot_num==761, 1, 0),
                             site23 = ifelse(plot_num==762, 1, 0),
                             site24 = ifelse(plot_num==763, 1, 0))
  return(to_return)
}

# make_compdens_matrix
# similar to above, but instead of `1`s, columns have the number of competitors that
# each row was competiting with at that site.
make_compdens_matrix <- function(df) {
  to_return <- df %>%
    mutate(site1 = ifelse(plot_num==740, N_per_neighborhood, 0), 
           site2 = ifelse(plot_num==741, N_per_neighborhood, 0),
           site3 = ifelse(plot_num==742, N_per_neighborhood, 0),
           site4 = ifelse(plot_num==743, N_per_neighborhood, 0),
           site5 = ifelse(plot_num==744, N_per_neighborhood, 0),
           site6 = ifelse(plot_num==745, N_per_neighborhood, 0),
           site7 = ifelse(plot_num==746, N_per_neighborhood, 0),
           site8 = ifelse(plot_num==747, N_per_neighborhood, 0),
           site9 = ifelse(plot_num==748, N_per_neighborhood, 0),
           site10 = ifelse(plot_num==749, N_per_neighborhood, 0),
           site11 = ifelse(plot_num==750, N_per_neighborhood, 0),
           site12 = ifelse(plot_num==751, N_per_neighborhood, 0),
           site13 = ifelse(plot_num==752, N_per_neighborhood, 0),
           site14 = ifelse(plot_num==753, N_per_neighborhood, 0),
           site15 = ifelse(plot_num==754, N_per_neighborhood, 0),
           site16 = ifelse(plot_num==755, N_per_neighborhood, 0),
           site17 = ifelse(plot_num==756, N_per_neighborhood, 0),
           site18 = ifelse(plot_num==757, N_per_neighborhood, 0),
           site19 = ifelse(plot_num==758, N_per_neighborhood, 0),
           site20 = ifelse(plot_num==759, N_per_neighborhood, 0),
           site21 = ifelse(plot_num==760, N_per_neighborhood, 0),
           site22 = ifelse(plot_num==761, N_per_neighborhood, 0),
           site23 = ifelse(plot_num==762, N_per_neighborhood, 0),
           site24 = ifelse(plot_num==763, N_per_neighborhood, 0))
  return(to_return)
}