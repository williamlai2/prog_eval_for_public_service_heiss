# edited version of https://www.andrewheiss.com/blog/2021/12/01/multilevel-models-panel-data-guide/

library(tidyverse)    # ggplot, dplyr, %>%, and friends
library(gapminder)    # Country-year panel data from the Gapminder Project
library(brms)         # Bayesian modeling through Stan
library(tidybayes)    # Manipulate Stan objects in a tidy way
library(broom)        # Convert model objects to data frames
library(broom.mixed)  # Convert brms model objects to data frames
library(emmeans)      # Calculate marginal effects in even fancier ways
library(ggh4x)        # For nested facets in ggplot
library(ggrepel)      # For nice non-overlapping labels in ggplot
library(ggdist)       # For distribution-related ggplot geoms
library(scales)       # For formatting numbers with comma(), dollar(), etc.
library(patchwork)    # For combining plots
library(ggokabeito)   # Colorblind-friendly color palette

# Make all the random draws reproducible
set.seed(1234)

# Bayes stuff
# Use the cmdstanr backend for Stan because it's faster and more modern than
# the default rstan. You need to install the cmdstanr package first
# (https://mc-stan.org/cmdstanr/) and then run cmdstanr::install_cmdstan() to
# install cmdstan on your computer.
# options(mc.cores = 4,  # Use 4 cores
#         brms.backend = "cmdstanr")

# run it this way if it crashes using cmdstanr
options(mc.cores = parallel::detectCores() - 1)

bayes_seed <- 1234


# Custom ggplot theme to make pretty plots
# Get Barlow Semi Condensed at https://fonts.google.com/specimen/Barlow+Semi+Condensed
theme_clean <- function() {
  theme_minimal(base_family = "Barlow Semi Condensed") +
    theme(panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = rel(0.8), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.title = element_text(face = "bold"))
}

# Make labels use Barlow by default
update_geom_defaults("label_repel", 
                     list(family = "Barlow Semi Condensed",
                          fontface = "bold"))
update_geom_defaults("label", 
                     list(family = "Barlow Semi Condensed",
                          fontface = "bold"))

# The ggh4x paackage includes a `facet_nested()` function for nesting facets
# (like countries in continents). Throughout this post, I want the
# continent-level facets to use bolder text and a lighter gray strip. I don't
# want to keep repeating all these settings, though, so I create a list of the
# settings here with `strip_nested()` and feed it to `facet_nested_wrap()` later
nested_settings <- strip_nested(
  text_x = list(element_text(family = "Barlow Semi Condensed Black", 
                             face = "plain"), NULL),
  background_x = list(element_rect(fill = "grey92"), NULL),
  by_layer_x = TRUE)


# gapminder data ---------------------------------

# Little dataset of 8 countries (2 for each of the 4 continents in the data)
# that are good examples of different trends and intercepts
countries <- tribble(
  ~country,       ~continent,
  "Egypt",        "Africa",
  "Sierra Leone", "Africa",
  "Pakistan",     "Asia",
  "Yemen, Rep.",  "Asia",
  "Bolivia",      "Americas",
  "Canada",       "Americas",
  "Italy",        "Europe",
  "Portugal",     "Europe"
)

# Clean up the gapminder data a little
gapminder <- gapminder::gapminder %>%
  # Remove Oceania since there are only two countries there and we want bigger
  # continent clusters
  filter(continent != "Oceania") %>%
  # Scale down GDP per capita so it's more interpretable ("a $1,000 increase in
  # GDP" vs. "a $1 increase in GDP")
  # Also log it
  mutate(gdpPercap_1000 = gdpPercap / 1000,
         gdpPercap_log = log(gdpPercap)) %>% 
  mutate(across(starts_with("gdp"), list("z" = ~scale(.)))) %>% 
  # Make year centered on 1952 (so we're counting the years since 1952). This
  # (1) helps with interpretability, since the intercept will show the average
  # at 1952 instead of the average at 0 CE, and (2) helps with estimation speed
  # since brms/Stan likes to work with small numbers
  mutate(year_orig = year,
         year = year - 1952) %>% 
  # Indicator for the 8 countries we're focusing on
  mutate(highlight = country %in% countries$country)

# Extract rows for the example countries
original_points <- gapminder %>% 
  filter(country %in% countries$country) %>% 
  # Use real years
  mutate(year = year_orig)

# The effect of continent, country, and time on life expectancy ------
ggplot(gapminder, aes(x = year_orig, y = lifeExp, 
                      group = country, color = continent)) +
  geom_line(aes(size = highlight)) +
  geom_smooth(method = "lm", aes(color = NULL, group = NULL), 
              color = "grey60", size = 1, linetype = "21",
              se = FALSE, show.legend = FALSE) +
  geom_label_repel(data = filter(gapminder, year == 0, highlight == TRUE), 
                   aes(label = country), direction = "y", size = 3, seed = 1234, 
                   show.legend = FALSE) +
  annotate(geom = "label", label = "Global trend", x = 1952, y = 50,
           size = 3, color = "grey60") +
  scale_size_manual(values = c(0.075, 1), guide = "none") +
  scale_color_okabe_ito(order = c(2, 3, 6, 1)) +
  labs(x = NULL, y = "Life expectancy", color = "Continent") +
  theme_clean() +
  theme(legend.position = "bottom")

# regular regression -----------
model_boring <- brm(
  bf(lifeExp ~ year),
  data = gapminder,
  chains = 4, seed = bayes_seed
)

tidy(model_boring)

pred_model_boring <- model_boring %>%
  epred_draws(newdata = expand_grid(country = countries$country,
                                    year = unique(gapminder$year))) %>% 
  mutate(year = year + 1952) %>% 
  left_join(countries, by = "country")

ggplot(pred_model_boring, aes(x = year, y = .epred)) +
  geom_point(data = original_points, aes(y = lifeExp), 
             color = "grey50", size = 3, alpha = 0.5) +
  stat_lineribbon(alpha = 0.5, size = 0.5) +
  scale_fill_brewer(palette = "Reds") +
  labs(title = "Global year trend with no country-based variation",
       subtitle = "lifeExp ~ year",
       x = NULL, y = "Predicted life expectancy") +
  guides(fill = "none") +
  facet_nested_wrap(vars(continent, country), nrow = 2, strip = nested_settings) +
  theme_clean() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(family = "Consolas"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

# Introduction to random effects: Intercepts for each continent -------
model_super_boring <- lm(lifeExp ~ year + continent, data = gapminder)
tidy(model_super_boring)

# Fixed (population-level) effects
# Each continent gets its own intercept - a multi level model
model_continent_only <- brm(
  bf(lifeExp ~ year + (1 | continent)),
  data = gapminder,
  control = list(adapt_delta = 0.95),
  chains = 4, seed = bayes_seed
)

tidy(model_continent_only)

summary(model_continent_only)

tidy(model_continent_only, effects = "fixed") 
tidy(model_continent_only, effects = "ran_pars")

# Continent-level random effects
continent_offsets <- ranef(model_continent_only)$continent %>% 
  as_tibble(rownames = "continent")
continent_offsets

coef(model_continent_only)$continent %>% 
  as_tibble(rownames = "continent") %>% 
  select(continent, starts_with("Estimate"))

# Visualise continent-specific trends
newdata_country_continent <- expand_grid(country = countries$country,
                                         year = unique(gapminder$year)) %>% 
  left_join(countries, by = "country")

pred_model_continent_only <- model_continent_only %>%
  epred_draws(newdata_country_continent, re_formula = NULL) %>% 
  mutate(year = year + 1952)

ggplot(pred_model_continent_only, aes(x = year, y = .epred)) +
  geom_point(data = original_points, aes(y = lifeExp),
             color = "grey50", size = 3, alpha = 0.5) +
  stat_lineribbon(alpha = 0.5) +
  scale_fill_brewer(palette = "Reds") +
  labs(title = "Intercepts for year trend vary by continent",
       subtitle = "lifeExp ~ year + (1 | continent)",
       x = NULL, y = "Predicted life expectancy") +
  guides(fill = "none") +
  facet_nested_wrap(vars(continent, country), nrow = 2, strip = nested_settings) +
  theme_clean() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(family = "Consolas"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

# Intercepts for each country -------------
# Each country gets its own intercept
model_country_only <- brm(
  bf(lifeExp ~ year + (1 | country)),
  data = gapminder,
  chains = 4, seed = bayes_seed,
  iter = 4000  # Double the number of iterations to help with convergence
)

tidy(model_country_only)

country_offsets <- ranef(model_country_only)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  select(country, starts_with("Estimate"))
country_offsets

coef(model_country_only)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  select(country, starts_with("Estimate"))

pred_model_country_only <- model_country_only %>%
  epred_draws(newdata = expand_grid(country = countries$country,
                                    year = unique(gapminder$year)),
              re_formula = NULL) %>% 
  mutate(year = year + 1952) %>% 
  left_join(countries, by = "country")

ggplot(pred_model_country_only, aes(x = year, y = .epred)) +
  geom_point(data = original_points, aes(y = lifeExp), 
             color = "grey50", size = 3, alpha = 0.5) +
  stat_lineribbon(alpha = 0.5) +
  scale_fill_brewer(palette = "Reds") +
  labs(title = "Intercepts for year trend vary by country",
       subtitle = "lifeExp ~ year + (1 | country)",
       x = NULL, y = "Predicted life expectancy") +
  guides(fill = "none") +
  facet_nested_wrap(vars(continent, country), nrow = 2, strip = nested_settings) +
  theme_clean() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(family = "Consolas"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

# Intercepts and slopes for each country -----
# Each country gets its own slope and intercept for the year trend
model_country_year <- brm(
  bf(lifeExp ~ year + (1 + year | country)),
  data = gapminder,
  chains = 4, seed = bayes_seed,
  iter = 4000  # Double the number of iterations to help with convergence
)

tidy(model_country_year) 

ggplot(data = tibble(x = seq(-0.3, 1, by = 0.1)), aes(x = x)) +
  stat_function(geom = "area",
                fun = dnorm, args = list(mean = 0.327, sd = 0.161),
                fill = palette_okabe_ito(order = 6)) +
  geom_vline(xintercept = 0) +
  labs(title = "Year trend across countries",
       x = "Annual increase in life expectancy") +
  theme_clean()

country_year_offsets <- ranef(model_country_year)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  select(country, starts_with("Estimate"))
country_year_offsets

coef(model_country_year)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  select(country, starts_with("Estimate"))


pred_model_country_year <- model_country_year %>%
  epred_draws(newdata = expand_grid(country = countries$country,
                                    year = unique(gapminder$year)),
              re_formula = NULL) %>% 
  mutate(year = year + 1952) %>% 
  left_join(countries, by = "country")

ggplot(pred_model_country_year, aes(x = year, y = .epred)) +
  geom_point(data = original_points, aes(y = lifeExp), 
             color = "grey50", size = 3, alpha = 0.5) +
  stat_lineribbon(alpha = 0.5) +
  scale_fill_brewer(palette = "Reds") +
  labs(title = "Intercepts and slopes for year trend vary by country",
       subtitle = "lifeExp ~ year + (1 + year | country)",
       x = NULL, y = "Predicted life expectancy") +
  guides(fill = "none") +
  facet_nested_wrap(vars(continent, country), nrow = 2, strip = nested_settings) +
  theme_clean() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(family = "Consolas"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

# Intercepts and slopes for each country + account for year-specific differences -------------
# Each country gets its own slope and intercept for the year trend
model_country_year_year <- brm(
  bf(lifeExp ~ year + (1 | year) + (1 + year | country)),
  data = gapminder,
  chains = 4, seed = bayes_seed,
  iter = 4000  # Double the number of iterations to help with convergence
)

tidy(model_country_year_year)

year_offsets <- ranef(model_country_year_year)$year %>%
  as_tibble(rownames = "year") %>% 
  mutate(year = as.numeric(year) + 1952) %>% 
  select(year, starts_with("Estimate"))
year_offsets

model_country_year_year %>% 
  emtrends(~ year + country,
           var = "year",
           at = list(year = c(0, 5), country = countries$country),
           epred = TRUE, re_formula = NULL, allow_new_levels = TRUE)

pred_model_country_year_year <- model_country_year_year %>%
  epred_draws(newdata = expand_grid(country = countries$country,
                                    year = unique(gapminder$year)),
              re_formula = NULL) %>% 
  mutate(year = year + 1952) %>% 
  left_join(countries, by = "country")

ggplot(pred_model_country_year_year, aes(x = year, y = .epred)) +
  geom_point(data = original_points, aes(y = lifeExp), 
             color = "grey50", size = 3, alpha = 0.5) +
  stat_lineribbon(alpha = 0.5) +
  scale_fill_brewer(palette = "Reds") +
  labs(title = "Intercepts and slopes for year trend vary by country and intercepts vary by year",
       subtitle = "lifeExp ~ year + (1 | year) + (1 + year | country)",
       x = NULL, y = "Predicted life expectancy") +
  guides(fill = "none") +
  facet_nested_wrap(vars(continent, country), nrow = 2, strip = nested_settings) +
  theme_clean() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(family = "Consolas"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


# Intercepts and slopes for each country and continent -------

# # This takes a while! It also has a bunch of divergent transitions, and pretty
# # much all the chains hit the maximum treedepth limit, but this is just a toy
# # example, so whatever
# model_country_continent_year <- brm(
#   bf(lifeExp ~ year + (1 + year | continent / country)),
#   data = gapminder,
#   chains = 4, seed = bayes_seed,
#   iter = 4000  # Double the number of iterations to help with convergence
# )

# You'd use this to calculate the continent/country effects
model_country_continent_year %>%
  emmeans(~ year + continent:country,
          at = list(year = c(0), country = countries$country),
          nesting = "country %in% continent",
          epred = TRUE, re_formula = NULL, allow_new_levels = TRUE)


# Quick digression on logging, scaling, and centreing ------------
# more of an explanation on the blog in this bit
# https://www.andrewheiss.com/blog/2021/12/01/multilevel-models-panel-data-guide/#2-scale-and-center
different_gdps <- gapminder %>% 
  select(gdpPercap, gdpPercap_1000, gdpPercap_z, gdpPercap_log, gdpPercap_log_z) %>% 
  pivot_longer(everything()) %>% 
  mutate(name_nice = recode(
    name, "gdpPercap" = "GDP per capita",
    "gdpPercap_1000" = "GDP per capita ($1,000)",
    "gdpPercap_z" = "GDP per capita (centered & scaled by one standard deviation)",
    "gdpPercap_log" = "GDP per capita (logged)",
    "gdpPercap_log_z" = "GDP per capita (logged; centered & scaled by one standard deviation)")) %>% 
  mutate(name_nice = fct_inorder(name_nice)) %>% 
  mutate(type = ifelse(str_detect(name, "log"), "Logged", "Original scale")) %>% 
  mutate(vline = ifelse(name == "gdpPercap_log", NA, 0))

ggplot(different_gdps, aes(x = value, fill = type)) +
  geom_density(color = NA) +
  geom_vline(data = different_gdps %>% drop_na(vline) %>% 
               group_by(name_nice) %>% slice(1),
             aes(xintercept = vline)) +
  scale_fill_viridis_d(option = "rocket", begin = 0.3, end = 0.6) +
  guides(fill = "none") +
  labs(x = NULL, y = NULL) +
  facet_wrap(vars(name_nice), nrow = 3, scales = "free", dir = "v") +
  theme_clean()

# The effect of wealth on health, accounting for country and time ------------
ggplot(gapminder, aes(x = gdpPercap, y = lifeExp, color = continent)) + 
  geom_point(size = 0.5, alpha = 0.25) +
  geom_smooth(method = "lm", aes(color = NULL), color = "grey60", size = 0.5, 
              se = FALSE, show.legend = FALSE) +
  annotate(geom = "label", label = "Global trend", x = 64000, y = 84, 
           size = 3, color = "grey60") +
  geom_path(aes(group = country, size = highlight),
            arrow = arrow(type = "open", angle = 30, length = unit(0.75, "lines")),
            show.legend = FALSE) +
  geom_label_repel(data = filter(gapminder, year == 15, highlight == TRUE), 
                   aes(label = country), size = 3, seed = 1234, 
                   show.legend = FALSE,
                   family = "Barlow Semi Condensed", fontface = "bold") +
  scale_size_manual(values = c(0.075, 1), guide = "none") +
  scale_color_okabe_ito(order = c(2, 3, 6, 1),
                        guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_x_log10(labels = dollar_format(accuracy = 1), breaks = 125 * (2^(1:10))) +
  labs(x = "GDP per capita (log)", y = "Life expectancy", color = "Continent") +
  theme_clean() +
  theme(legend.position = "bottom")

# used to unscale later
gdp_mean_sd <- attributes(gapminder$gdpPercap_z)
gdp_mean <- gdp_mean_sd$`scaled:center`
gdp_sd <- gdp_mean_sd$`scaled:scale`

gdp_log_mean_sd <- attributes(gapminder$gdpPercap_log_z)
gdp_log_mean <- gdp_log_mean_sd$`scaled:center`
gdp_log_sd <- gdp_log_mean_sd$`scaled:scale`

# Regular regression ----------
model_gdp_boring <- brm(
  bf(lifeExp ~ gdpPercap_z + year),
  data = gapminder,
  chains = 4, seed = bayes_seed
)


model_gdp_boring_log <- brm(
  bf(lifeExp ~ gdpPercap_log_z + year),
  data = gapminder,
  chains = 4, seed = bayes_seed
)

tidy(model_gdp_boring)
tidy(model_gdp_boring_log)

# back to original scales
tidy(model_gdp_boring) %>% 
  mutate(across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term == "gdpPercap_z", (. / gdp_sd) * 1000, .)))

tidy(model_gdp_boring_log) %>% 
  mutate(across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term == "gdpPercap_log_z", (. / gdp_log_sd), .)))

ame_model_gdp_boring <- model_gdp_boring %>% 
  emtrends(~ 1,
           var = "gdpPercap_z",
           at = list(year = 0),
           epred = TRUE, re_formula = NULL)

pred_ame_model_gdp_boring <- ame_model_gdp_boring %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = .value / gdp_sd * 1000) %>% 
  mutate(fake_facet_title = "GDP per capita")

plot_ame_gdp <- ggplot(pred_ame_model_gdp_boring, aes(x = .value)) +
  stat_halfeye(fill = palette_okabe_ito(5),
               .width = c(0.8, 0.95)) +
  labs(x = paste0("Average marginal effect on life expectancy", "\n", 
                  "of a $1,000 increase in GDP per capita"), 
       y = "Density") +
  facet_wrap(vars(fake_facet_title)) +
  theme_clean() +
  theme(strip.text = element_text(size = rel(1.1)))


ame_model_gdp_boring_log <- model_gdp_boring_log %>% 
  emtrends(~ 1,
           var = "gdpPercap_log_z",
           at = list(year = 0),
           epred = TRUE, re_formula = NULL)

pred_ame_model_gdp_boring_log <- ame_model_gdp_boring_log %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = .value / gdp_log_sd / 10) %>% 
  mutate(fake_facet_title = "Logged GDP per capita")

plot_ame_gdp_log <- ggplot(pred_ame_model_gdp_boring_log, aes(x = .value)) +
  stat_halfeye(fill = palette_okabe_ito(5),
               .width = c(0.8, 0.95)) +
  labs(x = paste0("Average marginal effect on life expectancy", "\n", 
                  "of a 10% increase in GDP per capita"), 
       y = "Density") +
  facet_wrap(vars(fake_facet_title)) +
  theme_clean() +
  theme(strip.text = element_text(size = rel(1.1)))

(plot_ame_gdp + plot_ame_gdp_log) +
  plot_annotation(title = paste0("Average marginal effect of GDP per capita", "\n",
                                 "(includes year trend with no country-based variation)"), 
                  subtitle = "lifeExp ~ gdpPercap_log_z + year",
                  caption = "80% and 95% credible intervals shown with error bar",
                  theme = theme_clean() + theme(plot.subtitle = element_text(family = "Consolas")))

# Each country gets its own intercept and GDP slope ----
model_gdp_country_only <- brm(
  bf(lifeExp ~ gdpPercap_z + year + (1 + gdpPercap_z | country),
     decomp = "QR"),
  data = gapminder,
  chains = 4, seed = bayes_seed,
  threads = threading(2)  # Two CPUs per chain to speed things up
)

model_gdp_country_only_log <- brm(
  bf(lifeExp ~ gdpPercap_log_z + year + (1 + gdpPercap_log_z | country),
     decomp = "QR"),
  data = gapminder,
  chains = 4, seed = bayes_seed,
  threads = threading(2)  # Two CPUs per chain to speed things up
)

# Model with regular GDP per capita
rstan::get_elapsed_time(model_gdp_country_only$fit) %>% 
  as_tibble(rownames = "chain") %>% mutate(total_seconds = warmup + sample)

# Model with logged GDP per capita 
rstan::get_elapsed_time(model_gdp_country_only_log$fit) %>% 
  as_tibble(rownames = "chain") %>% mutate(total_seconds = warmup + sample) # faster

# Unscale both the GDP coefficient and the GDP random variance coefficient
tidy(model_gdp_country_only) %>% 
  mutate(across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term %in% c("gdpPercap_z", "sd__gdpPercap_z"), 
                        (. / gdp_sd) * 1000, .)))

tidy(model_gdp_country_only_log) %>% 
  mutate(across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term %in% c("gdpPercap_log_z", "sd__gdpPercap_log_z"), 
                        (. / gdp_log_sd), .)))

ranef(model_gdp_country_only)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  select(country, starts_with("Estimate")) %>% 
  # Unscale the GDP offsets
  mutate(Estimate.gdpPercap_z = Estimate.gdpPercap_z / gdp_sd * 1000)

coef(model_gdp_country_only)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  select(country, starts_with("Estimate")) %>% 
  # Unscale the GDP offsets
  mutate(Estimate.gdpPercap_z = Estimate.gdpPercap_z / gdp_sd * 1000)

ame_model_gdp_country_only <- model_gdp_country_only %>% 
  emtrends(~ country,
           var = "gdpPercap_z",
           at = list(year = 0, country = countries$country),
           epred = TRUE, re_formula = NULL)

pred_ame_model_gdp_country_only <- ame_model_gdp_country_only %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = .value / gdp_sd * 1000) %>% 
  left_join(countries, by = "country")

ggplot(pred_ame_model_gdp_country_only, aes(x = .value)) +
  stat_halfeye(aes(fill = continent)) +
  geom_vline(xintercept = 0) +
  scale_fill_okabe_ito(order = c(2, 3, 6, 1), guide = "none") +
  labs(title = paste0("Average marginal effect of GDP per capita", "\n",
                      "(intercepts and slope of GDP per capita by country)"), 
       subtitle = "lifeExp ~ gdpPercap_z + year + (1 + gdpPercap_z | country)",
       x = paste0("Average marginal effect on life expectancy", "\n", 
                  "of a $1,000 increase in GDP per capita"), 
       y = "Density",
       caption = "80% and 95% credible intervals shown with error bar") +
  facet_nested_wrap(vars(continent, country), nrow = 2, strip = nested_settings) +
  theme_clean() +
  theme(plot.subtitle = element_text(family = "Consolas"))

ame_model_gdp_country_only_log <- model_gdp_country_only_log %>% 
  emtrends(~ country,
           var = "gdpPercap_log_z",
           at = list(year = 0, country = countries$country),
           epred = TRUE, re_formula = NULL)
ame_model_gdp_country_only_log

pred_ame_model_gdp_country_only_log <- ame_model_gdp_country_only_log %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = .value / gdp_log_sd / 10) %>% 
  left_join(countries, by = "country")

ggplot(pred_ame_model_gdp_country_only_log, aes(x = .value)) +
  stat_halfeye(aes(fill = continent)) +
  geom_vline(xintercept = 0) +
  scale_fill_okabe_ito(order = c(2, 3, 6, 1), guide = "none") +
  labs(title = paste0("Average marginal effect of GDP per capita", "\n",
                      "(intercepts and slope of GDP per capita by country)"), 
       subtitle = "lifeExp ~ gdpPercap_z + year + (1 + gdpPercap_z | country)",
       x = paste0("Average marginal effect on life expectancy", "\n", 
                  "of a 10% increase in GDP per capita"), 
       y = "Density",
       caption = "80% and 95% credible intervals shown with error bar") +
  facet_nested_wrap(vars(continent, country), nrow = 2, strip = nested_settings) +
  theme_clean() +
  theme(plot.subtitle = element_text(family = "Consolas"))

# Each country gets its own intercept and GDP and year slopes -------
model_gdp_country_year <- brm(
  bf(lifeExp ~ gdpPercap_z + year + (1 + gdpPercap_z + year | country),
     decomp = "QR"),
  data = gapminder,
  chains = 4, seed = bayes_seed,
  threads = threading(2)  # Two CPUs per chain to speed things up
)


model_gdp_country_year_log <- brm(
  bf(lifeExp ~ gdpPercap_log_z + year + (1 + gdpPercap_log_z + year | country),
     decomp = "QR"),
  data = gapminder,
  chains = 4, seed = bayes_seed,
  threads = threading(2)  # Two CPUs per chain to speed things up
)

# Model with regular GDP per capita
rstan::get_elapsed_time(model_gdp_country_year$fit) %>% 
  as_tibble(rownames = "chain") %>% mutate(total_seconds = warmup + sample)

# Model with logged GDP per capita 
rstan::get_elapsed_time(model_gdp_country_year_log$fit) %>% 
  as_tibble(rownames = "chain") %>% mutate(total_seconds = warmup + sample) # faster

# Bad R-hats and low effective sample sizes
bayestestR::diagnostic_posterior(model_gdp_country_year)

# Okay R-hats and okay-ish effective sample sizes
bayestestR::diagnostic_posterior(model_gdp_country_year_log)

# Unscale both the GDP coefficient and the GDP random variance coefficient
tidy(model_gdp_country_year) %>% 
  mutate(across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term %in% c("gdpPercap_z", "sd__gdpPercap_z"), 
                        (. / gdp_sd) * 1000, .)))

tidy(model_gdp_country_year_log) %>% 
  mutate(across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term %in% c("gdpPercap_log_z", "sd__gdpPercap_log_z"), 
                        (. / gdp_log_sd), .)))

coef(model_gdp_country_year)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  select(country, starts_with("Estimate")) %>% 
  # Unscale the GDP offsets
  mutate(Estimate.gdpPercap_z = Estimate.gdpPercap_z / gdp_sd * 1000)

ame_model_gdp_country_year <- model_gdp_country_year %>% 
  emtrends(~ year + country,
           var = "gdpPercap_z",
           at = list(year = 0, country = countries$country),
           epred = TRUE, re_formula = NULL)

pred_ame_model_gdp_country_year <- ame_model_gdp_country_year %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = .value / gdp_sd * 1000) %>% 
  left_join(countries, by = "country")

ggplot(pred_ame_model_gdp_country_year, aes(x = .value)) +
  stat_halfeye(aes(fill = continent)) +
  geom_vline(xintercept = 0) +
  scale_fill_okabe_ito(order = c(2, 3, 6, 1), guide = "none") +
  labs(title = paste0("Average marginal effect of GDP per capita", "\n",
                      "(intercepts and slopes of both GDP per capita and year vary by country)"), 
       subtitle = "lifeExp ~ gdpPercap_log_z + year + (1 + gdpPercap_z + year | country)",
       x = paste0("Average marginal effect on life expectancy", "\n", 
                  "of a $1,000 increase in GDP per capita"), 
       y = "Density",
       caption = "80% and 95% credible intervals shown with error bar") +
  facet_nested_wrap(vars(continent, country), nrow = 2, strip = nested_settings) +
  theme_clean() +
  theme(plot.subtitle = element_text(family = "Consolas"))

ame_model_gdp_country_year_log <- model_gdp_country_year_log %>% 
  emtrends(~ year + country,
           var = "gdpPercap_log_z",
           at = list(year = 0, country = countries$country),
           epred = TRUE, re_formula = NULL)

pred_ame_model_gdp_country_year_log <- ame_model_gdp_country_year_log %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = .value / gdp_log_sd / 10) %>% 
  left_join(countries, by = "country")

ggplot(pred_ame_model_gdp_country_year_log, aes(x = .value)) +
  stat_halfeye(aes(fill = continent)) +
  geom_vline(xintercept = 0) +
  scale_fill_okabe_ito(order = c(2, 3, 6, 1), guide = "none") +
  labs(title = paste0("Average marginal effect of GDP per capita", "\n",
                      "(intercepts and slopes of both GDP per capita and year vary by country)"), 
       subtitle = "lifeExp ~ gdpPercap_log_z + year + (1 + gdpPercap_log_z + year | country)",
       x = paste0("Average marginal effect on life expectancy", "\n", 
                  "of a 10% increase in GDP per capita"), 
       y = "Density",
       caption = "80% and 95% credible intervals shown with error bar") +
  facet_nested_wrap(vars(continent, country), nrow = 2, strip = nested_settings) +
  theme_clean() +
  theme(plot.subtitle = element_text(family = "Consolas"))

# Calculate the overall global effect
# `re_formula = NA` means that no random effects will be incorporated
ame_global_gdp_country_year <- model_gdp_country_year_log %>% 
  emtrends(~ 1,
           var = "gdpPercap_log_z",
           epred = TRUE, re_formula = NA)
ame_global_gdp_country_year
##  1       gdpPercap_log_z.trend lower.HPD upper.HPD
##  overall                  5.27      4.32      6.29
## 
## Point estimate displayed: median 
## HPD interval probability: 0.95

# Get the posterior distribution of this global effect and make a plot
pred_global_gdp_country_year <- ame_global_gdp_country_year %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = (.value / gdp_log_sd) / 10) %>% 
  mutate(fake_facet_title = "Global grand mean")

plot_global_gdp_country_year <- ggplot(pred_global_gdp_country_year, aes(x = .value)) +
  stat_halfeye(fill = palette_okabe_ito(5)) +
  geom_vline(xintercept = 0) +
  labs(x = paste0("Average marginal effect on life expectancy", "\n", 
                  "of a 10% increase in GDP per capita"), 
       y = "Density") +
  facet_wrap(vars(fake_facet_title)) +
  theme_clean() +
  theme(strip.text = element_text(size = rel(1.1)))

# Calculate the effect for Atlantis
# `re_formula = NULL` means the full random effects structure will be
# incorporated. `allow_new_levels` lets R deal with a new country, and
# `sample_new_levels = "gaussian"` means that the characteristics for this new
# country will be based on random draws from the model
ame_hypo_gdp_country_year <- model_gdp_country_year_log %>% 
  emtrends(~ 1 + country,
           var = "gdpPercap_log_z",
           at = list(country = "Atlantis"),
           epred = TRUE, re_formula = NULL, 
           allow_new_levels = TRUE, sample_new_levels = "gaussian")
ame_hypo_gdp_country_year
##  country  gdpPercap_log_z.trend lower.HPD upper.HPD
##  Atlantis                  5.14     -2.93      13.5
## 
## Point estimate displayed: median 
## HPD interval probability: 0.95

# Get the posterior distribution of this Atlantis effect and make a plot
pred_hypo_gdp_country_year <- ame_hypo_gdp_country_year %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = (.value / gdp_log_sd) / 10) %>% 
  mutate(fake_facet_title = "AME for hypothetical Atlantis")

plot_hypo_gdp_country_year <- ggplot(pred_hypo_gdp_country_year, aes(x = .value)) +
  stat_halfeye(fill = palette_okabe_ito(7)) +
  geom_vline(xintercept = 0) +
  labs(x = paste0("Average marginal effect on life expectancy", "\n", 
                  "of a 10% increase in GDP per capita"), 
       y = "Density") +
  facet_wrap(vars(fake_facet_title)) +
  theme_clean() +
  theme(strip.text = element_text(size = rel(1.1)))

# Show the two AME distributions side-by-side
(plot_global_gdp_country_year | plot_hypo_gdp_country_year) +
  plot_annotation(title = paste0("Average marginal effect of GDP per capita", "\n",
                                 "(intercepts and slopes of both GDP per capita and year vary by country)"), 
                  subtitle = "lifeExp ~ gdpPercap_log_z + year + (1 + gdpPercap_log_z + year | country)",
                  caption = "80% and 95% credible intervals shown with error bar",
                  theme = theme_clean() + theme(plot.subtitle = element_text(family = "Consolas")))

model_gdp_country_year_year <- brm(
  bf(lifeExp ~ gdpPercap_log_z + year + (1 + gdpPercap_log_z | year) + 
       (1 + gdpPercap_log_z + year | country),
     decomp = "QR"),
  data = gapminder,
  cores = 4, chains = 4, seed = bayes_seed,
  threads = threading(2)  # Two CPUs per chain to speed things up
)

rstan::get_elapsed_time(model_gdp_country_year_year$fit) %>% 
  as_tibble(rownames = "chain") %>% mutate(total_seconds = warmup + sample)

# Unscale both the GDP coefficient and the GDP random variance coefficient
tidy(model_gdp_country_year_year) %>% 
  mutate(across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term %in% c("gdpPercap_log_z", "sd__gdpPercap_log_z"), 
                        (. / gdp_log_sd), .)))

coef(model_gdp_country_year)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  select(country, starts_with("Estimate")) %>% 
  # Unscale the GDP offsets
  mutate(Estimate.gdpPercap_z = Estimate.gdpPercap_z / gdp_sd * 1000)

# Different slopes in each year!
model_gdp_country_year_year %>% 
  emtrends(~ year + country,
           var = "gdpPercap_log_z",
           at = list(year = c(0, 5), country = countries$country),
           epred = TRUE, re_formula = NULL)

ame_model_gdp_country_year_year <- model_gdp_country_year_year %>% 
  emtrends(~ year + country,
           var = "gdpPercap_log_z",
           at = list(year = seq(0, 55, by = 5), country = countries$country),
           epred = TRUE, re_formula = NULL)

pred_model_gdp_country_year_year <- ame_model_gdp_country_year_year %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = (.value / gdp_log_sd) / 10) %>% 
  mutate(year = year + 1952,
         year = fct_inorder(factor(year))) %>%
  left_join(countries, by = "country")

ggplot(pred_model_gdp_country_year_year, aes(x = .value, fill = year)) +
  stat_slab(aes(y = fct_rev(year), fill = year,
                fill_ramp = stat(cut_cdf_qi(cdf, .width = c(0.02, 0.8, 0.95, 1)))),
            height = 2, color = "white", slab_size = 0.5) +
  scale_fill_viridis_d(option = "rocket", guide = "none", end = 0.9) +
  scale_fill_ramp_discrete(range = c(1, 0.2), guide = "none") +
  geom_vline(xintercept = 0) +
  labs(title = paste0("Average marginal effect of GDP per capita", "\n",
                      "(intercepts and slopes of both GDP per capita and year ",
                      "vary by country + intercepts vary by year)"), 
       subtitle = paste0("lifeExp ~ gdpPercap_log_z + year + ", 
                         "(1 + gdpPercap_log_z | year) + ", 
                         "(1 + gdpPercap_log_z + year | country)"),
       x = paste0("Average marginal effect on life expectancy", "\n", 
                  "of a 10% increase in GDP per capita"), 
       y = "Year",
       caption = "80% and 95% credible intervals shown with shading") +
  facet_nested_wrap(vars(continent, country), nrow = 2, strip = nested_settings,
                    scales = "free_x") +
  theme_clean() +
  theme(panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(family = "Consolas"))

# Calculate the overall global effect
# `re_formula = NA` means that no random effects will be incorporated
ame_global_gdp_country_year_year <- model_gdp_country_year_year %>% 
  emtrends(~ 1,
           var = "gdpPercap_log_z",
           epred = TRUE, re_formula = NA)

# Get the posterior distribution of this global effect and make a plot
pred_global_gdp_country_year_year <- ame_global_gdp_country_year_year %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = (.value / gdp_log_sd) / 10) %>% 
  mutate(fake_facet_title = "Global grand mean")

plot_global_gdp_country_year_year <- ggplot(pred_global_gdp_country_year_year, aes(x = .value)) +
  stat_slab(aes(fill_ramp = stat(cut_cdf_qi(cdf, .width = c(0.02, 0.8, 0.95, 1)))),
            fill = palette_okabe_ito(5)) +
  scale_fill_viridis_d(option = "plasma", guide = "none", begin = 0.5) +
  scale_fill_ramp_discrete(range = c(1, 0.2), guide = "none") +
  geom_vline(xintercept = 0) +
  labs(x = paste0("Average marginal effect on life expectancy", "\n", 
                  "of a 10% increase in GDP per capita"), 
       y = "Density") +
  facet_wrap(vars(fake_facet_title)) +
  theme_clean() +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(size = rel(1.1)))

# Calculate the effect for Atlantis in 2020
# `re_formula = NULL` means the full random effects structure will be
# incorporated. `allow_new_levels` lets R deal with a new country, and
# `sample_new_levels = "gaussian"` means that the characteristics for this new
# country will be based on random draws from the model
ame_hypo1_gdp_country_year_year <- model_gdp_country_year_year %>% 
  emtrends(~ country + year,
           var = "gdpPercap_log_z",
           at = list(year = (2020 - 1952), country = "Atlantis"),
           epred = TRUE, re_formula = NULL, 
           allow_new_levels = TRUE, sample_new_levels = "gaussian")

# Atlantis median 2020 effect
ame_hypo1_gdp_country_year_year %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = (.value / gdp_log_sd) / 10) %>% 
  median_hdci()
## # A tibble: 1 × 8
##   country   year .value .lower .upper .width .point .interval
##   <fct>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <chr>  <chr>    
## 1 Atlantis    68  0.331 -0.220  0.903   0.95 median hdci

# Get the posterior distribution of this Atlantis effect and make a plot
pred_hypo1_gdp_country_year_year <- ame_hypo1_gdp_country_year_year %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = (.value / gdp_log_sd) / 10) %>% 
  mutate(fake_facet_title = "AME for hypothetical Atlantis in 2020")

plot_hypo1_gdp_country_year_year <- ggplot(pred_hypo1_gdp_country_year_year, aes(x = .value)) +
  stat_slab(aes(fill_ramp = stat(cut_cdf_qi(cdf, .width = c(0.02, 0.8, 0.95, 1)))),
            fill = palette_okabe_ito(7)) +
  scale_fill_ramp_discrete(range = c(1, 0.2), guide = "none") +
  geom_vline(xintercept = 0) +
  labs(x = paste0("Average marginal effect on life expectancy", "\n", 
                  "of a 10% increase in GDP per capita"), 
       y = "Density") +
  facet_wrap(vars(fake_facet_title)) +
  theme_clean() +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(size = rel(1.1)))

# Calculate the effect for Atlantis across all existing years
ame_hypo_gdp_country_year_year <- model_gdp_country_year_year %>% 
  emtrends(~ country + year,
           var = "gdpPercap_log_z",
           at = list(year = seq(0, 55, by = 5), country = "Atlantis"),
           epred = TRUE, re_formula = NULL, 
           allow_new_levels = TRUE, sample_new_levels = "gaussian")

# Get the posterior distribution of this Atlantis effect and make a plot
pred_hypo_gdp_country_year_year <- ame_hypo_gdp_country_year_year %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = (.value / gdp_log_sd) / 10) %>% 
  mutate(year = year + 1952,
         year = fct_inorder(factor(year))) %>%
  mutate(fake_facet_title = "AME for hypothetical Atlantis across time")

plot_hypo_gdp_country_year_year <- ggplot(pred_hypo_gdp_country_year_year, aes(x = .value)) +
  stat_slab(aes(y = fct_rev(year), fill = year,
                fill_ramp = stat(cut_cdf_qi(cdf, .width = c(0.02, 0.8, 0.95, 1)))),
            height = 2, color = "white") +
  scale_fill_viridis_d(option = "rocket", guide = "none", end = 0.9) +
  scale_fill_ramp_discrete(range = c(1, 0.2), guide = "none") +
  geom_vline(xintercept = 0) +
  labs(x = paste0("Average marginal effect on life expectancy", "\n", 
                  "of a 10% increase in GDP per capita"), 
       y = "Year") +
  facet_wrap(vars(fake_facet_title)) +
  theme_clean() +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(size = rel(1.1)))

# Show all the AME distributions in a mega plot
((plot_global_gdp_country_year_year / plot_hypo1_gdp_country_year_year) | 
    plot_hypo_gdp_country_year_year) +
  plot_annotation(title = paste0("Average marginal effect of GDP per capita", "\n",
                                 "(intercepts and slopes of both GDP per capita and year ",
                                 "vary by country + intercepts and GDP slopes vary by year)"), 
                  subtitle = paste0("lifeExp ~ gdpPercap_log_z + year + ", 
                                    "(1 + gdpPercap_log_z | year) + ", 
                                    "(1 + gdpPercap_log_z + year | country)"),
                  caption = "80% and 95% credible intervals shown with shading",
                  theme = theme_clean() + theme(plot.subtitle = element_text(family = "Consolas")))

model_gdp_continent_country_year <- brm(
  bf(lifeExp ~ gdpPercap_log_z + year + 
       (1 + gdpPercap_log_z + year | continent / country),
     decomp = "QR"),
  data = gapminder,
  cores = 4, chains = 4, seed = bayes_seed,
  threads = threading(2)  # Two CPUs per chain to speed things up
)

# You'd use this to calculate the continent/country specific effects
model_gdp_continent_country_year %>%
  emtrends(~ year + continent:country,
           var = "gdpPercap_log_z",
           at = list(year = c(0), country = countries$country),
           nesting = "country %in% continent",
           epred = TRUE, re_formula = NULL, allow_new_levels = TRUE)