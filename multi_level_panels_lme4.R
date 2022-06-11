# edited version of https://www.andrewheiss.com/blog/2021/12/01/multilevel-models-panel-data-guide/
# but with frequentist packages (mostly for speed)

library(tidyverse)    # ggplot, dplyr, %>%, and friends
library(gapminder)    # Country-year panel data from the Gapminder Project
library(broom)        # Convert model objects to data frames
library(broom.mixed)  # Convert lme4 model objects to data frames
library(emmeans)      # Calculate marginal effects in even fancier ways
library(ggh4x)        # For nested facets in ggplot
library(ggrepel)      # For nice non-overlapping labels in ggplot
library(ggdist)       # For distribution-related ggplot geoms
library(scales)       # For formatting numbers with comma(), dollar(), etc.
library(patchwork)    # For combining plots
library(ggokabeito)   # Colourblind-friendly colour palette
library(lme4)         # For the mixed effect models
library(janitor)      # tidy column names

# Make all the random draws reproducible
set.seed(1234)

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
model_boring <- lm(lifeExp ~ year, data = gapminder)

tidy(model_boring)

pred_model_boring <- augment(model_boring,
                             newdata = expand_grid(country = countries$country,
                                                   year = unique(gapminder$year))) %>% 
  mutate(year = year + 1952) %>% 
  left_join(countries, by = "country")
  
ggplot(pred_model_boring, aes(x = year, y = .fitted)) +
  geom_point(data = original_points, aes(y = lifeExp), 
             color = "grey50", size = 3, alpha = 0.5) +
  facet_nested_wrap(vars(continent, country), nrow = 2, strip = nested_settings) + # this moved up
  scale_fill_brewer(palette = "Reds") + 
  labs(title = "Global year trend with no country-based variation",
       subtitle = "lifeExp ~ year",
       x = NULL, y = "Predicted life expectancy") +
  guides(fill = "none") +
  theme_clean() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(family = "Consolas"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

# Introduction to random effects: Intercepts for each continent -------
model_super_boring <- lm(lifeExp ~ year + continent, data = gapminder)
tidy(model_super_boring)

# Fixed (population-level) effects
# Each continent gets its own intercept - a multi level model
model_continent_only <- lmer(lifeExp ~ year + (1 | continent), data = gapminder)

summary(model_continent_only)

fixef(model_continent_only)
ranef(model_continent_only)

coef(model_continent_only)$continent %>% 
  as_tibble(rownames = "continent") %>% 
  clean_names()

# Visualise continent-specific trends
newdata_country_continent <- expand_grid(country = countries$country,
                                         year = unique(gapminder$year)) %>% 
  left_join(countries, by = "country")

pred_model_continent_only <- augment(model_continent_only,
                                     newdata = newdata_country_continent) %>% 
  mutate(year = year + 1952) 

ggplot(pred_model_continent_only, aes(x = year, y = .fitted)) +
  geom_point(data = original_points, aes(y = lifeExp),
             color = "grey50", size = 3, alpha = 0.5) +
  facet_nested_wrap(vars(continent, country), nrow = 2, strip = nested_settings) +
  scale_fill_brewer(palette = "Reds") +
  labs(title = "Intercepts for year trend vary by continent",
       subtitle = "lifeExp ~ year + (1 | continent)",
       x = NULL, y = "Predicted life expectancy") +
  guides(fill = "none") +
  theme_clean() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(family = "Consolas"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


# Intercepts for each country -------------
# Each country gets its own intercept
model_country_only <- lmer(lifeExp ~ year + (1 | country), data = gapminder)

tidy(model_country_only)

country_offsets <- ranef(model_country_only)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  clean_names()
country_offsets

coef(model_country_only)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  clean_names()

pred_model_country_only <- augment(model_country_only,
                                   newdata = expand_grid(country = countries$country,
                                                         year = unique(gapminder$year))) %>% 
  mutate(year = year + 1952) %>% 
  left_join(countries, by = "country")

ggplot(pred_model_country_only, aes(x = year, y = .fitted)) +
  geom_point(data = original_points, aes(y = lifeExp), 
             color = "grey50", size = 3, alpha = 0.5) +
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
model_country_year <- lmer(lifeExp ~ year + (1 + year | country), data = gapminder)

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
  clean_names()
country_year_offsets

coef(model_country_year)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  clean_names()


pred_model_country_year <- augment(model_country_year,
                                   newdata = expand_grid(country = countries$country,
                                                         year = unique(gapminder$year))) %>% 
  mutate(year = year + 1952) %>% 
  left_join(countries, by = "country")

ggplot(pred_model_country_year, aes(x = year, y = .fitted)) +
  geom_point(data = original_points, aes(y = lifeExp), 
             color = "grey50", size = 3, alpha = 0.5) +
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
model_country_year_year <- lmer(lifeExp ~ year + (1 | year) + (1 + year | country), data = gapminder)

tidy(model_country_year_year) 

year_offsets <- ranef(model_country_year_year)$year %>%
  as_tibble(rownames = "year") %>% 
  mutate(year = as.numeric(year) + 1952) %>% 
  clean_names()
year_offsets

# # doesn't seem to work
# model_country_year_year %>%
#   emtrends(~ year + country,
#            var = "year",
#            at = list(year = c(0, 5), country = countries$country),
#            epred = TRUE, re_formula = NULL)

pred_model_country_year_year <- augment(model_country_year_year,
                                        newdata = expand_grid(country = countries$country,
                                                              year = unique(gapminder$year))) %>% 
  mutate(year = year + 1952) %>% 
  left_join(countries, by = "country")

ggplot(pred_model_country_year_year, aes(x = year, y = .fitted)) +
  geom_point(data = original_points, aes(y = lifeExp), 
             color = "grey50", size = 3, alpha = 0.5) +
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

model_country_continent_year <- lmer(lifeExp ~ year + (1 + year | continent / country), data = gapminder)

tidy(model_country_continent_year)
summary(model_country_continent_year)

# # doesn't seem to work
# # You'd use this to calculate the continent/country effects
# model_country_continent_year %>%
#   emmeans(~ year + continent:country,
#           at = list(year = c(0), country = countries$country),
#           nesting = "country %in% continent",
#           allow_new_levels = TRUE,
#           epred = TRUE, re_formula = NULL)

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
model_gdp_boring <- lm(lifeExp ~ gdpPercap_z + year, data = gapminder)

model_gdp_boring_log <- lm(lifeExp ~ gdpPercap_log_z + year, data = gapminder)

tidy(model_gdp_boring)
tidy(model_gdp_boring_log)

# back to original scales
tidy(model_gdp_boring) %>% 
  mutate(conf.low = estimate - 1.96 * std.error,
         conf.high = estimate + 1.96 * std.error,
         across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term == "gdpPercap_z", (. / gdp_sd) * 1000, .)))

tidy(model_gdp_boring_log) %>% 
  mutate(conf.low = estimate - 1.96 * std.error,
         conf.high = estimate + 1.96 * std.error,
         across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term == "gdpPercap_log_z", (. / gdp_log_sd), .)))

ame_model_gdp_boring <- model_gdp_boring %>% 
  emtrends(~ 1,
           var = "gdpPercap_z",
           at = list(year = 0),
           epred = TRUE, re_formula = NULL)

# not bayesian so can't do next bits (omitted)

ame_model_gdp_boring_log <- model_gdp_boring_log %>% 
  emtrends(~ 1,
           var = "gdpPercap_log_z",
           at = list(year = 0),
           epred = TRUE, re_formula = NULL)


# Each country gets its own intercept and GDP slope ----
model_gdp_country_only <- lmer(lifeExp ~ gdpPercap_z + year + (1 + gdpPercap_z | country), data = gapminder)

model_gdp_country_only_log <- lmer(lifeExp ~ gdpPercap_log_z + year + (1 + gdpPercap_log_z | country), data = gapminder)

# Unscale both the GDP coefficient and the GDP random variance coefficient
tidy(model_gdp_country_only) %>% 
  mutate(conf.low = estimate - 1.96 * std.error,
         conf.high = estimate + 1.96 * std.error,
         across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term %in% c("gdpPercap_z", "sd__gdpPercap_z"), 
                        (. / gdp_sd) * 1000, .)))

tidy(model_gdp_country_only_log) %>% 
  mutate(conf.low = estimate - 1.96 * std.error,
         conf.high = estimate + 1.96 * std.error,
         across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term %in% c("gdpPercap_log_z", "sd__gdpPercap_log_z"), 
                        (. / gdp_log_sd), .)))

ranef(model_gdp_country_only)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  clean_names() %>% 
  # Unscale the GDP offsets
  mutate(gdp_percap_z = gdp_percap_z / gdp_sd * 1000)

coef(model_gdp_country_only)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  clean_names() %>% 
  # Unscale the GDP offsets
  mutate(gdp_percap_z = gdp_percap_z / gdp_sd * 1000)

# # doesn't seem to work
# ame_model_gdp_country_only <- model_gdp_country_only %>%
#   emtrends(~ country,
#            var = "gdpPercap_z",
#            at = list(year = 0, country = countries$country),
#            epred = TRUE, re_formula = NULL)
# 

# the rest of it is bayesian so omitted

# Each country gets its own intercept and GDP and year slopes -------
model_gdp_country_year <- lmer(lifeExp ~ gdpPercap_z + year + (1 + gdpPercap_z + year | country), data = gapminder)

model_gdp_country_year_log <- lmer(lifeExp ~ gdpPercap_log_z + year + (1 + gdpPercap_log_z + year | country),
  data = gapminder)

# Unscale both the GDP coefficient and the GDP random variance coefficient
tidy(model_gdp_country_year) %>% 
  mutate(conf.low = estimate - 1.96 * std.error,
         conf.high = estimate + 1.96 * std.error,
         across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term %in% c("gdpPercap_z", "sd__gdpPercap_z"), 
                        (. / gdp_sd) * 1000, .)))

tidy(model_gdp_country_year_log) %>% 
  mutate(conf.low = estimate - 1.96 * std.error,
         conf.high = estimate + 1.96 * std.error,
         across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term %in% c("gdpPercap_log_z", "sd__gdpPercap_log_z"), 
                        (. / gdp_log_sd), .)))

coef(model_gdp_country_year)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  clean_names() %>% 
  # Unscale the GDP offsets
  mutate(gdp_percap_z = gdp_percap_z / gdp_sd * 1000)

# # doesn't seem to work
# ame_model_gdp_country_year <- model_gdp_country_year %>% 
#   emtrends(~ year + country,
#            var = "gdpPercap_z",
#            at = list(year = 0, country = countries$country),
#            epred = TRUE, re_formula = NULL)


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

# the rest of it is bayesian

model_gdp_country_year_year <- lmer(lifeExp ~ gdpPercap_log_z + year + (1 + gdpPercap_log_z | year) + 
       (1 + gdpPercap_log_z + year | country),
  data = gapminder)


# Unscale both the GDP coefficient and the GDP random variance coefficient
tidy(model_gdp_country_year_year) %>% 
  mutate(conf.low = estimate - 1.96 * std.error,
         conf.high = estimate + 1.96 * std.error,
         across(c(estimate, std.error, conf.low, conf.high),
                ~ifelse(term %in% c("gdpPercap_log_z", "sd__gdpPercap_log_z"), 
                        (. / gdp_log_sd), .)))

coef(model_gdp_country_year)$country %>%
  as_tibble(rownames = "country") %>% 
  filter(country %in% countries$country) %>% 
  clean_names() %>%  
  # Unscale the GDP offsets
  mutate(gdp_percap_z = gdp_percap_z / gdp_sd * 1000)

# # doesn't seem to work
# # Different slopes in each year!
# model_gdp_country_year_year %>% 
#   emtrends(~ year + country,
#            var = "gdpPercap_log_z",
#            at = list(year = c(0, 5), country = countries$country),
#            epred = TRUE, re_formula = NULL)

# the rest of it has been removed as it is bayesian