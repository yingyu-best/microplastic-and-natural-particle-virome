###partial regression testing individual effects of environmental factors on viral diversity ####

library(lme4)
library(dplyr)
library(purrr)
library(tibble)
library(broom)
library(ggplot2)

#function
mk_form <- function(lhs, rhs, rand = NULL) {
  fixed <- if (length(rhs)) paste(rhs, collapse = " + ") else "1"
  if (!is.null(rand)) as.formula(paste(lhs, "~", fixed, "+", rand))
  else as.formula(paste(lhs, "~", fixed))
}

compute_partial_set <- function(dat, env_vars, dataset_name = "MP") {
  partial_r2_one <- function(var) {
    others <- setdiff(env_vars, var)
    
    mod_y <- lmer(mk_form("Evenness", others, "(1|site)"), data = dat, REML = TRUE)
    res_y <- resid(mod_y)
    
    # Residualize focal x against the other predictors (standard AV step)
    mod_x <- lm(mk_form(var, others), data = dat)
    res_x <- resid(mod_x)
    
    fit_av <- lm(res_y ~ res_x)
    s <- summary(fit_av)
    
    stats_row <- tibble(
      dataset    = dataset_name,
      variable   = var,
      partial_R2 = unname(s$r.squared),
      slope      = unname(coef(s)[2, "Estimate"]),
      t_value    = unname(coef(s)[2, "t value"]),
      p_value    = unname(coef(s)[2, "Pr(>|t|)"]),
      n          = length(res_y),
      df_resid   = unname(s$df[2])
    )
    
    av_df <- tibble(dataset = dataset_name, variable = var, res_x = res_x, res_y = res_y)
    list(stats = stats_row, av_df = av_df)
  }
  
  res_list <- map(env_vars, partial_r2_one)
  list(
    stats = bind_rows(map(res_list, "stats")),
    av_df = bind_rows(map(res_list, "av_df"))
  )
}

partial.input<-read.csv('partial.diversity.lifestyle.csv',row.names = 1)

# MP: rows 1:60 (lytic),and 121:180 (lysogenic) (includes totalMP)
partial.input<-data.frame(partial.input[,c(1:4,17)],scale(partial.input[,-c(1:4,17)]))
dat_MP <- data.frame(partial.input[1:60,]) #lytic
#dat_MP <- data.frame(partial.input[121:180,]) #lysogenic
env_MP <- c("Rainfall","Turbidity","Salinity","pH","TN","TP",
            "DO","COD","SO2","velocity","Temperature","totalMP")

# NP: rows 61:120 (lytic),and 181:240 (lysogenic)  (exclude totalMP)
dat_NP <- data.frame(partial.input[61:120, ])
#dat_NP <- data.frame(partial.input[181:240, ])
env_NP <- c("Rainfall","Turbidity","Salinity","pH","TN","TP",
            "DO","COD","SO2","velocity","Temperature")

# ---------- compute per set ----------
res_MP <- compute_partial_set(dat_MP, env_MP, dataset_name = "MP")
res_NP <- compute_partial_set(dat_NP, env_NP, dataset_name = "NP")

# ---------- merge ----------
partial_summary <- bind_rows(res_MP$stats, res_NP$stats)
av_df_all       <- bind_rows(res_MP$av_df, res_NP$av_df)

# ---------- sample groups & colors (applied within each dataset separately) ----------
av_df_all <- av_df_all %>%
  dplyr::group_by(dataset, variable) %>%
  dplyr::mutate(
    sample_id = dplyr::row_number(),
    group = dplyr::case_when(
      dataset == "MP" & sample_id <= 15 ~ "mp1",
      dataset == "MP" & sample_id <= 30 ~ "mp2",
      dataset == "MP" & sample_id <= 45 ~ "mp3",
      dataset == "MP"                   ~ "mp4",
      dataset == "NP" & sample_id <= 15 ~ "np1",
      dataset == "NP" & sample_id <= 30 ~ "np2",
      dataset == "NP" & sample_id <= 45 ~ "np3",
      TRUE                              ~ "np4"
    )
  ) %>%
  dplyr::ungroup()

group_cols <- c("mp1"   = "#993e36",
                "mp2"  = "#eb746a",
                "mp3"  = "#f6e2e1",
                "mp4"  = "#e9ada9",
                "np1"   = "#135e32",
                "np2"  = "#86c13a",
                "np3"  = "#dbebd1",
                "np4"  = "#bad79d")
sig_df <- partial_summary %>%
  transmute(
    dataset,
    variable,
    signif = ifelse(p_value <= 0.05, "p ≤ 0.05", "p > 0.05")
  )
# ---------- per-panel labels ----------
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))
labels_df <- partial_summary %>%
  transmute(
    dataset, variable,
    label = paste0("Partial R^2 = ", sprintf("%.3f", partial_R2),
                   "\nP ", fmt_p(p_value))
  )

# ---------- facet order (your custom order) ----------
var_order <- c("Temperature","Rainfall","Turbidity","velocity",
               "Salinity","pH","TN","TP","DO","COD","SO2","totalMP")

av_df_all  <- av_df_all  %>% mutate(variable = factor(variable, levels = var_order),
                                    dataset  = factor(dataset, levels = c("MP","NP")))
labels_df  <- labels_df  %>% mutate(variable = factor(variable, levels = var_order),
                                    dataset  = factor(dataset, levels = c("MP","NP")))
sig_df     <- sig_df     %>% mutate(variable = factor(variable, levels = var_order),
                                    dataset  = factor(dataset, levels = c("MP","NP")))

# ---------- colors for datasets (your request) ----------
dataset_cols <- c("MP" = "#FB8072FF", "NP" = "#18b327")

# ---------- plot: MP vs NP rows, variables as columns ----------
p <- ggplot(av_df_all, aes(res_x, res_y)) +
  # Points = sample groups
  geom_point(aes(color = group), alpha = 0.85, size = 2) +
  
  # Regression lines = dataset color + significance linetype
  geom_smooth(
    data = av_df_all %>%
      dplyr::left_join(sig_df[, c("dataset","variable","signif")],
                       by = c("dataset","variable")),
    aes(color = dataset,fill=dataset),alpha=0.15,
    method = "lm"#, se = FALSE
  ) +
  
  # Manual color scales
  scale_color_manual(values = c(group_cols, dataset_cols),
                     breaks = c(names(group_cols), names(dataset_cols)),
                     name = "Legend") +
  scale_fill_manual(values = dataset_cols, name = "Dataset")+
  scale_linetype_manual(values = c("p ≤ 0.05" = "solid",
                                   "p > 0.05" = "dashed"),
                        name = "Significance") +
  
  facet_wrap(~ variable, scales = "free") +
  
  # Labels (colored by dataset)
  geom_text(
    data = labels_df %>% filter(dataset == "MP"),
    aes(x = -Inf, y = Inf, label = label, color = dataset),
    hjust = -0.05, vjust = 1.1, inherit.aes = FALSE, size = 3
  ) +
  geom_text(
    data = labels_df %>% filter(dataset == "NP"),
    aes(x = -Inf, y = Inf, label = label, color = dataset),
    hjust = -0.05, vjust = 2.4, inherit.aes = FALSE, size = 3
  ) +
  
  labs(
    x = "Residuals of predictor ~ other predictors",
    y = "Residuals of indicator.rb ~ other predictors"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  theme(
    strip.text.x   = element_text(face = "bold"),
    plot.title     = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p)

ggsave("richness.lysogenic.scaled.partiallinear.pdf",width = 15,height = 8)