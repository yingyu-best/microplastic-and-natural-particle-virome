#####indicator fisher test####
fisher<-read.csv('fisher.indicator.test.csv')

library(dplyr)
library(tidyr)
#  MP/NP ==> Case/Ctrl
out <- fisher %>%
  mutate(group = if_else(category == "MP", "Case", "Ctrl")) %>%
  select(func, season, group, carried, noncarried) %>%
  pivot_wider(
    names_from = group,
    values_from = c(carried, noncarried),
    values_fill = 0
  ) %>%
  transmute(
    func, season,
    Case_yes = carried_Case,
    Case_no  = noncarried_Case,
    Ctrl_yes = carried_Ctrl,
    Ctrl_no  = noncarried_Ctrl
  )

out2 <- out %>%
  tidyr::unite(Category, func, season, sep = " | ", remove = FALSE)


# Haldane–Anscombe correction
recalc <- out %>%
  rowwise() %>%
  do({
    a <- .$Case_yes; b <- .$Case_no; c <- .$Ctrl_yes; d <- .$Ctrl_no
    
    # 原始 Fisher p 值（不做 +0.5 修正）
    p <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))$p.value
    
    # Haldane–Anscombe 连续性校正：四格都 +0.5
    a2 <- a + 0.5; b2 <- b + 0.5; c2 <- c + 0.5; d2 <- d + 0.5
    
    # 校正后的 OR 和 Wald 型 95%CI
    OR   <- (a2 * d2) / (b2 * c2)
    se   <- sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
    lcl  <- exp(log(OR) - 1.96 * se)
    ucl  <- exp(log(OR) + 1.96 * se)
    
    tibble(OR = OR, CI_low = lcl, CI_high = ucl, p_value = p)
  }) %>%
  ungroup()

res2 <- out %>%
  bind_cols(recalc) %>%
  unite(Category, func, season, sep = " | ", remove = FALSE)

res_plot <- res2 %>%
  mutate(
    Category = factor(Category, levels = Category[order(OR)]),
    neglog10p = -log10(p_value),
    sig = p_value < 0.05
  )

res_plot2 <- res_plot %>%
  mutate(
    log10OR   = log10(OR),
    log10_low = log10(CI_low),
    log10_high= log10(CI_high),
    side = case_when(
      log10OR > 0 ~ "MP",   
      log10OR < 0 ~ "NP",   
      TRUE       ~ "= 0"
    ),
    neglog10p = -log10(pmax(p_value, .Machine$double.xmin))
  )

#visualization：
library(dplyr)
library(ggplot2)
res_plot2$func<-factor(res_plot2$func,levels = unique(res_plot2$func))
ggplot(res_plot2, aes(y = func, x = OR)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high,color=side), height = 0.15) +
  geom_point(aes(size = 5, fill =side,color=side,shape=sig),alpha=0.9) +
  scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1),
                     name = "p < 0.05",
                     labels = c(`TRUE` = "Yes", `FALSE` = "No")) +
  scale_fill_manual(values = c(  "#F8766D", "#00BA38")) +
  scale_color_manual(values = c(  "#F8766D", "#00BA38")) +
  scale_x_log10() +
  scale_size_continuous(name = "-log10(p)") +
  labs(x = "Odds Ratio (log scale)", y = NULL, title = "Fisher's Exact Test – Forest Plot") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )+facet_grid(~ season, scales = "free_y")
