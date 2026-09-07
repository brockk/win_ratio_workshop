
# Load ----
library(tidyverse)
library(glue)
library(ggpattern)
library(gridExtra)

# Cumulative analysis ----
run_cumulative_analysis <- function(data_package, priority, file_name = NULL) {
  if(!is.null(file_name)) {
    # Analysis might already exist
    if(file.exists(file_name)) {
      # Analysis does exist! Use it instead of re-running
      x <- readRDS(file_name)
      # Immediately return
      return(x)
    }
  }
  
  # Cached analysis could not be found. Run
  n_arm <- count(data_package$df, arm)
  testthat::expect_equal(nrow(n_arm), 2)
  
  x <- map(
    seq_along(priority),
    function(i) {
      
      # Attention. In WINS, priority refers to the position of columns whose 
      # names start with Y_, rather than columns names.
      # I.e. data = tibble(Y_3, Y_1, Y_2) and priority = 1:2 means
      # Y_3 > Y_1
      # and not
      # Y_1 > Y_2.
      # That in itself is not too much hassle. But there's more.
      # win.stat also requires that ep_type is the same length as priority.
      # This means you cannot just provide your full data-frame with a 
      # full-length ep_type vector and specify just priority = 1:2 to use the 
      # first two outcomes. This will error because priority and ep_type have
      # different lengths.
      # A workaround is to reorder the columns in the data-frame, subset ep_type
      # and tau congruently, and specify priority = 1:i, where i is the depth in
      # the cumulative analyse. Frustrating.
      
      these_indices <- seq(1, i)
      this_priority <- priority[these_indices]
      
      x <- win.stat(
        data = data_package$df %>% 
          # Note - this reorders outcomes so ep_type and tau must reorder too
          select(id, arm, all_of(paste0("Y_", this_priority))),
        ep_type = data_package$ep_type[this_priority],
        arm.name = c("T", "C"),
        priority = these_indices,
        tau = data_package$tau[this_priority],
        weight = "unstratified",
        alpha = 0.05,
        digit = 3,
        pvalue = "two-sided",
        summary.print = FALSE
      )
      
      # Manually calculate the number of comparisons because win.stat does not
      x$N <- prod(n_arm$n)
      # And record priority vector
      x$priority <- this_priority
      return(x)
    }
  )
  names(x) <- data_package$labels[priority]
  
  # Save if desired
  if(!is.null(file_name)) {
    saveRDS(x, file_name)
  }
  
  # Return
  return(x)
}

# Helpers ----
print_p <- function(p, dp = 3) {
  if(p < 10^-dp) {
    return(glue("p < { 10^-dp }"))
  } else {
    return(
      glue("p = { format(round(p, dp), nsmall = dp, scientific = FALSE) }")
    )
  }
}
# E.g.
# print_p(p = 0.0005)
# print_p(p = 0.0005, dp = 2)
# print_p(p = 0.0005, dp = 4)
# print_p(p = 0.0214563125, dp = 10)
# print_p(p = 0.0214563125, dp = 3)

wr_string_gpc <- function(x, p_value = TRUE, dp = 2, p_dp = 3) {
  # E.g. "Win ratio = 1.12, 95% CI:[1.03, 1.22], p = 0.008"
  wr <- format(
    x$Win_statistic$Win_Ratio[["WR"]],
    digits = dp,
    nsmall = dp,
    scientific = FALSE
  )
  wr_l <- format(
    x$Win_statistic$Win_Ratio[["WR_L"]],
    digits = dp,
    nsmall = dp,
    scientific = FALSE
  )
  wr_u <- format(
    x$Win_statistic$Win_Ratio[["WR_U"]],
    digits = dp,
    nsmall = dp,
    scientific = FALSE
  )
  template <-
    "Win ratio = { wr }, 95% CI:[{ wr_l }, { wr_u }]"
  msg <- glue(template)

  if(p_value) {
    p <- print_p(x$p_value$pvalue_WR, dp = p_dp)
    msg <- paste0(msg, ", ", p)
  }

  return(msg)
}
# E.g.
# wr_string_gpc(adrenal_1234_1d)
# wr_string_gpc(adrenal_1234_1d, dp = 3)
# wr_string_gpc(adrenal_1234_1d, p_value = FALSE)
# wr_string_gpc(adrenal_1234_1d, p_value = TRUE, p_dp = 8)

wo_string_gpc <- function(x, p_value = TRUE, dp = 2, p_dp = 3) {
  # E.g. "Win odds = 1.12, 95% CI:[1.03, 1.22], p = 0.008"
  wo <- format(
    x$Win_statistic$Win_Odds[["WO"]],
    digits = dp,
    nsmall = dp,
    scientific = FALSE
  )
  wo_l <- format(
    x$Win_statistic$Win_Odds[["WO_L"]],
    digits = dp,
    nsmall = dp,
    scientific = FALSE
  )
  wo_u <- format(
    x$Win_statistic$Win_Odds[["WO_U"]],
    digits = dp,
    nsmall = dp,
    scientific = FALSE
  )
  template <-
    "Win odds = { wo }, 95% CI:[{ wo_l }, { wo_u }]"
  msg <- glue(template)

  if(p_value) {
    p <- print_p(x$p_value$pvalue_WO, dp = p_dp)
    msg <- paste0(msg, ", ", p)
  }

  return(msg)
}
# E.g.
# wo_string_gpc(adrenal_1234_1d)
# wo_string_gpc(adrenal_1234_1d, dp = 3)
# wo_string_gpc(adrenal_1234_1d, p_value = FALSE)
# wo_string_gpc(adrenal_1234_1d, p_value = TRUE, p_dp = 8)

nb_string_gpc <- function(x, p_value = TRUE, dp = 1, p_dp = 3) {
  # E.g. "Net benefit = 1.9%, 95% CI:[0.5%, 3.4%], p = 0.008"
  nb <- format(
    100 * x$Win_statistic$Net_Benefit[["NB"]],
    digits = dp,
    nsmall = dp,
    scientific = FALSE
  )
  nb_l <- format(
    100 * x$Win_statistic$Net_Benefit[["NB_L"]],
    digits = dp,
    nsmall = dp,
    scientific = FALSE
  )
  nb_u <- format(
    100 * x$Win_statistic$Net_Benefit[["NB_U"]],
    digits = dp,
    nsmall = dp,
    scientific = FALSE
  )
  template <-
    "Net benefit = { nb }%, 95% CI:[{ nb_l }%, { nb_u }%]"
  msg <- glue(template)

  if(p_value) {
    p <- print_p(x$p_value$pvalue_NB, dp = p_dp)
    msg <- paste0(msg, ", ", p)
  }

  return(msg)
}
# E.g.
# nb_string_gpc(adrenal_1234_1d)
# nb_string_gpc(adrenal_1234_1d, dp = 3)
# nb_string_gpc(adrenal_1234_1d, p_value = FALSE)
# nb_string_gpc(adrenal_1234_1d, p_value = TRUE, p_dp = 8)

# Nice tabular outputs
as_tibble_gpc <- function(x) {
  
  summ_df <-
    x$summary_ep %>% 
    imap(
      ~ as_tibble(.x) %>% 
        mutate(series = .y)
    ) %>% 
    reduce(bind_rows) %>% 
    mutate(
      arm = str_extract(series, "Trt|Con"),
      Tier = str_extract(series, "[0-9]*$")
    ) %>% 
    select(Count, arm, Tier) %>% 
    pivot_wider(names_from = arm, values_from = Count, values_fill = 0) %>% 
    mutate(
      Ties = 0,
      N = x$N,
      Ties = N - cumsum(Trt) - cumsum(Con)
    )
  
  bind_rows(
    summ_df %>% 
      mutate(
        WR = Trt / Con,
        WO = (Trt + 0.5 * Ties) / (Con + 0.5 * Ties),
        NB = (Trt - Con) / N,
        NNT = 1 / NB,
        `Trt%` = Trt / N,
        `Con%` = Con / N,
        `Ties%` = Ties / N,
        NetDecisions = abs(Trt - Con),
        Contribution = NetDecisions / sum(NetDecisions)
      ),
    summ_df %>% 
      summarise(
        Tier = "Total",
        Trt = sum(Trt, na.rm = TRUE),
        Con = sum(Con, na.rm = TRUE),
        Ties = tail(Ties, 1),
        N = x$N
      ) %>% 
      mutate(
        WR = Trt / Con,
        WO = (Trt + 0.5 * Ties) / (Con + 0.5 * Ties),
        NB = (Trt - Con) / N,
        NNT = 1 / NB,
        `Trt%` = Trt / N,
        `Con%` = Con / N,
        `Ties%` = Ties / N,
        NetDecisions = abs(Trt - Con),
        Contribution = 1
      )
  )
}

as_tibble_row_gpc <- function(x) {
  tibble(
    WR = x$Win_statistic$Win_Ratio[["WR"]],
    WR2.5 = x$Win_statistic$Win_Ratio[["WR_L"]],
    WR97.5 = x$Win_statistic$Win_Ratio[["WR_U"]],
    p = x$p_value$pvalue_WR,
    NB = x$Win_statistic$Net_Benefit[["NB"]],
    NB2.5 = x$Win_statistic$Net_Benefit[["NB_L"]],
    NB97.5 = x$Win_statistic$Net_Benefit[["NB_U"]],
    Ties = 1 - x$Win_prop$P_trt - x$Win_prop$P_con
  )
}
gpc_to_table_row <- as_tibble_row_gpc

# Visualisation functions ----
wr_breaks <- function(lims) {
  # lims is a numeric vector c(min, max)

  lo <- floor(lims[1] * 10) / 10
  hi <- ceiling(lims[2] * 10) / 10
  symmetrical_delta <- max(abs(1 - lo), abs(1 - hi)) # Force 1.0 as centre

  seq(
    1 - symmetrical_delta,
    1 + symmetrical_delta,
    by = 0.5
  )
}
# wr_breaks(c(0.92, 1.15))
# wr_breaks(c(0.92, 1.25))

# Main two visualisation techniques used in ADVANCE re-analysis.
# Written by Severine and tweaked by Kristian
wr_cumulative_panel_plot <- function(
    x, 
    font_size = 13,
    panel_widths = c(5, 3, 2, 1.5),
    control_fill = "#56B4E9",
    treatment_fill = "#0072B2"
) {
  
  df1 <-
    imap(
      x, 
      ~ as_tibble_row_gpc(.x) %>% 
        mutate(LABEL = .y)
    ) %>% 
    reduce(bind_rows) %>% 
    rename(
      WR_L = `WR2.5`,
      WR_U = `WR97.5`,
      NB_L = `NB2.5`,
      NB_U = `NB97.5`,
      PVAL = p
    ) %>% 
    mutate(
      NB_DENOM = max(NB), 
      NB_PCT = NB / NB_DENOM * 100
    )
  
  df2 <- 
    imap(
      x, 
      ~ as_tibble_gpc(.x) %>% 
        mutate(LABEL = .y) %>% 
        filter(Tier == "Total")
    ) %>% 
    reduce(bind_rows) %>% 
    rename(
      P_trt = `Trt%`,
      P_con = `Con%`
    ) %>% 
    select(P_trt, P_con, LABEL)
  
  df3 <- 
    df2 %>%
    pivot_longer(
      cols = c("P_trt", "P_con"),
      names_to = c("Tmt"),
      values_to = "PROP"
    ) %>%
    mutate(
      Treatment = case_when(
        Tmt == "P_trt" ~ "Treatment", 
        Tmt == "P_con" ~ "Control", 
        .default = NA
      )
    ) %>%
    group_by(LABEL) %>%
    mutate(P_ties = 1 - sum(PROP))
  
  
  level_order <- rev(names(x))
  
  p1 <- ggplot(
    data = df3, 
    aes(x = factor(LABEL, level = level_order), y = PROP * 100)
  ) +
    geom_bar_pattern(
      aes(
        pattern_fill = Treatment,
        pattern = Treatment,
        fill = Treatment
      ),
      colour = 'white',
      pattern_angle = 40,
      pattern_density = 0.01,
      pattern_spacing = 0.015,
      pattern_color = "white",
      stat = "identity",
      width = 0.5,
      position = position_dodge()
    ) +
    coord_flip() +
    geom_vline(
      xintercept = 0.6 + length(x),
      linewidth = 2
    ) +
    geom_text(aes(
      y = 0,
      hjust = 1.1,
      label = paste0(round(P_ties * 100, 1), "%*")
    )) +
    labs(
      x = NULL, 
      title = "Percentage of wins", 
      caption = "*: Percentage ties."
    ) +
    geom_text(
      aes(label = paste0(round(PROP * 100, 1), "%"), group = Treatment),
      hjust = -0.2,
      # size = 4.5,
      position = position_dodge(width = .5)
    ) +
    scale_y_continuous(
      limits = function(lims) c(-4, lims[2] + 4)
    ) +
    scale_fill_manual(values = c(
      "Control" = control_fill,
      "Treatment" = treatment_fill
    )) +
    scale_pattern_discrete(choices = c("stripe", "none")) +
    guides(
      fill = guide_legend(
        override.aes = list(pattern = c("none", "stripe")), 
        reverse = TRUE
      ),
      pattern = "none",
      # suppress separate pattern legend
      pattern_fill = "none"  # suppress separate pattern_fill legend
    ) +
    # theme_classic(base_size = 14) +
    theme_classic(base_size = font_size) +
    theme(
      plot.margin = unit(c(2, 1, 1, 1), "lines"),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.y = element_line(colour = "white"),
      axis.line.x = element_line(colour = "white", linewidth = 0.6),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = font_size),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = c(0.9, 0.7),
      # legend.position = c(0.9, 0.9),
      legend.title = element_blank(),
      legend.text = element_text(size = font_size)
    )
  # p1
  
  
  p2 <- ggplot(
    df1,
    aes(
      x = WR,
      y = factor(LABEL, level = level_order),
      shape = LABEL
    )
  ) +
    #Add dot plot and error bars
    geom_errorbar(aes(xmin = WR_L, xmax = WR_U), width = 0.25) +
    geom_point(size = 2.5) +
    ggtitle("Control better       Treatment better") + # Attention!
    # Add a line above graph
    geom_hline(
      yintercept = 0.6 + length(x), 
      linewidth = 2
    ) +
    #Add a reference dashed line at 0
    geom_vline(xintercept = 1, linetype = "longdash") +
    labs(x = "", y = "Component") +
    scale_x_continuous(
      breaks = wr_breaks
    ) +
    scale_shape_manual(values = c(18, 18, 18, 18, 18)) +
    theme_gray(base_size = font_size) +
    # Remove legend and remove y-axis line and ticks
    theme(
      legend.position = "none",
      plot.margin = unit(c(2, 1, 1, 1), "lines"),
      plot.title = element_text(hjust = 0.5),
      axis.line.x = element_line(linewidth = 0.6),
      axis.ticks.length = unit(0.3, "cm"),
      axis.text.y  = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y  = element_blank()
    )
  # p2
  
  
  t1 <- ggplot(data = df1) +
    geom_text(aes(
      y = factor(LABEL, level = level_order),
      x = 1,
      label = paste0(
        format(WR, digits = 2, nsmall = 2),
        " (",
        format(WR_L, digits = 2, nsmall = 2),
        "; ",
        format(WR_U, digits = 2, nsmall = 2),
        ")  ",
        if_else(PVAL >= 0.001, paste("p =", round(PVAL, digits = 3)), "p < 0.001")
      )
    ), vjust = 0) +
    #Add a line above graph
    geom_hline(
      yintercept = 0.6 + length(x), 
      linewidth = 2
    ) +
    ggtitle("Win ratio (95% CI)") +
    xlab("  ") +
    theme_classic(base_size = font_size) +
    theme(
      plot.margin = unit(c(2, 1, 1, 1), "lines"),
      axis.line.y = element_blank(),
      axis.line.x = element_line(color = "white"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(color = "white"),
      axis.ticks.length = unit(0.3, "cm"),
      axis.title.y = element_blank(),
      axis.text.x = element_text(color = "white"),
      plot.title = element_text(hjust = 0.5)
    )
  # t1
  
  
  t2 <- ggplot(data = df1) +
    geom_text(aes(
      y = factor(LABEL, level = level_order),
      x = 1,
      label = paste0(format(
        NB_PCT, digits = 1, nsmall = 1
      )),
      vjust = 0
    )) +
    # Add a line above graph
    geom_hline(
      yintercept = 0.6 + length(x), 
      linewidth = 2
    ) +
    ggtitle("% of decisions") +
    xlab("  ") +
    theme_classic(base_size = font_size) +
    theme(
      plot.margin = unit(c(2, 1, 1, 1), "lines"),
      axis.line.y = element_blank(),
      axis.line.x = element_line(color = "white"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(color = "white"),
      axis.ticks.length = unit(0.3, "cm"),
      axis.title.y = element_blank(),
      axis.text.x = element_text(color = "white"),
      plot.title = element_text(hjust = 0.5)
    )
  # t2
  
  # Put the individual components of the forest plot together
  fplt <- grid.arrange(p1, p2, t1, t2, widths = panel_widths)
  return(fplt)
}

# Written by Severine and tweaked by Kristian
wr_cumulative_flow_plot <- function(
    x,
    plot_title,
    fill = "#0A86C8"
) {
  
  # The final analysis is the last in the list of cumulative analyses
  final_analysis <- x[[length(x)]]
  
  df2 <- 
    imap(
      x, 
      ~ as_tibble_gpc(.x) %>% 
        mutate(LABEL = .y) %>% 
        filter(Tier == "Total")
    ) %>% 
    reduce(bind_rows) %>% 
    rename(
      P_trt = `Trt%`,
      P_con = `Con%`
    ) %>% 
    select(P_trt, P_con, LABEL)
  
  df3 <- 
    df2 %>%
    pivot_longer(
      cols = c("P_trt", "P_con"),
      names_to = c("Tmt"),
      values_to = "PROP"
    ) %>%
    mutate(
      Treatment = case_when(
        Tmt == "P_trt" ~ "Treatment", 
        Tmt == "P_con" ~ "Control", 
        .default = NA
      )
    ) %>%
    group_by(LABEL) %>%
    mutate(P_ties = 1 - sum(PROP))
  
  overall_wide_pct <- 
    as_tibble_gpc(final_analysis) %>% 
    filter(Tier == "Total") %>% 
    rename(
      P_trt = `Trt%`,
      P_con = `Con%`
    ) %>% 
    select(P_trt, P_con) %>% 
    mutate(LABEL = "Overall") %>% 
    pivot_longer(
      cols = c("P_trt", "P_con"),
      names_to = c("Tmt"),
      values_to = "PROP"
    ) %>%
    mutate(
      Treatment = case_when(
        Tmt == "P_trt" ~ "Treatment", 
        Tmt == "P_con" ~ "Control", 
        .default = NA
      )
    ) %>%
    group_by(LABEL) %>%
    mutate(P_ties = 1 - sum(PROP))
  
  
  df_plot <- 
    # bind_rows(BP_wide_pct, overall_wide_pct) %>% 
    df3 %>% 
    
    # 1) keep the hierarchy order (use your row id if you have it)
    # If you have an explicit row number column (like the 1..10 shown), use it:
    # arrange(row_id) %>%
    # arrange(LABEL, Treatment) %>%  # fallback if no row id; better to use row_id
    
    # 2) clean endpoint names
    # mutate(endpoint = str_remove(LABEL, "^\\+\\s*")) %>%
    rename(endpoint = LABEL) %>% 
    
    # 3) wide format: one row per endpoint with cum values by arm
    select(endpoint, Treatment, PROP, P_ties) %>%
    # distinct() %>%
    pivot_wider(
      names_from  = Treatment,
      values_from = PROP,
      names_prefix = "cum_"
    ) %>%
    # P_ties is duplicated across Treatment/Control rows -> keep one
    group_by(endpoint) %>%
    summarise(
      cum_Treatment = first(cum_Treatment),
      cum_Control   = first(cum_Control),
      ties_cum      = first(P_ties),
      .groups = "drop"
    ) %>%
    # 4) enforce the intended hierarchy order
    # If your table is already in hierarchy order, you can supply it explicitly:
    mutate(endpoint = factor(
      endpoint,
      # levels = c(names(x), "Overall")
      levels = names(x)
    )) %>%
    arrange(endpoint) %>%
    
    # 5) cumulative -> incremental + conditional ties
    mutate(
      wins_treat = cum_Treatment - lag(cum_Treatment, default = 0),
      wins_ctrl  = cum_Control   - lag(cum_Control, default = 0),
      # ties_step  = ties_cum / lag(ties_cum, default = 1),
      
      # Optional: win difference at each step (incremental)
      win_diff   = wins_treat - wins_ctrl
    ) %>%
    
    bind_rows(
      as_tibble_gpc(final_analysis) %>% 
        filter(Tier == "Total") %>% 
        mutate(endpoint = "Overall") %>% 
        select(
          endpoint,
          # cum_Treatment = `Trt%`,
          wins_treat = `Trt%`,
          # cum_Control = `Con%`,
          wins_ctrl = `Con%`,
          win_diff = NB,
          # ties_step = `Ties%`,
          ties_cum = `Ties%`,
        )
    ) %>% 
    
    # 6) convert to percent (like your plot)
    mutate(across(
      c(
        wins_treat, wins_ctrl, 
        # ties_step, 
        win_diff, ties_cum
      ),
      ~ .x * 100
    )) %>% 
    mutate(
      # y = rev(row_number())
      endpoint = factor(
        endpoint,
        levels = rev(c(names(x), "Overall"))
      )
    )
  
  # ---- Layout coordinates (x positions for the 3 columns) ----
  x_left  <- 1
  x_mid   <- 2
  x_right <- 3
  x_diff  <- 4.2   # right-hand text column
  
  # Box geometry
  box_w <- 0.72
  box_h <- 0.55
  
  # Helper to build box rectangles for each column
  boxes <- 
    df_plot %>%
    mutate(y = as.numeric(endpoint)) %>%
    pivot_longer(
      # cols = c(wins_treat, ties_step, wins_ctrl),
      cols = c(wins_treat, ties_cum, wins_ctrl),
      names_to = "type",
      values_to = "pct"
    ) %>%
    mutate(
      x = case_when(
        type == "wins_treat" ~ x_left,
        # type == "ties_step"  ~ x_mid,
        type == "ties_cum"  ~ x_mid,
        type == "wins_ctrl"  ~ x_right
      ),
      xmin = x - box_w / 2,
      xmax = x + box_w / 2,
      ymin = y - box_h / 2,
      ymax = y + box_h / 2
    ) %>% 
    mutate(label = sprintf("%.1f%%", pct))
  
  # Make sure your y values follow the factor order you want top->bottom
  df_y <- 
    df_plot %>%
    mutate(y = as.numeric(endpoint)) %>%
    arrange(desc(y))
  
  segs_down3 <- 
    bind_rows(
      # to left box center below
      df_y %>%
        mutate(endpoint_next = lead(endpoint), y_next = lead(y)) %>%
        filter(!is.na(y_next)) %>%
        filter(endpoint != "Overall", endpoint_next != "Overall") %>%
        transmute(
          x = x_mid,
          y = y - box_h / 2,
          xend = x_left,
          yend = y_next  + box_h / 2,
          kind = "to_left_center"
        ),
      
      # to middle box top below
      df_y %>%
        mutate(endpoint_next = lead(endpoint), y_next = lead(y)) %>%
        filter(!is.na(y_next)) %>%
        filter(endpoint != "Overall", endpoint_next != "Overall") %>%
        transmute(
          x = x_mid,
          y = y - box_h / 2,
          xend = x_mid,
          yend = y_next + box_h / 2,
          kind = "to_mid_top"
        ),
      
      # to right box top below
      df_y %>%
        mutate(endpoint_next = lead(endpoint), y_next = lead(y)) %>%
        filter(!is.na(y_next)) %>%
        filter(endpoint != "Overall", endpoint_next != "Overall") %>%
        transmute(
          x = x_mid,
          y = y - box_h / 2,
          xend = x_right,
          yend = y_next + box_h / 2,
          kind = "to_right_top"
        )
    )
  
  
  y_overall <- 
    df_y %>% 
    filter(endpoint == "Overall") %>% 
    pull(y)
  y_prev <- 
    df_y %>%
    filter(y > y_overall) %>%
    summarise(y_prev = min(y)) %>%
    pull(y_prev)
  y_sep <- (y_overall + y_prev) / 2
  
  
  # ---- Plot ----
  p <- ggplot() +
    # Connectors
    
    # --- NEW: 3 downward connectors (per row) ---
    geom_segment(
      data = segs_down3,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey50"
    ) +
    
    # Boxes: Treatment wins (solid blue)
    geom_rect(
      data = boxes %>% filter(type == "wins_treat"),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      fill = fill,
      colour = "black",
      linewidth = 0.3
    ) +
    
    # Boxes: Ties (solid grey)
    # geom_rect(
    #   data = boxes %>% filter(type == "ties_step"),
    #   aes(
    #     xmin = xmin,
    #     xmax = xmax,
    #     ymin = ymin,
    #     ymax = ymax
    #   ),
    #   fill = "#BFBFBF",
    #   colour = "black",
    #   linewidth = 0.3
    # ) +
    geom_rect(
      data = boxes %>% filter(type == "ties_cum"),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      fill = "#BFBFBF",
      colour = "black",
      linewidth = 0.3
    ) +
    
    # Boxes: Control wins (hatched light-blue)
    ggpattern::geom_rect_pattern(
      data = boxes %>% dplyr::filter(type == "wins_ctrl"),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      # rectangle styling
      fill = "white",
      colour = fill,
      linewidth = 0.5,
      
      # pattern styling (45° diagonal stripes)
      pattern = "stripe",
      pattern_angle = 135,
      # <-- 45 degrees
      pattern_fill = fill,
      # stripe color
      pattern_colour = fill,
      pattern_alpha = 1,
      pattern_density = 0.05,
      # stripe coverage (0–1). Higher = more stripes
      pattern_spacing = 0.02,
      # distance between stripes (smaller = denser)
      pattern_size = 0.2         # stripe thickness
    ) +
    
    # Percent labels
    geom_text(
      data = boxes,
      aes(x = x, y = y, label = label),
      size = 6,
      colour = "black",
      fontface = "plain"
    ) +
    
    # Left-side endpoint names
    geom_text(
      data = df_plot %>% mutate(y = as.numeric(endpoint)),
      aes(
        x = 0.15,
        y = y,
        label = as.character(endpoint)
      ),
      hjust = 0,
      size = 6
    ) +
    
    # Column headings
    annotate(
      "text",
      x = x_left,
      y = max(as.numeric(df_plot$endpoint)) + 0.9,
      label = "Treatment wins",
      fontface = "bold",
      size = 6
    ) +
    annotate(
      "text",
      x = x_mid,
      y = max(as.numeric(df_plot$endpoint)) + 0.9,
      label = "Ties",
      fontface = "bold",
      size = 6
    ) +
    annotate(
      "text",
      x = x_right,
      y = max(as.numeric(df_plot$endpoint)) + 0.9,
      label = "Control wins",
      fontface = "bold",
      size = 6
    ) +
    annotate(
      "text",
      x = x_diff,
      y = max(as.numeric(df_plot$endpoint)) + 0.9,
      label = "Win difference",
      fontface = "bold",
      size = 6
    ) +
    
    # Win difference values
    geom_text(
      data = df_plot %>% mutate(y = as.numeric(endpoint)),
      aes(
        x = x_diff,
        y = y,
        label = sprintf("%.1f%%", win_diff)
      ),
      size = 6
    ) +
    
    # Plus signs between lines (to mimic the stacked + on the right)
    annotate(
      "text",
      x = x_diff,
      # y = c(2.5, 3.5, 4.5), # TODO
      y = 0.5 + 
        df_y %>% 
        filter(y < max(y) & y > min(y)) %>%  # Exclude first and last
        pull(y),
      label = "+",
      size = 6
    ) + 
    
    # Title
    annotate(
      "text",
      x = 0.15,
      y = max(as.numeric(df_plot$endpoint)) + 1.7,
      label = plot_title,
      hjust = 0,
      size = 10,
      colour = fill,
      fontface = "plain"
    ) +
    
    # Top right stats box
    geom_label(
      data = data.frame(
        x = 3.2,
        y = max(as.numeric(df_plot$endpoint)) + 1.65,
        label = paste0(
          wr_string_gpc(final_analysis, p_value = TRUE), "\n",
          wo_string_gpc(final_analysis, p_value = FALSE), "\n",
          nb_string_gpc(final_analysis, p_value = FALSE)
        )
      ),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      fill = "grey92",
      linewidth = 0,
      size = 4
    ) +
    geom_hline(
      yintercept = y_sep,
      colour = "grey30",
      linewidth = 1
    ) +
    
    # Styling
    coord_cartesian(
      # xlim = c(0, 4.6),
      # ylim = c(0.5, 7.8), # TODO danger!
      clip = "off"
    ) +
    
    theme_void() +
    theme(plot.margin = margin(20, 20, 20, 20))
  
  return(p)
}
