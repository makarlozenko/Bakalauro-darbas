install.packages("patchwork")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(purrr)

set.seed(123)

#-------------------------Simuliacijos parametrai-------------------------------

pd_values <- c(0.001, 0.005, 0.01, 0.02, 0.05) # PD reikšmės

n_values <- c(5000, 10000, 25000, 50000, 100000) # Portfelio dydis

rho_values <- c(0.01, 0.05, 0.10, 0.20, 0.25) # Koreliacijos reikšmės

alpha <- 0.05 # reikšmingumo lygmuo

B <- 10000 # Simuliacijų kiekis

# Kursim po 200 portfelių, kad greičiau veiktų skriptas
chunk_size <- 200

# Limituojam galimus intervalus iki - [0, 1]
clip_01 <- function(x) {
  pmin(1, pmax(0, x))
}

# -------------------Funkcijos pasikliovimo intervalamas------------------------

# Wald intervalas
wald_ci <- function(x, n, alpha = 0.05) {
  p_hat <- x / n
  z <- qnorm(1 - alpha / 2)

  lower <- clip_01(p_hat - z * sqrt(p_hat * (1 - p_hat) / n))
  upper <- clip_01(p_hat + z * sqrt(p_hat * (1 - p_hat) / n))
  
  data.frame(lower = lower, upper = upper)
}

# Clopper - Pearson intervalas
clopper_pearson_ci <- function(x, n, alpha = 0.05) {
  lower <- numeric(length(x))
  upper <- numeric(length(x))
  
  lower[x == 0] <- 0
  lower[x > 0] <- qbeta(alpha / 2, x[x > 0], n - x[x > 0] + 1)
  
  upper[x == n] <- 1
  upper[x < n] <- qbeta(1 - alpha / 2, x[x < n] + 1, n - x[x < n])
  
  data.frame(lower = lower, upper = upper)
}

# Wilson intervalas
wilson_ci <- function(x, n, alpha = 0.05) {
  p_hat <- x / n
  z <- qnorm(1 - alpha / 2)
  
  denominator <- 1 + z^2 / n
  
  center <- (p_hat + z^2 / (2 * n)) / denominator
  
  half_width <- (z / denominator) *
    sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
  
  lower <- clip_01(center - half_width)
  upper <- clip_01(center + half_width)
  
  data.frame(lower = lower, upper = upper)
}

# Agresti - Coull intervalas
agresti_coull_ci <- function(x, n, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  
  n_tilde <- n + z^2
  p_tilde <- (x + z^2 / 2) / n_tilde
  
  se <- sqrt(p_tilde * (1 - p_tilde) / n_tilde)
  
  lower <- clip_01(p_tilde - z * se)
  upper <- clip_01(p_tilde + z * se)
  
  data.frame(lower = lower, upper = upper)
}

# Visi naudojami intervalai
ci_methods <- list(
  Wald = wald_ci,
  Clopper_Pearson = clopper_pearson_ci,
  Wilson = wilson_ci,
  Agresti_Coull = agresti_coull_ci
)

#----------------------Nepriklausomų įvykių modelis-----------------------------

simulate_independent_counts <- function(B, n, p) {
  rbinom(n = B, size = n, prob = p)
}

#-------------------Priklausomų įvykių modelis Vašiček--------------------------

simulate_vasicek_counts <- function(B, n, p, rho, chunk_size = 200,
                                    arma_ar = c(0.5, -0.2),
                                    arma_ma = 0.3) {
  
  threshold <- qnorm(p)
  x <- numeric(B)
  
  start_indices <- seq(1, B, by = chunk_size)
  
  for (start in start_indices) {
    
    end <- min(start + chunk_size - 1, B)
    current_B <- end - start + 1
    
    # Bendras sisteminis rizikos veiksnys kiekvienam portfeliui
    y <- rnorm(current_B)
    
    # Individualus veiksnys kiekvienam klientui
    z <- matrix(
      rnorm(current_B * n),
      nrow = current_B,
      ncol = n
    )
    
    # Latentinis kintamasis
    asset_values <- sqrt(1 - rho) * z
    
    asset_values <- sweep(
      asset_values,
      MARGIN = 1,
      STATS = sqrt(rho) * y,
      FUN = "+"
    )
    
    # Tikrinama, ar įvyko įsipareigojimų nevykdymas
    defaults <- asset_values <= threshold
    
    # Defoltų skaičius kiekviename portfelyje
    x[start:end] <- rowSums(defaults)
  }
  
  x
}

#-----------------------------Generuojam defoltų kiekį--------------------------

simulate_default_counts <- function(B, n, p, rho = 0,model = "independent", chunk_size = 200) {
  
  if (model == "independent") {
    return(simulate_independent_counts(B = B, n = n, p = p))
  }
  if (model == "vasicek") {
    return(simulate_vasicek_counts( B = B, n = n, p = p,rho = rho, chunk_size = chunk_size))
  }
  stop("Nežinomas įvestas modelis")
}

#--------------------------Tikrinam vieną scenarijų-----------------------------

evaluate_scenario <- function(B, n, p, rho = 0, model = "independent", alpha = 0.05, chunk_size = 200) {
  # Gaunam defoltų kiekį
  x <- simulate_default_counts(B = B, n = n, p = p, rho = rho, model = model, chunk_size = chunk_size)
  # Gaunam PD įvertį
  p_hat <- x / n
  
  scenario_results <- list()
  
  for (method_name in names(ci_methods)) {
    # Einam per visus pasikliovimo intervalus
    ci <- ci_methods[[method_name]](x = x, n = n, alpha = alpha)
    
    interval_width <- ci$upper - ci$lower
    # Skaičiuojam palyginimo metrikas
    coverage <- mean(ci$lower <= p & p <= ci$upper)
    avg_width <- mean(interval_width)
    median_width <- median(interval_width)
    lower_zero_rate <- mean(ci$lower == 0)
    
    # Kuriam lentelę kiekvinam pasikliovimo intervalui
    scenario_results[[method_name]] <- data.frame(
      model = model, 
      method = method_name,
      p = p,
      n = n,
      rho = ifelse(model == "independent", 0, rho),
      B = B,
      coverage = coverage,
      avg_width = avg_width,
      median_width = median_width,
      lower_zero_rate = lower_zero_rate,
      mean_p_hat = mean(p_hat),
      sd_p_hat = sd(p_hat),
      mean_defaults = mean(x),
      sd_defaults = sd(x)
    )
  }
  # Sujungiam viską kartu
  do.call(rbind, scenario_results)
}

#------------------------------Testavimas---------------------------------------
# Nepriklausomas atvejis
test_independent <- evaluate_scenario(
  B = 1000,
  n = 1000,
  p = 0.01,
  model = "independent",
  alpha = alpha,
  chunk_size = chunk_size
)
# Priklausomas atvejis
test_vasicek <- evaluate_scenario(
  B = 1000,
  n = 1000,
  p = 0.01,
  rho = 0.10,
  model = "vasicek",
  alpha = alpha,
  chunk_size = chunk_size
)
# Rezultatai
print(test_independent)
print(test_vasicek)

#----------------------Pradedame Monte Carlo simuliaciją------------------------

all_results <- list()
counter <- 1
# Kiek nusimato skirtingų scenarijų
total_scenarios <- length(pd_values) * length(n_values) *
  (1 + length(rho_values))

scenario_number <- 1

for (p in pd_values) {
  for (n in n_values) {
    # Kuriam vieną nepriklausomą scenarijų
    cat(
      "Scenarijus", scenario_number, "iš", total_scenarios,
      "| modelis = Nepriklausomas",
      "| p =", p,
      "| n =", n,
      "| rho = 0\n"
    )
    
    all_results[[counter]] <- evaluate_scenario(
      B = B,
      n = n,
      p = p,
      rho = 0,
      model = "independent",
      alpha = alpha,
      chunk_size = chunk_size
    )
    
    counter <- counter + 1
    scenario_number <- scenario_number + 1
    # Kiekvienai koreliacijai kuriam priklausomą scenarijų
    for (rho in rho_values) {
      
      cat(
        "Scenarijus", scenario_number, "iš", total_scenarios,
        "| modelis = vasicek",
        "| p =", p,
        "| n =", n,
        "| rho =", rho, "\n"
      )
      
      all_results[[counter]] <- evaluate_scenario(
        B = B,
        n = n,
        p = p,
        rho = rho,
        model = "vasicek",
        alpha = alpha,
        chunk_size = chunk_size
      )
      
      counter <- counter + 1
      scenario_number <- scenario_number + 1
    }
  }
}

results <- do.call(rbind, all_results)
# Simuliacijos rezultatai
write.csv(
  results,
  file = "results/ci_simulation_results.csv",
  row.names = FALSE
)
#----------------------------Rezultatų nuskaitymas------------------------------

results <- read.csv("results/ci_simulation_results.csv")

#--------------------------Duomenų apdorojimas----------------------------------

results_plot <- results %>%
  mutate(
    model_label = ifelse(
      model == "independent",
      "Nepriklausomas modelis",
      "Vasicek modelis"
    ),
    p_label = paste0("PD = ", p),
    n_label = paste0("n = ", n),
    rho_label = paste0("rho = ", rho)
  )

# Empirinė pasikliovimo intervalo padengimo tikimybė
coverage_table <- results %>%
  select(model, method, p, n, rho, coverage) %>%
  arrange(model, p, n, rho, method)

write.csv(
  coverage_table,
  file = "results/coverage_table.csv",
  row.names = FALSE
)

# Padengimo tikimybė nepriklausomam atvejui - pasirinkti variantai
coverage_independent_selected <- results %>%
  filter(
    model == "independent",
    p %in% c(0.001, 0.01, 0.05),
    n %in% c(5000, 50000, 100000)
  ) %>%
  select(p, n, method, coverage) %>%
  pivot_wider(names_from = method, values_from = coverage) %>%
  arrange(p, n)

write.csv(
  coverage_independent_selected,
  "results/coverage_independent_selected.csv",
  row.names = FALSE
)

# Pilnoji lentelė
coverage_independent_full <- results %>%
  filter(model == "independent") %>%
  select(p, n, method, coverage) %>%
  pivot_wider(names_from = method, values_from = coverage) %>%
  arrange(p, n)

write.csv(
  coverage_independent_full,
  "results/coverage_independent_full.csv",
  row.names = FALSE
)

# Vidutinis pasikliovimo intervalo plotis
width_table <- results %>%
  select(model, method, p, n, rho, avg_width, median_width) %>%
  arrange(model, p, n, rho, method)

write.csv(
  width_table,
  file = "results/width_table.csv",
  row.names = FALSE
)

# Supaprastinta vidutinio intervalo pločio lentelė nepriklausomam atvejui
width_independent_selected <- results %>%
  filter(
    model == "independent",
    p %in% c(0.001, 0.01, 0.05),
    n %in% c(5000, 50000, 100000)
  ) %>%
  select(p, n, method, avg_width) %>%
  pivot_wider(names_from = method, values_from = avg_width) %>%
  arrange(p, n)

write.csv(
  width_independent_selected,
  "results/width_independent_selected.csv",
  row.names = FALSE
)

# Pilna vidutinio intervalo pločio lentelė nepriklausomam atvejui
width_independent_full <- results %>%
  filter(model == "independent") %>%
  select(p, n, method, avg_width) %>%
  pivot_wider(names_from = method, values_from = avg_width) %>%
  arrange(p, n)

write.csv(
  width_independent_full,
  "results/width_independent_full.csv",
  row.names = FALSE
)

# Duomenys nepriklausomo modelio grafikams
ind_data <- results_plot %>%
  filter(model == "independent") %>%
  mutate(
    method = factor(
      method,
      levels = c("Wald", "Clopper_Pearson", "Wilson", "Agresti_Coull"),
      labels = c("Wald", "Clopper-Pearson", "Wilson", "Agresti-Coull")
    ),
    p_label = paste0("PD = ", format(p, decimal.mark = ",", scientific = FALSE))
  )


#---------Padengimo tikimybė nepriklausomam atvejui pagal portfolio dydį--------

make_coverage_plot <- function(data, pd_value) {
  
  data %>%
    filter(p == pd_value) %>%
    ggplot(aes(x = n, y = coverage, color = method, group = method)) +
    geom_hline(
      yintercept = 1 - alpha,
      linetype = "dashed",
      linewidth = 0.6,
      color = "gray35"
    ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.4)+
    scale_x_continuous(
      breaks = sort(unique(data$n)),
      labels = function(x) format(x, scientific = FALSE, big.mark = " ")
    ) +
    scale_y_continuous(
      limits = c(0.7, 1.00),
      breaks = seq(0.4, 1.0, by = 0.1)
    ) +
    labs(
      title = paste0("PD = ", format(pd_value, decimal.mark = ",", scientific = FALSE)),
      x = "Portfelio dydis n",
      y = "Empirinė padengimo tikimybė",
      color = "Metodas"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "gray85"),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "gray85"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

p1 <- make_coverage_plot(ind_data, 0.001)
p2 <- make_coverage_plot(ind_data, 0.005)
p3 <- make_coverage_plot(ind_data, 0.01)
p4 <- make_coverage_plot(ind_data, 0.02)
p5 <- make_coverage_plot(ind_data, 0.05)

plot_ind_coverage <- (
  (p1 | p2 | p3) /
    (plot_spacer() | p4 | p5 | plot_spacer())
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

plot_ind_coverage <- plot_ind_coverage +
  plot_annotation(
    title = "Pasikliovimo intervalų padengimo tikimybė nepriklausomų įvykių atveju",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  )

print(plot_ind_coverage) # žiūrim

#Išaugojam
ggsave(
  filename = "plots/coverage_independent.png",
  plot = plot_ind_coverage,
  width = 12,
  height = 7.5,
  dpi = 300
)

#-----Vidutinis intervalo plotis nepriklausomam atvejui pagal portfolio dydį----

make_width_plot <- function(data, pd_value) {
  
  data %>%
    filter(p == pd_value) %>%
    ggplot(aes(x = n, y = avg_width, color = method, group = method)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.4) +
    scale_x_continuous(
      breaks = sort(unique(data$n)),
      labels = function(x) format(x, scientific = FALSE, big.mark = " ")
    ) +
    scale_y_continuous(
      labels = function(x) format(x, scientific = FALSE, decimal.mark = ",")
    ) +
    labs(
      title = paste0("PD = ", format(pd_value, decimal.mark = ",", scientific = FALSE)),
      x = "Portfelio dydis n",
      y = "Vidutinis intervalo plotis",
      color = "Metodas"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "gray85"),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "gray85"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

w1 <- make_width_plot(ind_data, 0.001)
w2 <- make_width_plot(ind_data, 0.005)
w3 <- make_width_plot(ind_data, 0.01)
w4 <- make_width_plot(ind_data, 0.02)
w5 <- make_width_plot(ind_data, 0.05)

plot_ind_width <- (
  (w1 | w2 | w3) /
    (plot_spacer() | w4 | w5 | plot_spacer())
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

plot_ind_width <- plot_ind_width +
  plot_annotation(
    title = "Vidutinis pasikliovimo intervalų plotis nepriklausomų įvykių atveju",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  )

print(plot_ind_width)

ggsave(
  filename = "plots/width_independent.png",
  plot = plot_ind_width,
  width = 12,
  height = 7.5,
  dpi = 300
)
#--------------Duomenys standartinio nuokrypio analizei-------------------------

sd_phat_data <- results %>%
  select(model, p, n, rho, mean_p_hat, sd_p_hat, mean_defaults, sd_defaults) %>%
  distinct() %>%
  mutate(
    model_label = ifelse(
      model == "independent",
      "Nepriklausomas modelis",
      "Vasicek modelis"
    ),
    p_label = paste0("PD = ", format(p, decimal.mark = ",", scientific = FALSE)),
    n_label = paste0("n = ", format(n, scientific = FALSE, big.mark = " ")),
    rho_label = paste0(
      "rho = ",
      gsub("\\.", ",", sprintf("%.2f", rho))
    )
  )

# Supaprastinta standartinio nuokrypio lentelė pagrindiniam tekstui
sd_phat_selected <- sd_phat_data %>%
  filter(
    p %in% c(0.001, 0.01, 0.05),
    n %in% c(5000, 50000, 100000),
    (model == "independent" & rho == 0) |
      (model == "vasicek" & rho %in% c(0.01, 0.05, 0.10, 0.20, 0.25))
  ) %>%
  mutate(
    rho_column = case_when(
      model == "independent" ~ "rho_0",
      rho == 0.01 ~ "rho_001",
      rho == 0.05 ~ "rho_005",
      rho == 0.10 ~ "rho_010",
      rho == 0.20 ~ "rho_020",
      rho == 0.25 ~ "rho_025"
    )
  ) %>%
  select(p, n, rho_column, sd_p_hat) %>%
  pivot_wider(names_from = rho_column, values_from = sd_p_hat) %>%
  arrange(p, n)

write.csv(
  sd_phat_selected,
  "results/sd_phat_selected.csv",
  row.names = FALSE
)

# Pilna standartinio nuokrypio lentelė priedams
sd_phat_full <- sd_phat_data %>%
  filter(
    (model == "independent" & rho == 0) |
      (model == "vasicek")
  ) %>%
  mutate(
    rho_column = case_when(
      model == "independent" ~ "rho_0",
      rho == 0.01 ~ "rho_001",
      rho == 0.05 ~ "rho_005",
      rho == 0.10 ~ "rho_010",
      rho == 0.20 ~ "rho_020",
      rho == 0.25 ~ "rho_025"
    )
  ) %>%
  select(p, n, rho_column, sd_p_hat) %>%
  pivot_wider(names_from = rho_column, values_from = sd_p_hat) %>%
  arrange(p, n)

write.csv(
  sd_phat_full,
  "results/sd_phat_full.csv",
  row.names = FALSE
)

#---------Vašiček modelis: standartinis nuokrypis pagal koreliaciją-------------

sd_vasicek_data <- sd_phat_data %>%
  filter(model == "vasicek") %>%
  mutate(
    n_label = paste0("n = ", format(n, scientific = FALSE, big.mark = " ")),
    p_label = paste0("PD = ", format(p, decimal.mark = ",", scientific = FALSE))
  )


pd_values <- c(0.001, 0.005, 0.01, 0.02, 0.05)

y_limits <- sd_vasicek_data %>%
  filter(p %in% pd_values) %>%
  summarise(
    ymin = min(sd_p_hat, na.rm = TRUE),
    ymax = max(sd_p_hat, na.rm = TRUE)
  ) %>%
  unlist()


make_sd_phat_plot <- function(data, pd_value) {
  
  data %>%
    filter(p == pd_value) %>%
    ggplot(aes(x = rho, y = sd_p_hat, color = n_label, group = n_label)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.3) +
    scale_x_continuous(
      breaks = sort(unique(data$rho)),
      labels = function(x) gsub("\\.", ",", sprintf("%.2f", x))
    ) +
    scale_y_continuous(
      labels = function(x) format(
        x,
        scientific = FALSE,
        decimal.mark = ",",
        trim = TRUE
      )
    ) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = paste0("PD = ", format(pd_value, decimal.mark = ",", scientific = FALSE)),
      x = expression("Turto koreliacija " * rho),
      y = expression(SD(hat(p))),
      color = "Portfelio dydis"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "gray85"),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "gray85"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

s1 <- make_sd_phat_plot(sd_vasicek_data, 0.001)
s2 <- make_sd_phat_plot(sd_vasicek_data, 0.005)
s3 <- make_sd_phat_plot(sd_vasicek_data, 0.01)
s4 <- make_sd_phat_plot(sd_vasicek_data, 0.02)
s5 <- make_sd_phat_plot(sd_vasicek_data, 0.05)

plot_sd_phat <- (
  (s1 | s2 | s3) /
    (plot_spacer() | s4 | s5 | plot_spacer())
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

plot_sd_phat <- plot_sd_phat +
  plot_annotation(
    title = "Stebimo nemokumo dažnio standartinis nuokrypis Vašiček modelyje",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  )

print(plot_sd_phat)

ggsave(
  filename = "plots/sd_p_hat_vasicek.png",
  plot = plot_sd_phat,
  width = 12,
  height = 7.5,
  dpi = 300
)

#------------------Duomenys Vašiček modelio grafikams---------------------------

vasicek_data <- results_plot %>%
  filter(model == "vasicek") %>%
  mutate(
    method = factor(
      method,
      levels = c("Wald", "Clopper_Pearson", "Wilson", "Agresti_Coull"),
      labels = c("Wald", "Clopper-Pearson", "Wilson", "Agresti-Coull")
    ),
    p_label = paste0("PD = ", format(p, decimal.mark = ",", scientific = FALSE)),
    n_label = paste0("n = ", format(n, scientific = FALSE, big.mark = " ")),
    rho_label = paste0(
      "rho = ",
      gsub("\\.", ",", sprintf("%.2f", rho))
    )
  )

# Padengimo tikimybė - supaprastinta
coverage_vasicek_selected <- results %>%
  filter(
    model == "vasicek",
    p %in% c(0.001, 0.01, 0.05),
    n == 10000
  ) %>%
  select(p, n, rho, method, coverage) %>%
  pivot_wider(names_from = method, values_from = coverage) %>%
  arrange(p, rho)

write.csv(
  coverage_vasicek_selected,
  "results/coverage_vasicek_selected.csv",
  row.names = FALSE
)

# Padengimo tikimybė - pilnoji
coverage_vasicek_full <- results %>%
  filter(model == "vasicek") %>%
  select(p, n, rho, method, coverage) %>%
  pivot_wider(names_from = method, values_from = coverage) %>%
  arrange(p, n, rho)

write.csv(
  coverage_vasicek_full,
  "results/coverage_vasicek_full.csv",
  row.names = FALSE
)

#  Vidutinio pločio supaprastinta lentelė
width_vasicek_selected <- results %>%
  filter(
    model == "vasicek",
    p %in% c(0.001, 0.01, 0.05),
    n == 10000
  ) %>%
  select(p, n, rho, method, avg_width) %>%
  pivot_wider(names_from = method, values_from = avg_width) %>%
  arrange(p, rho)

write.csv(
  width_vasicek_selected,
  "results/width_vasicek_selected.csv",
  row.names = FALSE
)

#  Vidutinio pločio pilnoji lentelė
width_vasicek_full <- results %>%
  filter(model == "vasicek") %>%
  select(p, n, rho, method, avg_width) %>%
  pivot_wider(names_from = method, values_from = avg_width) %>%
  arrange(p, n, rho)

write.csv(
  width_vasicek_full,
  "results/width_vasicek_full.csv",
  row.names = FALSE
)

#------Vašiček modelis: padengimo tikimybė pagal koreliaciją pasirinktam n------

make_vasicek_coverage_plot <- function(data, pd_value, n_value = 10000) {
  
  data %>%
    filter(p == pd_value, n == n_value) %>%
    ggplot(aes(x = rho, y = coverage, color = method, group = method)) +
    geom_hline(
      yintercept = 1 - alpha,
      linetype = "dashed",
      linewidth = 0.6,
      color = "gray35"
    ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.4) +
    scale_x_continuous(
      breaks = sort(unique(data$rho)),
      labels = function(x) gsub("\\.", ",", sprintf("%.2f", x))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1),
      labels = function(x) format(x, scientific = FALSE, decimal.mark = ",")
    ) +
    labs(
      title = paste0("PD = ", format(pd_value, decimal.mark = ",", scientific = FALSE)),
      x = expression("Turto koreliacija " * rho),
      y = "Empirinė padengimo tikimybė",
      color = "Metodas"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "gray85"),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "gray85"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

vc1 <- make_vasicek_coverage_plot(vasicek_data, 0.001, 10000)
vc2 <- make_vasicek_coverage_plot(vasicek_data, 0.005, 10000)
vc3 <- make_vasicek_coverage_plot(vasicek_data, 0.01, 10000)
vc4 <- make_vasicek_coverage_plot(vasicek_data, 0.02, 10000)
vc5 <- make_vasicek_coverage_plot(vasicek_data, 0.05, 10000)

plot_vasicek_coverage <- (
  (vc1 | vc2 | vc3) /
    (plot_spacer() | vc4 | vc5 | plot_spacer())
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

plot_vasicek_coverage <- plot_vasicek_coverage +
  plot_annotation(
    title = "Pasikliovimo intervalų padengimo tikimybė Vašiček modelyje, kai n = 10 000",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  )

print(plot_vasicek_coverage)

ggsave(
  filename = "plots/coverage_vasicek_n10000.png",
  plot = plot_vasicek_coverage,
  width = 12,
  height = 7.5,
  dpi = 300
)

#--------Vašiček modelis: vid plotis pagal koreliaciją pasirinktam n------------

make_vasicek_width_plot <- function(data, pd_value, n_value = 10000) {
  
  data %>%
    filter(p == pd_value, n == n_value) %>%
    ggplot(aes(x = rho, y = avg_width, color = method, group = method)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.4) +
    scale_x_continuous(
      breaks = sort(unique(data$rho)),
      labels = function(x) gsub("\\.", ",", sprintf("%.2f", x))
    ) +
    scale_y_continuous(
      labels = function(x) format(
        x,
        scientific = FALSE,
        decimal.mark = ",",
        trim = TRUE
      )
    ) +
    labs(
      title = paste0("PD = ", format(pd_value, decimal.mark = ",", scientific = FALSE)),
      x = expression("Turto koreliacija " * rho),
      y = "Vidutinis intervalo plotis",
      color = "Metodas"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "gray85"),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "gray85"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

vw1 <- make_vasicek_width_plot(vasicek_data, 0.001, 10000)
vw2 <- make_vasicek_width_plot(vasicek_data, 0.005, 10000)
vw3 <- make_vasicek_width_plot(vasicek_data, 0.01, 10000)
vw4 <- make_vasicek_width_plot(vasicek_data, 0.02, 10000)
vw5 <- make_vasicek_width_plot(vasicek_data, 0.05, 10000)

plot_vasicek_width <- (
  (vw1 | vw2 | vw3) /
    (plot_spacer() | vw4 | vw5 | plot_spacer())
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

plot_vasicek_width <- plot_vasicek_width +
  plot_annotation(
    title = "Vidutinis pasikliovimo intervalų plotis Vašiček modelyje, kai n = 10 000",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  )

print(plot_vasicek_width)

ggsave(
  filename = "plots/width_vasicek_n10000.png",
  plot = plot_vasicek_width,
  width = 12,
  height = 7.5,
  dpi = 300
)


#-------Nepriklausomo ir Vašiček atvejų padengimo tikimybės palyginimas---------

coverage_ind <- results %>%
  filter(model == "independent", rho == 0) %>%
  select(method, p, n, coverage_ind = coverage)

coverage_vas <- results %>%
  filter(model == "vasicek") %>%
  select(method, p, n, rho, coverage_vasicek = coverage)

coverage_comparison <- coverage_vas %>%
  left_join(
    coverage_ind,
    by = c("method", "p", "n")
  ) %>%
  mutate(
    coverage_drop = coverage_ind - coverage_vasicek
  ) %>%
  arrange(p, n, rho, method)

write.csv(
  coverage_comparison,
  "results/coverage_comparison_ind_vs_vasicek.csv",
  row.names = FALSE
)

# Supaprastinta lentelė
coverage_drop_selected <- coverage_comparison %>%
  filter(
    n == 10000,
    p %in% c(0.001, 0.01, 0.05),
    rho %in% c(0.01, 0.10, 0.25)
  ) %>%
  select(p, rho, method, coverage_drop) %>%
  pivot_wider(names_from = method, values_from = coverage_drop) %>%
  arrange(p, rho)

write.csv(
  coverage_drop_selected,
  "results/coverage_drop_selected.csv",
  row.names = FALSE
)

# Pilnoji lentelė
coverage_drop_full <- coverage_comparison %>%
  select(p, n, rho, method, coverage_ind, coverage_vasicek, coverage_drop) %>%
  arrange(p, n, rho, method)

write.csv(
  coverage_drop_full,
  "results/coverage_drop_full.csv",
  row.names = FALSE
)

#-------Nepriklausomo ir Vašiček atvejų padengimo tikimybės palyginimo grafikas---------

comparison_plot_data <- coverage_comparison %>%
  filter(n == 10000) %>%
  mutate(
    method = factor(
      method,
      levels = c("Wald", "Clopper_Pearson", "Wilson", "Agresti_Coull"),
      labels = c("Wald", "Clopper-Pearson", "Wilson", "Agresti-Coull")
    )
  )

make_coverage_drop_plot <- function(data, pd_value) {
  
  data %>%
    filter(p == pd_value) %>%
    ggplot(aes(x = rho, y = coverage_drop, color = method, group = method)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.5,
      color = "gray35"
    ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.4) +
    scale_x_continuous(
      breaks = sort(unique(data$rho)),
      labels = function(x) gsub("\\.", ",", sprintf("%.2f", x))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1),
      labels = function(x) format(x, scientific = FALSE, decimal.mark = ",")
    ) +
    labs(
      title = paste0("PD = ", format(pd_value, decimal.mark = ",", scientific = FALSE)),
      x = expression("Turto koreliacija " * rho),
      y = expression(Delta * " Coverage"),
      color = "Metodas"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "gray85"),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "gray85"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

cd1 <- make_coverage_drop_plot(comparison_plot_data, 0.001)
cd2 <- make_coverage_drop_plot(comparison_plot_data, 0.005)
cd3 <- make_coverage_drop_plot(comparison_plot_data, 0.01)
cd4 <- make_coverage_drop_plot(comparison_plot_data, 0.02)
cd5 <- make_coverage_drop_plot(comparison_plot_data, 0.05)

plot_coverage_drop <- (
  (cd1 | cd2 | cd3) /
    (plot_spacer() | cd4 | cd5 | plot_spacer())
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

plot_coverage_drop <- plot_coverage_drop +
  plot_annotation(
    title = "Padengimo tikimybės sumažėjimas pereinant nuo nepriklausomo prie Vašiček modelio, kai n = 10 000",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  )

print(plot_coverage_drop)

ggsave(
  filename = "plots/coverage_drop_n10000.png",
  plot = plot_coverage_drop,
  width = 12,
  height = 7.5,
  dpi = 300
)

#------Apatinės pasikliovimo intervalo ribos, lygios nuliui, dažnio lentelės----

lower_zero_full <- results %>%
  select(model, method, p, n, rho, lower_zero_rate) %>%
  arrange(model, p, n, rho, method)

write.csv(
  lower_zero_full,
  "results/lower_zero_full.csv",
  row.names = FALSE
)

# Supaprastinta
lower_zero_selected <- results %>%
  filter(
    p %in% c(0.001, 0.005, 0.01),
    n == 10000,
    (model == "independent" & rho == 0) |
      (model == "vasicek" & rho %in% c(0.01, 0.10, 0.25))
  ) %>%
  select(model, p, n, rho, method, lower_zero_rate) %>%
  arrange(model, p, rho, method)

write.csv(
  lower_zero_selected,
  "results/lower_zero_selected.csv",
  row.names = FALSE
)

#------Apatinės pasikliovimo intervalo ribos, lygios nuliui, grafikas-----------
lower_zero_plot_data <- results_plot %>%
  filter(
    p %in% c(0.001, 0.005, 0.01),
    n == 10000,
    (model == "independent" & rho == 0) |
      (model == "vasicek")
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c("Wald", "Clopper_Pearson", "Wilson", "Agresti_Coull"),
      labels = c("Wald", "Clopper-Pearson", "Wilson", "Agresti-Coull")
    ),
    p_label = paste0("PD = ", format(p, decimal.mark = ",", scientific = FALSE)),
    model_label = ifelse(model == "independent", "Nepriklausomas", "Vasicek")
  )

plot_lower_zero_selected <- lower_zero_plot_data %>%
  filter(model == "vasicek") %>%
  ggplot(aes(x = rho, y = lower_zero_rate, color = method, group = method)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.4) +
  facet_wrap(~ p_label) +
  scale_x_continuous(
    breaks = sort(unique(lower_zero_plot_data$rho)),
    labels = function(x) gsub("\\.", ",", sprintf("%.2f", x))
  ) +
  labs(
    title = "Apatinės intervalo ribos, lygios nuliui, dažnis Vašiček modelyje, kai n = 10 000",
    x = expression("Turto koreliacija " * rho),
    y = "Dažnis",
    color = "Metodas"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

print(plot_lower_zero_selected)

ggsave(
  filename = "plots/lower_zero_rate_vasicek_n10000.png",
  plot = plot_lower_zero_selected,
  width = 10,
  height = 5.5,
  dpi = 300
)


#-------Nepriklausomo ir Vašiček atvejų vid intervalo pločio palyginimas--------

width_ind <- results %>%
  filter(model == "independent", rho == 0) %>%
  select(method, p, n, width_ind = avg_width)

width_vas <- results %>%
  filter(model == "vasicek") %>%
  select(method, p, n, rho, width_vasicek = avg_width)

width_comparison <- width_vas %>%
  left_join(
    width_ind,
    by = c("method", "p", "n")
  ) %>%
  mutate(
    width_diff = width_vasicek - width_ind
  ) %>%
  arrange(p, n, rho, method)

write.csv(
  width_comparison,
  "results/width_comparison_ind_vs_vasicek.csv",
  row.names = FALSE
)

# Supaprastinta
width_diff_selected <- width_comparison %>%
  filter(
    n == 10000,
    p %in% c(0.001, 0.01, 0.05),
    rho %in% c(0.01, 0.10, 0.25)
  ) %>%
  select(p, rho, method, width_diff) %>%
  pivot_wider(names_from = method, values_from = width_diff) %>%
  arrange(p, rho)

write.csv(
  width_diff_selected,
  "results/width_diff_selected.csv",
  row.names = FALSE
)

#-------Nepriklausomo ir Vašiček atvejų vid intervalo pločio grafikas-----------

comparison_width_plot_data <- width_comparison %>%
  filter(n == 10000) %>%
  mutate(
    method = factor(
      method,
      levels = c("Wald", "Clopper_Pearson", "Wilson", "Agresti_Coull"),
      labels = c("Wald", "Clopper-Pearson", "Wilson", "Agresti-Coull")
    )
  )

make_width_diff_plot <- function(data, pd_value) {
  
  data %>%
    filter(p == pd_value) %>%
    ggplot(aes(x = rho, y = width_diff, color = method, group = method)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.5,
      color = "gray35"
    ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.4) +
    expand_limits(y = 0) +
    scale_x_continuous(
      breaks = sort(unique(data$rho)),
      labels = function(x) gsub("\\.", ",", sprintf("%.2f", x))
    ) +
    scale_y_continuous(
      labels = function(x) format(
        x,
        scientific = FALSE,
        decimal.mark = ",",
        trim = TRUE
      )
    ) +
    labs(
      title = paste0("PD = ", format(pd_value, decimal.mark = ",", scientific = FALSE)),
      x = expression("Turto koreliacija " * rho),
      y = expression(Delta * " Width"),
      color = "Metodas"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "gray85"),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "gray85"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

wd1 <- make_width_diff_plot(comparison_width_plot_data, 0.001)
wd2 <- make_width_diff_plot(comparison_width_plot_data, 0.005)
wd3 <- make_width_diff_plot(comparison_width_plot_data, 0.01)
wd4 <- make_width_diff_plot(comparison_width_plot_data, 0.02)
wd5 <- make_width_diff_plot(comparison_width_plot_data, 0.05)

plot_width_diff <- (
  (wd1 | wd2 | wd3) /
    (plot_spacer() | wd4 | wd5 | plot_spacer())
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

plot_width_diff <- plot_width_diff +
  plot_annotation(
    title = "Vidutinio intervalo pločio pokytis pereinant nuo nepriklausomo prie Vašiček modelio, kai n = 10 000",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  )

print(plot_width_diff)

ggsave(
  filename = "plots/width_diff_n10000.png",
  plot = plot_width_diff,
  width = 12,
  height = 7.5,
  dpi = 300
)

#---------Efektyvus imties dydis: empirinė defaultų koreliacija-----------------

#Vašiček simuliacija su defaultų skaičiumi ir dviem porų stulpeliais
simulate_vasicek_counts_and_pairs <- function(B, n, p, rho,
                                              pairs_per_sim = 1000,
                                              chunk_size = 100) {
  
  threshold <- qnorm(p)
  counts <- numeric(B)
  
  total_pairs <- B * pairs_per_sim
  D1 <- integer(total_pairs)
  D2 <- integer(total_pairs)
  
  pos <- 1
  start_indices <- seq(1, B, by = chunk_size)
  
  for (start in start_indices) {
    
    end <- min(start + chunk_size - 1, B)
    current_B <- end - start + 1
    
    # Bendras sisteminis rizikos veiksnys kiekvienam portfeliui
    y <- rnorm(current_B)
    
    # Individualūs veiksniai kiekvienam skolininkui
    z <- matrix(
      rnorm(current_B * n),
      nrow = current_B,
      ncol = n
    )
    
    # Latentinis kintamasis pagal Vašiček modelį
    asset_values <- sqrt(1 - rho) * z
    
    asset_values <- sweep(
      asset_values,
      MARGIN = 1,
      STATS = sqrt(rho) * y,
      FUN = "+"
    )
    
    # Įsipareigojimų nevykdymo indikatoriai
    defaults <- asset_values <= threshold
    
    # Defaultų skaičius kiekviename portfelyje
    counts[start:end] <- rowSums(defaults)
    
    # Atsitiktinės skolininkų poros kiekvienoje simuliacijoje
    for (b in 1:current_B) {
      
      i <- sample.int(n, pairs_per_sim, replace = TRUE)
      
      # Parenkame j taip, kad j nelygu i
      j <- sample.int(n - 1, pairs_per_sim, replace = TRUE)
      j <- ifelse(j >= i, j + 1, j)
      
      idx <- pos:(pos + pairs_per_sim - 1)
      
      D1[idx] <- as.integer(defaults[b, i])
      D2[idx] <- as.integer(defaults[b, j])
      
      pos <- pos + pairs_per_sim
    }
  }
  
  list(
    counts = counts,
    D1 = D1,
    D2 = D2
  )
}

#---------Empirinė defaultų koreliacija iš dviejų stulpelių---------------------

estimate_default_corr_two_columns <- function(D1, D2) {
  
  rho_default <- suppressWarnings(cor(D1, D2))
  
  # Apsauga, jei dėl labai retų įvykių koreliacija negali būti apskaičiuota
  if (is.na(rho_default)) {
    rho_default <- 0
  }
  
  # Efektyvaus imties dydžio korekcijai neigiamos reikšmės nėra naudojamos
  rho_default <- max(rho_default, 0)
  
  rho_default
}

#--------------Efektyvus imties dydis-------------------------------------------

effective_sample_size_from_corr <- function(n, rho_default) {
  
  deff <- 1 + (n - 1) * rho_default
  n_eff <- n / deff
  
  list(
    rho_default = rho_default,
    deff = deff,
    n_eff = n_eff
  )
}

#-----------------Koreguoti pasikliovimo intervalai su n_eff--------------------

corrected_intervals_eff_n <- function(x, n, n_eff, alpha = 0.05) {
  
  z <- qnorm(1 - alpha / 2)
  p_hat <- x / n
  
  # Kad formulėse būtų galima naudoti pseudo įvykių skaičių
  x_eff <- p_hat * n_eff
  

  # Wald
  se_eff <- sqrt(p_hat * (1 - p_hat) / n_eff)
  
  wald_l <- max(0, p_hat - z * se_eff)
  wald_u <- min(1, p_hat + z * se_eff)
  

  # Wilson
  denom <- 1 + z^2 / n_eff
  
  wilson_center <- (p_hat + z^2 / (2 * n_eff)) / denom
  
  wilson_half <- (
    z * sqrt(
      p_hat * (1 - p_hat) / n_eff +
        z^2 / (4 * n_eff^2)
    )
  ) / denom
  
  wilson_l <- max(0, wilson_center - wilson_half)
  wilson_u <- min(1, wilson_center + wilson_half)
  

  # Agresti-Coull
  n_tilde <- n_eff + z^2
  p_tilde <- (x_eff + z^2 / 2) / n_tilde
  
  ac_half <- z * sqrt(p_tilde * (1 - p_tilde) / n_tilde)
  
  ac_l <- max(0, p_tilde - ac_half)
  ac_u <- min(1, p_tilde + ac_half)
  

  # Clopper-Pearson pseudo korekcija
  cp_l <- ifelse(
    x_eff <= 0,
    0,
    qbeta(alpha / 2, x_eff, n_eff - x_eff + 1)
  )
  
  cp_u <- ifelse(
    x_eff >= n_eff,
    1,
    qbeta(1 - alpha / 2, x_eff + 1, n_eff - x_eff)
  )
  
  tibble(
    method = c(
      "Wald",
      "Clopper_Pearson",
      "Wilson",
      "Agresti_Coull"
    ),
    lower = c(wald_l, cp_l, wilson_l, ac_l),
    upper = c(wald_u, cp_u, wilson_u, ac_u),
    width = upper - lower
  )
}

# Įvertinimas
evaluate_corrected_scenario_eff_n <- function(counts, n, p, rho,rho_default, deff, n_eff, alpha = 0.05) {
  
  map_dfr(
    counts,
    ~ corrected_intervals_eff_n(
      x = .x,
      n = n,
      n_eff = n_eff,
      alpha = alpha
    ),
    .id = "simulation"
  ) %>%
    group_by(method) %>%
    summarise(
      coverage = mean(lower <= p & upper >= p),
      avg_width = mean(width),
      median_width = median(width),
      lower_zero_rate = mean(lower == 0),
      .groups = "drop"
    ) %>%
    mutate(
      p = p,
      n = n,
      rho = rho,
      rho_default = rho_default,
      deff = deff,
      n_eff = n_eff,
      model = "vasicek_eff_n"
    )
}

#--------------------Korekcijos analizė, kai n = 10000--------------------------

n_correction <- 10000 
pairs_per_sim <- 1000

corrected_results_eff_n <- expand.grid(
  p = pd_values,
  rho = rho_values
) %>%
  as_tibble() %>%
  pmap_dfr(function(p, rho) {
    
    cat("Efektyvus n: p =", p, "n =", n_correction, "rho =", rho, "\n")
    
    sim <- simulate_vasicek_counts_and_pairs(
      B = B,
      n = n_correction,
      p = p,
      rho = rho,
      pairs_per_sim = pairs_per_sim,
      chunk_size = 100
    )
    
    rho_default <- estimate_default_corr_two_columns(
      D1 = sim$D1,
      D2 = sim$D2
    )
    
    eff <- effective_sample_size_from_corr(
      n = n_correction,
      rho_default = rho_default
    )
    
    evaluate_corrected_scenario_eff_n(
      counts = sim$counts,
      n = n_correction,
      p = p,
      rho = rho,
      rho_default = eff$rho_default,
      deff = eff$deff,
      n_eff = eff$n_eff,
      alpha = alpha
    )
  })

write.csv(
  corrected_results_eff_n,
  "results/corrected_results_eff_n_n10000.csv",
  row.names = FALSE
)


# Efektyvaus imties dydžio lentelė
eff_n_table <- corrected_results_eff_n %>%
  select(p, n, rho, rho_default, deff, n_eff) %>%
  distinct() %>%
  arrange(p, rho)

write.csv(
  eff_n_table,
  "results/effective_sample_size_table.csv",
  row.names = FALSE
)
# Supaprastinta
eff_n_selected <- eff_n_table %>%
  filter(
    p %in% c(0.001, 0.01, 0.05),
    rho %in% c(0.01, 0.10, 0.25)
  )

write.csv(
  eff_n_selected,
  "results/effective_sample_size_selected.csv",
  row.names = FALSE
)

#-----------Efektyvus imties dydis pagal koreliaciją grafikas-------------------

plot_eff_n <- eff_n_table %>%
  mutate(
    p_label = paste0("PD = ", format(p, decimal.mark = ",", scientific = FALSE))
  ) %>%
  ggplot(aes(x = rho, y = n_eff, color = p_label, group = p_label)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  geom_hline(
    yintercept = n_correction,
    linetype = "dashed",
    linewidth = 0.5,
    color = "gray35"
  ) +
  scale_x_continuous(
    breaks = sort(unique(eff_n_table$rho)),
    labels = function(x) gsub("\\.", ",", sprintf("%.2f", x))
  ) +
  scale_y_continuous(
    labels = function(x) format(x, scientific = FALSE, big.mark = " ")
  ) +
  labs(
    title = "Efektyvus imties dydis taikant empirinę defaultų koreliaciją",
    x = expression("Turto koreliacija " * rho),
    y = expression(n[eff]),
    color = "Scenarijus"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

print(plot_eff_n)

ggsave(
  filename = "plots/effective_sample_size_n10000.png",
  plot = plot_eff_n,
  width = 10,
  height = 6,
  dpi = 300
)

#--------------------Klasikinių ir koreguotų intervalų palyginimas--------------

classic_vasicek_n10000 <- results %>%
  filter(
    model == "vasicek",
    n == n_correction
  ) %>%
  select(
    p, n, rho, method,
    coverage_classic = coverage,
    width_classic = avg_width
  )

corrected_n10000 <- corrected_results_eff_n %>%
  select(
    p, n, rho, method,
    coverage_corrected = coverage,
    width_corrected = avg_width,
    rho_default,
    deff,
    n_eff
  )

correction_comparison <- classic_vasicek_n10000 %>%
  left_join(
    corrected_n10000,
    by = c("p", "n", "rho", "method")
  ) %>%
  mutate(
    coverage_change = coverage_corrected - coverage_classic,
    width_change = width_corrected - width_classic,
    width_ratio = width_corrected / width_classic
  ) %>%
  arrange(p, rho, method)

write.csv(
  correction_comparison,
  "results/correction_comparison_n10000.csv",
  row.names = FALSE
)

# Supaprastinta
correction_comparison_selected <- correction_comparison %>%
  filter(
    p %in% c(0.001, 0.01, 0.05),
    rho %in% c(0.01, 0.10, 0.25)
  )

write.csv(
  correction_comparison_selected,
  "results/correction_comparison_selected.csv",
  row.names = FALSE
)

# Supaprastinta
coverage_correction_selected <- correction_comparison %>%
  filter(
    p %in% c(0.001, 0.01, 0.05),
    rho %in% c(0.01, 0.10, 0.25)
  ) %>%
  select(p, rho, method, coverage_classic, coverage_corrected) %>%
  mutate(
    method = recode(
      method,
      "Clopper_Pearson" = "Clopper--Pearson",
      "Agresti_Coull" = "Agresti--Coull"
    )
  ) %>%
  arrange(p, rho, method)

write.csv(
  coverage_correction_selected,
  "results/coverage_correction_selected.csv",
  row.names = FALSE
)

#-------------Padengimo tikymybė prieš ir po korekcijos grafikas----------------

# Duomenys grafikui
coverage_correction_plot_data <- correction_comparison %>%
  mutate(
    method = factor(
      method,
      levels = c("Wald", "Clopper_Pearson", "Wilson", "Agresti_Coull"),
      labels = c("Wald", "Clopper-Pearson", "Wilson", "Agresti-Coull")
    )
  ) %>%
  select(p, rho, method, coverage_classic, coverage_corrected) %>%
  pivot_longer(
    cols = c(coverage_classic, coverage_corrected),
    names_to = "interval_type",
    values_to = "coverage"
  ) %>%
  mutate(
    interval_type = recode(
      interval_type,
      "coverage_classic" = "Klasikinis",
      "coverage_corrected" = "Koreguotas"
    ),
    p_label = paste0("PD = ", format(p, decimal.mark = ",", scientific = FALSE))
  )

make_correction_coverage_plot <- function(data, pd_value) {
  
  data %>%
    filter(p == pd_value) %>%
    ggplot(
      aes(
        x = rho,
        y = coverage,
        color = method,
        linetype = interval_type,
        group = interaction(method, interval_type)
      )
    ) +
    geom_hline(
      yintercept = 1 - alpha,
      linetype = "dashed",
      linewidth = 0.5,
      color = "gray35"
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.3) +
    scale_x_continuous(
      breaks = sort(unique(data$rho)),
      labels = function(x) gsub("\\.", ",", sprintf("%.2f", x))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1),
      labels = function(x) format(x, scientific = FALSE, decimal.mark = ",")
    ) +
    labs(
      title = paste0("PD = ", format(pd_value, decimal.mark = ",", scientific = FALSE)),
      x = expression("Turto koreliacija " * rho),
      y = "Empirinė padengimo tikimybė",
      color = "Metodas",
      linetype = "Intervalas"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "gray85"),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "gray85"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

cc1 <- make_correction_coverage_plot(coverage_correction_plot_data, 0.001)
cc2 <- make_correction_coverage_plot(coverage_correction_plot_data, 0.005)
cc3 <- make_correction_coverage_plot(coverage_correction_plot_data, 0.01)
cc4 <- make_correction_coverage_plot(coverage_correction_plot_data, 0.02)
cc5 <- make_correction_coverage_plot(coverage_correction_plot_data, 0.05)

plot_correction_coverage <- (
  (cc1 | cc2 | cc3) /
    (plot_spacer() | cc4 | cc5 | plot_spacer())
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

plot_correction_coverage <- plot_correction_coverage +
  plot_annotation(
    title = "Klasikinių ir efektyviu imties dydžiu koreguotų intervalų padengimo tikimybė, kai n = 10 000",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  )

print(plot_correction_coverage)

ggsave(
  filename = "plots/correction_coverage_n10000.png",
  plot = plot_correction_coverage,
  width = 12,
  height = 7.5,
  dpi = 300
)


#------------Vid intervalo plotis prieš ir po korekcijos grafikas---------------

# Duomenys grafikui
width_correction_plot_data <- correction_comparison %>%
  mutate(
    method = factor(
      method,
      levels = c("Wald", "Clopper_Pearson", "Wilson", "Agresti_Coull"),
      labels = c("Wald", "Clopper-Pearson", "Wilson", "Agresti-Coull")
    )
  ) %>%
  select(p, rho, method, width_classic, width_corrected) %>%
  pivot_longer(
    cols = c(width_classic, width_corrected),
    names_to = "interval_type",
    values_to = "width"
  ) %>%
  mutate(
    interval_type = recode(
      interval_type,
      "width_classic" = "Klasikinis",
      "width_corrected" = "Koreguotas"
    ),
    p_label = paste0("PD = ", format(p, decimal.mark = ",", scientific = FALSE))
  )

make_correction_width_plot <- function(data, pd_value) {
  
  data %>%
    filter(p == pd_value) %>%
    ggplot(
      aes(
        x = rho,
        y = width,
        color = method,
        linetype = interval_type,
        group = interaction(method, interval_type)
      )
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.3) +
    scale_x_continuous(
      breaks = sort(unique(data$rho)),
      labels = function(x) gsub("\\.", ",", sprintf("%.2f", x))
    ) +
    scale_y_continuous(
      labels = function(x) format(
        x,
        scientific = FALSE,
        decimal.mark = ",",
        trim = TRUE
      )
    ) +
    labs(
      title = paste0("PD = ", format(pd_value, decimal.mark = ",", scientific = FALSE)),
      x = expression("Turto koreliacija " * rho),
      y = "Vidutinis intervalo plotis",
      color = "Metodas",
      linetype = "Intervalas"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "gray85"),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "gray85"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

cw1 <- make_correction_width_plot(width_correction_plot_data, 0.001)
cw2 <- make_correction_width_plot(width_correction_plot_data, 0.005)
cw3 <- make_correction_width_plot(width_correction_plot_data, 0.01)
cw4 <- make_correction_width_plot(width_correction_plot_data, 0.02)
cw5 <- make_correction_width_plot(width_correction_plot_data, 0.05)

plot_correction_width <- (
  (cw1 | cw2 | cw3) /
    (plot_spacer() | cw4 | cw5 | plot_spacer())
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

plot_correction_width <- plot_correction_width +
  plot_annotation(
    title = "Klasikinių ir efektyviu imties dydžiu koreguotų intervalų plotis, kai n = 10 000",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  )

print(plot_correction_width)

ggsave(
  filename = "plots/correction_width_n10000.png",
  plot = plot_correction_width,
  width = 12,
  height = 7.5,
  dpi = 300
)

#------------------------Korekcijos analizė, kai n = 100000---------------------

n_correction <- 100000
pairs_per_sim <- 1000

corrected_results_eff_n <- expand.grid(
  p = pd_values,
  rho = rho_values
) %>%
  as_tibble() %>%
  pmap_dfr(function(p, rho) {
    
    cat("Efektyvus n: p =", p, "n =", n_correction, "rho =", rho, "\n")
    
    sim <- simulate_vasicek_counts_and_pairs(
      B = B,
      n = n_correction,
      p = p,
      rho = rho,
      pairs_per_sim = pairs_per_sim,
      chunk_size = 100
    )
    
    rho_default <- estimate_default_corr_two_columns(
      D1 = sim$D1,
      D2 = sim$D2
    )
    
    eff <- effective_sample_size_from_corr(
      n = n_correction,
      rho_default = rho_default
    )
    
    evaluate_corrected_scenario_eff_n(
      counts = sim$counts,
      n = n_correction,
      p = p,
      rho = rho,
      rho_default = eff$rho_default,
      deff = eff$deff,
      n_eff = eff$n_eff,
      alpha = alpha
    )
  })

write.csv(
  corrected_results_eff_n,
  "results/corrected_results_eff_n_n100000.csv",
  row.names = FALSE
)

print(corrected_results_eff_n)
corrected_results_eff_n <- read.csv("results/corrected_results_eff_n_n100000.csv")

#---------------------Efektyvaus imties dydžio lentelė--------------------------

eff_n_table <- corrected_results_eff_n %>%
  select(p, n, rho, rho_default, deff, n_eff) %>%
  distinct() %>%
  arrange(p, rho)

write.csv(
  eff_n_table,
  "results/effective_sample_size_table100000.csv",
  row.names = FALSE
)

# Supaprastinta
eff_n_selected <- eff_n_table %>%
  filter(
    p %in% c(0.001, 0.01, 0.05),
    rho %in% c(0.01, 0.10, 0.25)
  )

write.csv(
  eff_n_selected,
  "results/effective_sample_size_selected100000.csv",
  row.names = FALSE
)

#-----------Efektyvus imties dydis pagal koreliaciją grafikas-------------------

plot_eff_n <- eff_n_table %>%
  mutate(
    p_label = paste0("PD = ", format(p, decimal.mark = ",", scientific = FALSE))
  ) %>%
  ggplot(aes(x = rho, y = n_eff, color = p_label, group = p_label)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  geom_hline(
    yintercept = n_correction,
    linetype = "dashed",
    linewidth = 0.5,
    color = "gray35"
  ) +
  scale_x_continuous(
    breaks = sort(unique(eff_n_table$rho)),
    labels = function(x) gsub("\\.", ",", sprintf("%.2f", x))
  ) +
  scale_y_continuous(
    labels = function(x) format(x, scientific = FALSE, big.mark = " ")
  ) +
  labs(
    title = "Efektyvus imties dydis taikant empirinę defaultų koreliaciją",
    x = expression("Turto koreliacija " * rho),
    y = expression(n[eff]),
    color = "Scenarijus"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

print(plot_eff_n)

ggsave(
  filename = "plots/effective_sample_size_n100000.png",
  plot = plot_eff_n,
  width = 10,
  height = 6,
  dpi = 300
)

#------------Klasikinių ir koreguotų intervalų palyginimas----------------------

classic_vasicek_n100000 <- results %>%
  filter(
    model == "vasicek",
    n == n_correction
  ) %>%
  select(
    p, n, rho, method,
    coverage_classic = coverage,
    width_classic = avg_width
  )

corrected_n100000 <- corrected_results_eff_n %>%
  select(
    p, n, rho, method,
    coverage_corrected = coverage,
    width_corrected = avg_width,
    rho_default,
    deff,
    n_eff
  )

correction_comparison <- classic_vasicek_n100000 %>%
  left_join(
    corrected_n100000,
    by = c("p", "n", "rho", "method")
  ) %>%
  mutate(
    coverage_change = coverage_corrected - coverage_classic,
    width_change = width_corrected - width_classic,
    width_ratio = width_corrected / width_classic
  ) %>%
  arrange(p, rho, method)

write.csv(
  correction_comparison,
  "results/correction_comparison_n100000.csv",
  row.names = FALSE
)

correction_comparison_selected <- correction_comparison %>%
  filter(
    p %in% c(0.001, 0.01, 0.05),
    rho %in% c(0.01, 0.10, 0.25)
  )

write.csv(
  correction_comparison_selected,
  "results/correction_comparison_selected.csv",
  row.names = FALSE
)

coverage_correction_selected <- correction_comparison %>%
  filter(
    p %in% c(0.001, 0.01, 0.05),
    rho %in% c(0.01, 0.10, 0.25)
  ) %>%
  select(p, rho, method, coverage_classic, coverage_corrected) %>%
  mutate(
    method = recode(
      method,
      "Clopper_Pearson" = "Clopper-Pearson",
      "Agresti_Coull" = "Agresti-Coull"
    )
  ) %>%
  arrange(p, rho, method)

write.csv(
  coverage_correction_selected,
  "results/coverage_correction_selected100000.csv",
  row.names = FALSE
)

#-------------Padengimo tikymybė prieš ir po korekcijos grafikas----------------

# Duomenys grafikui
coverage_correction_plot_data <- correction_comparison %>%
  mutate(
    method = factor(
      method,
      levels = c("Wald", "Clopper_Pearson", "Wilson", "Agresti_Coull"),
      labels = c("Wald", "Clopper-Pearson", "Wilson", "Agresti-Coull")
    )
  ) %>%
  select(p, rho, method, coverage_classic, coverage_corrected) %>%
  pivot_longer(
    cols = c(coverage_classic, coverage_corrected),
    names_to = "interval_type",
    values_to = "coverage"
  ) %>%
  mutate(
    interval_type = recode(
      interval_type,
      "coverage_classic" = "Klasikinis",
      "coverage_corrected" = "Koreguotas"
    ),
    p_label = paste0("PD = ", format(p, decimal.mark = ",", scientific = FALSE))
  )

make_correction_coverage_plot <- function(data, pd_value) {
  
  data %>%
    filter(p == pd_value) %>%
    ggplot(
      aes(
        x = rho,
        y = coverage,
        color = method,
        linetype = interval_type,
        group = interaction(method, interval_type)
      )
    ) +
    geom_hline(
      yintercept = 1 - alpha,
      linetype = "dashed",
      linewidth = 0.5,
      color = "gray35"
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.3) +
    scale_x_continuous(
      breaks = sort(unique(data$rho)),
      labels = function(x) gsub("\\.", ",", sprintf("%.2f", x))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1),
      labels = function(x) format(x, scientific = FALSE, decimal.mark = ",")
    ) +
    labs(
      title = paste0("PD = ", format(pd_value, decimal.mark = ",", scientific = FALSE)),
      x = expression("Turto koreliacija " * rho),
      y = "Empirinė padengimo tikimybė",
      color = "Metodas",
      linetype = "Intervalas"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "gray85"),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "gray85"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

cc1 <- make_correction_coverage_plot(coverage_correction_plot_data, 0.001)
cc2 <- make_correction_coverage_plot(coverage_correction_plot_data, 0.005)
cc3 <- make_correction_coverage_plot(coverage_correction_plot_data, 0.01)
cc4 <- make_correction_coverage_plot(coverage_correction_plot_data, 0.02)
cc5 <- make_correction_coverage_plot(coverage_correction_plot_data, 0.05)

plot_correction_coverage <- (
  (cc1 | cc2 | cc3) /
    (plot_spacer() | cc4 | cc5 | plot_spacer())
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

plot_correction_coverage <- plot_correction_coverage +
  plot_annotation(
    title = "Klasikinių ir efektyviu imties dydžiu koreguotų intervalų padengimo tikimybė, kai n = 100 000",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  )

print(plot_correction_coverage)

ggsave(
  filename = "plots/correction_coverage_n100000.png",
  plot = plot_correction_coverage,
  width = 12,
  height = 7.5,
  dpi = 300
)


#------------Vid intervalo plotis prieš ir po korekcijos grafikas---------------

# Duomenys grafikui
width_correction_plot_data <- correction_comparison %>%
  mutate(
    method = factor(
      method,
      levels = c("Wald", "Clopper_Pearson", "Wilson", "Agresti_Coull"),
      labels = c("Wald", "Clopper-Pearson", "Wilson", "Agresti-Coull")
    )
  ) %>%
  select(p, rho, method, width_classic, width_corrected) %>%
  pivot_longer(
    cols = c(width_classic, width_corrected),
    names_to = "interval_type",
    values_to = "width"
  ) %>%
  mutate(
    interval_type = recode(
      interval_type,
      "width_classic" = "Klasikinis",
      "width_corrected" = "Koreguotas"
    ),
    p_label = paste0("PD = ", format(p, decimal.mark = ",", scientific = FALSE))
  )

make_correction_width_plot <- function(data, pd_value) {
  
  data %>%
    filter(p == pd_value) %>%
    ggplot(
      aes(
        x = rho,
        y = width,
        color = method,
        linetype = interval_type,
        group = interaction(method, interval_type)
      )
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.3) +
    scale_x_continuous(
      breaks = sort(unique(data$rho)),
      labels = function(x) gsub("\\.", ",", sprintf("%.2f", x))
    ) +
    scale_y_continuous(
      labels = function(x) format(
        x,
        scientific = FALSE,
        decimal.mark = ",",
        trim = TRUE
      )
    ) +
    labs(
      title = paste0("PD = ", format(pd_value, decimal.mark = ",", scientific = FALSE)),
      x = expression("Turto koreliacija " * rho),
      y = "Vidutinis intervalo plotis",
      color = "Metodas",
      linetype = "Intervalas"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "gray85"),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "gray85"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

cw1 <- make_correction_width_plot(width_correction_plot_data, 0.001)
cw2 <- make_correction_width_plot(width_correction_plot_data, 0.005)
cw3 <- make_correction_width_plot(width_correction_plot_data, 0.01)
cw4 <- make_correction_width_plot(width_correction_plot_data, 0.02)
cw5 <- make_correction_width_plot(width_correction_plot_data, 0.05)

plot_correction_width <- (
  (cw1 | cw2 | cw3) /
    (plot_spacer() | cw4 | cw5 | plot_spacer())
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

plot_correction_width <- plot_correction_width +
  plot_annotation(
    title = "Klasikinių ir efektyviu imties dydžiu koreguotų intervalų plotis, kai n = 100 000",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  )

print(plot_correction_width)

ggsave(
  filename = "plots/correction_width_n100000.png",
  plot = plot_correction_width,
  width = 12,
  height = 7.5,
  dpi = 300
)