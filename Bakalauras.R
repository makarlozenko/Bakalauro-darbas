library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(123)

################################################################################
#                         Simuliacijos parametrai


pd_values <- c(0.001, 0.005, 0.01, 0.02, 0.05)

n_values <- c(500, 1000, 2000, 5000, 10000) # dar bus pataisymas

rho_values <- c(0.01, 0.05, 0.10, 0.20, 0.25)

alpha <- 0.05

B <- 10000

# Kursim po 200 portfelių, kad greičiau veiktų skriptas
chunk_size <- 200


# Limituojam galimus intervalus iki - [0, 1]
clip_01 <- function(x) {
  pmin(1, pmax(0, x))
}

################################################################################
#                Funkcijos pasikliovimo intervalamas

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


################################################################################
#                     Nepriklausomų įvykių modelis

simulate_independent_counts <- function(B, n, p) {
  rbinom(n = B, size = n, prob = p)
}

################################################################################
#                  Priklausomų įvykių modelis Vasicek su ARMA(2,1)

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
    
    # Individualus veiksnys kiekvienam klientui generuojamas pagal ARMA(2,1)
    z <- matrix(NA, nrow = current_B, ncol = n)
    
    for (b in 1:current_B) {
      z[b, ] <- as.numeric(arima.sim(
        model = list(ar = arma_ar, ma = arma_ma),
        n = n
      ))
    }
    
    # Kad Z_i būtų palyginamas su standartiniu normaliuoju veiksniu
    z <- scale(z)
    
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


################################################################################
#                        Generuojam defoltų kiekį

simulate_default_counts <- function(B, n, p, rho = 0,model = "independent", chunk_size = 200) {
  
  if (model == "independent") {
    return(simulate_independent_counts(B = B, n = n, p = p))
  }
  if (model == "vasicek") {
    return(simulate_vasicek_counts( B = B, n = n, p = p,rho = rho, chunk_size = chunk_size))
  }
  stop("Nežinomas įvestas modelis")
}
################################################################################
#                         Tikrinam vieną scenarijų


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
    # Skaičiuojam palyginimo metrikas--------------
    coverage <- mean(ci$lower <= p & p <= ci$upper)
    avg_width <- mean(interval_width)
    median_width <- median(interval_width)
    lower_zero_rate <- mean(ci$lower == 0)
    #----------------------------------------------
    
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


################################################################################
################################################################################
#  Testuojame
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
################################################################################
################################################################################
#                  Pradedame Monte Carlo simuliaciją


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
# Simuliacijos rezultatai lieka faile
write.csv(
  results,
  file = "results/ci_simulation_results.csv",
  row.names = FALSE
)
print(head(results))


