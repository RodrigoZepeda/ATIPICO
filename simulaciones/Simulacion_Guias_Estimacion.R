#Script para simular la estructura de datos y verificar
#si el modelo funciona para recuperarlos
rm(list = ls())

library(tidyverse)
library(cmdstanr)
library(scales)
library(posterior)
library(kdensity)
library(bayesplot)
library(mixdist)

set.seed(23476785)
stan_seed <- 23749
chains = 4; iter_warmup = 250; nsim = 500; pchains = 4; 
threads_per_chain = 4; threads = 8; iter_variational = 10000
method = "variational" #faster (less accurate option) = "variational"

#Flags for my compiler faster
compiler_path_cxx <- "/usr/local/opt/llvm/bin/clang++"
options(mc.cores = parallel::detectCores())

#Simulation distribution
#gamma para que el modelo esté bien especificado;
#weibull para que no esté bien especificado
simdist <- "gamma" #weibull ó gamma
media    <- 10
varianza <- 20

#No todos los datos están redondeados así que se establece un porcentaje
#de cuántos tienen redondeo
perc_redondeado <- 0.75  #Porcentaje de datos redondeados a 7
perc_outliers   <- 0.05  #Se agrega el 1% de outliers
nsims           <- 10000 #Número de simulaciones para el modelo

#Generamos las simulaciones:
if (simdist == "gamma"){
  message("Distribución Gamma")
  scale        <- varianza/media  %>% round(.,2)
  shape        <- media/scale     %>% round(.,2)
  col_redondeo <- "orange"
  tipo_modelo  <- paste0("bien especificado pues T~Gamma(", 
                         shape, ",", scale,")")
  sample_dist  <- function(){rgamma(nsims, shape = shape, scale = scale)}
  dist_fun     <- function(x){dgamma(x, shape = shape, scale = scale)}
} else if (simdist == "weibull"){
  message("Distribución Weibull")
  params       <- weibullpar(media, sqrt(varianza))
  shape        <- params["shape"] %>% as.numeric() %>% round(.,2)
  scale        <- params["scale"] %>% as.numeric() %>% round(.,2)
  col_redondeo <- "purple"
  tipo_modelo  <- paste0("mal especificado pues T~Weibull(", 
                         shape, ",", scale,")")
  sample_dist  <- function(x){rweibull(nsims, shape = shape, scale = scale)}
  dist_fun     <- function(x){dweibull(x, shape = shape, scale = scale)}
} else {
  stop("Distribución inválida selecciona 'gamma' o 'weibull'.")
}

#Simulamos los verdaderos datos filtrados como días (enteros)
datos_distribucion <- sample_dist()
datos_distribucion <- round(datos_distribucion,0)

#Agregamos outliers
outliers_id <- sample(1:nsims, ceiling(perc_outliers*nsims))
outliers    <- exp(rnorm(ceiling(perc_outliers*nsims), log(50), 1))
datos_distribucion[outliers_id] <- outliers

#Redondeamos la mayoría de los datos a múltiplos de 7 a partir de 3
redondear   <- datos_distribucion[datos_distribucion > 3]
a_redondear <- sample(ceiling(perc_redondeado*length(redondear)))
redondear[a_redondear]                     <- round(redondear[a_redondear]/7)*7
datos_distribucion[datos_distribucion > 3] <- redondear

#Modelo
message("Fitting model. Go grab a coffee this will take A LOT")
if (!is.null(compiler_path_cxx)){
  cpp_options <- list(cxx_flags = "-O3 -march=native", 
                      cxx = compiler_path_cxx, stan_threads = TRUE)
} else {
  cpp_options <- list(cxx_flags = "-O3 -march=native", 
                      stan_threads = TRUE)
}

model_gamma_v1 <- cmdstan_model("models/Modelo_Gamma_Outliers.stan",
                                cpp_options = cpp_options)

datos <- list(N = length(datos_distribucion), Tobs = datos_distribucion)

initf2 <- function(chain_id = 1) {
  list(error1       = runif(length(datos_distribucion), 0, 3.5), 
       error2       = runif(length(datos_distribucion), 0, 1000), 
       alpha        = rnorm(2, 2.5, 1) %>% abs(),
       beta         = rnorm(2, 2.5, 1) %>% abs(),
       theta        = runif(1,0,1)
  )}
init_ll <- lapply(1:chains, function(id) initf2(chain_id = id))

if (!dir.exists("cmdstan")){dir.create("cmdstan")}
if (method == "HMC"){
  model_sample <- model_gamma_v1$sample(data = datos, 
                                        chains = chains, 
                                        seed = stan_seed, 
                                        iter_warmup = iter_warmup,
                                        adapt_delta = 0.95, 
                                        iter_sampling = nsim - iter_warmup,
                                        init = initf2,
                                        output_dir = "cmdstan", 
                                        max_treedepth = 2^(11),
                                        threads_per_chain = threads_per_chain)
} else if (method == "variational"){
  model_sample <- model_gamma_v1$variational(data = datos, 
                                        seed = stan_seed, 
                                        iter= iter_variational,
                                        init = initf2,
                                        threads = threads,
                                        output_dir = "cmdstan")
} else {
  message(paste0("Method ", method, " not found. Try 'HMC' or 'variational'"))
}

# Herramienta de diagnósitco para verificar el ajuste
#model_sample$cmdstan_diagnose()

#Obtenemos la distribución posterior
ppdist    <- model_sample$draws(variables="Tpred") %>% as_draws_df()
prop_true <- model_sample$draws(variables="theta") %>% summarise_draws()

#Ajustamos densidades para graficar
#Si tenemos demasiados datos reducimos para la gráfica
if (nsims > 1000){datos_distribucion <- sample(datos_distribucion, 1000)}
x             <- seq(0, max(datos_distribucion), length.out = 1000) 
densidad_real <- function(x){dist_fun(x)}
densidad_obs  <- kdensity(datos_distribucion, kernel = "gaussian")
densidad_pred <- kdensity(ppdist$Tpred + 0.01, start = "gamma", 
                          kernel = "gaussian", support = c(0, Inf))

#Gráfica de los datos
data.frame(x = x, Real = densidad_real(x), Observada = densidad_obs(x), 
           Predicha = densidad_pred(x)) %>%
ggplot() +
  geom_line(aes(x = x, y = Real, color = "Real"), size = 1) +
  geom_line(aes(x = x, y = Predicha, color = "Modelo"), size = 1.5,
            linetype = "dotted") +
  geom_histogram(aes(x = x, y = ..density.., fill = "Redondeado (observado)"),
                 breaks = seq(0,100, by = 1), 
                 data = data.frame(x = datos_distribucion)) +
  annotate("text", x = 55, y = 0.1, label = "Valores atípicos") +
  geom_segment(aes(x = 55, y = 0.09, xend = 55, yend = 0.01),
               arrow = arrow(length = unit(0.1, "cm"))) +
  geom_segment(aes(x = 55, y = 0.09, xend = 49, yend = 0.01),
               arrow = arrow(length = unit(0.1, "cm"))) +
  theme_classic() +
  scale_color_manual("Tiempo de recuperación",
                     values = c("Modelo" = "#BF1363",
                                "Real" = "#39A6A3",
                                "Redondeado (observado)" = col_redondeo)) +
  scale_fill_manual("Tiempo de recuperación",
                    values = c("Modelo" = "#BF1363",
                               "Real" = "#39A6A3",
                               "Redondeado (observado)" = col_redondeo)) +
  labs(
    x = "Tiempo de recuperación", 
    y = "",
    title = "Modelo bayesiano para ajuste de redondeo semanal",
    subtitle = paste0("Ajuste Gamma para outliers con un modelo ", tipo_modelo),
    caption = paste0("Muestra tamaño ", comma(nsims), ".\nAjuste de ", 
                     comma(nsim), " simulaciones con ",
                     percent(perc_redondeado), " de los datos redondeados y ", 
                     percent(perc_outliers), " de valores atípicos.")
  ) +
  coord_cartesian(xlim = c(0, 75)) +
  scale_y_continuous(labels = scales::percent)
ggsave(paste0("images/Atipicos",simdist,".png"), width = 8, 
       height = 4, dpi = 750)
