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
threads_per_chain = 4; threads = 12; iter_variational = 50000;
adapt_iter = 1000;
method = "variational" #faster (less accurate option) = "variational"

#Flags for my compiler faster
compiler_path_cxx <- "/usr/local/opt/llvm/bin/clang++"
options(mc.cores = parallel::detectCores())

#Simulation distribution
#gamma para que el modelo esté bien especificado;
#weibull para que no esté bien especificado
simdist <- "gamma" #weibull ó gamma
media_baseline <- 5
varianza       <- 10
gamma_edad     <- 0.1
beta_sexo      <- -1

#Edades a interpolar
edades_interpol <- c(30,40,50)


#No todos los datos están redondeados así que se establece un porcentaje
#de cuántos tienen redondeo
perc_redondeado <- 0.75  #Porcentaje de datos redondeados a 7
perc_outliers   <- 0.05  #Se agrega el 1% de outliers
nsims           <- 10000 #Número de simulaciones para el modelo

#Simulamos las edades 
edades <- rnorm(nsims, mean = 40, sd = 10)
edades[(edades < 18 | edades > 70)] <- runif(sum((edades < 18 | edades > 70)), 20, 70)

#Simulamos el sexo
psexo  <- 0.4
sexo   <- sample(c(0,1), nsims, replace = T, prob = c(psexo, 1-psexo))

#Generamos los parámetros
genera_media <- function(edad, sexo, sims = nsims, random = T){
  if (random){
    mu <- media_baseline + rnorm(sims, gamma_edad, 0.1)*log(edad) + rnorm(sims, beta_sexo, 0.1)*sexo
    mu[mu < 0] <- media_baseline
  } else {
    mu <- media_baseline + gamma_edad*log(edad) + beta_sexo*sexo
    mu[mu < 0] <- media_baseline
  }
  return(mu)
}
media <- genera_media(edades, sexo)

if (min(media < 0)){
  stop("Distribución no está bien especificada la media debe ser > 0")
}

#Generamos las simulaciones:
if (simdist == "gamma"){
  message("Distribución Gamma")
  scale        <- varianza/media  %>% round(.,2)
  shape        <- media/scale     %>% round(.,2)
  col_redondeo <- "orange"
  tipo_modelo  <- paste0("bien especificado pues T~Gamma(", 
                         shape, ",", scale,")")
  sample_dist  <- function(){
    rgamma(nsims, shape = shape, scale = scale)
  }
  dist_fun     <- function(x, edad, sexo){
    scale <- varianza/genera_media(edad, sexo, 1, random = F) 
    shape <- genera_media(edad, sexo, 1, random = F)/scale
    dgamma(x, shape = shape, scale = scale)
  }
} else if (simdist == "weibull"){
  message("Distribución Weibull")
  params       <- weibullpar(media, sqrt(varianza))
  shape        <- params["shape"]  %>% round(.,2)
  scale        <- params["scale"]  %>% round(.,2)
  col_redondeo <- "purple"
  tipo_modelo  <- paste0("mal especificado pues T~Weibull")
  sample_dist  <- function(x){
    rweibull(nsims, shape = unlist(shape), scale =  unlist(scale))
  }
  dist_fun     <- function(x, edad, sexo){
    params       <- weibullpar(genera_media(edad,sexo, 1, random = F), sqrt(varianza))
    shape        <- params["shape"] %>% round(.,2)
    scale        <- params["scale"] %>% round(.,2)
    dweibull(x, shape = unlist(shape), scale = unlist(scale))
  }
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

model_gamma_v1 <- cmdstan_model("models/Modelo_Gamma_Outliers_Edad.stan",
                                cpp_options = cpp_options)

datos <- list(N = length(datos_distribucion), 
              Tobs = datos_distribucion, Edad = edades, Sexo = sexo,
              EdadesResultado = edades_interpol, M = length(edades_interpol))

initf2 <- function(chain_id = 1) {
  list(error1       = runif(length(datos_distribucion), 0, 3.5), 
       error2       = runif(length(datos_distribucion), 0, 1000), 
       beta         = (rnorm(2, 2.5, 1) %>% abs()) + 1,
       gamma        = rnorm(2, 0, 0.1),
       eta          = rnorm(2, 0, 0.1),
       nu           = rnorm(2, 100, 0.1) %>% abs(),
       theta        = runif(1,0,1)
  )}

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
                                             adapt_iter = adapt_iter,
                                             adapt_engaged = T,
                                             output_samples	= 1000,
                                             threads = threads,
                                             output_dir = "cmdstan")
} else {
  message(paste0("Method ", method, " not found. Try 'HMC' or 'variational'"))
}

# Herramienta de diagnósitco para verificar el ajuste
#model_sample$cmdstan_diagnose()

#Obtenemos la distribución posterior
ppdist_0  <- model_sample$draws(variables="Tpred_0") %>% as_draws_df()
ppdist_1  <- model_sample$draws(variables="Tpred_1") %>% as_draws_df()
ppdist_0  <- ppdist_0 %>% select(-starts_with("."))
ppdist_1  <- ppdist_1 %>% select(-starts_with("."))
colnames(ppdist_0) <- edades_interpol
colnames(ppdist_1) <- edades_interpol
ppdist_0  <- ppdist_0 %>% mutate("Sexo" = "Hombre")
ppdist_1  <- ppdist_1 %>% mutate("Sexo" = "Mujer")
ppdist    <- rbind(ppdist_0, ppdist_1)
ppdist    <- ppdist %>% 
  pivot_longer(cols =  as.character(edades_interpol), 
               values_to = "Tiempo", names_to = "Edad")
ppdist    <- ppdist %>% mutate(Edad = paste0("Edad = ", Edad)) 

#Proporción estimada de outliers
prop_true <- model_sample$draws(variables="theta") %>% summarise_draws()

#Ajustamos densidades para graficar
#Si tenemos demasiados datos reducimos para la gráfica
x             <- seq(0, 100, length.out = 1000) 
for (edad in edades_interpol){
  for (sexo in c(0,1)){
    distribucion <- dist_fun(x, edad, sexo)
    if (edad == edades_interpol[1] & sexo == 0){
      datos_dist <- data.frame(Edad = paste0("Edad = ",edad), Sexo = "Hombre", 
                               Tiempo = distribucion, x = x)
    } else {
      if (sexo == 0){sexname = "Hombre"}else{sexname="Mujer"}
      datos_dist <- data.frame(Edad = paste0("Edad = ",edad), Sexo = sexname, 
                               Tiempo = distribucion, x = x) %>% 
        bind_rows(datos_dist)
    }
  }
}

#Gráfica de los datos
ggplot(datos_dist) +
  geom_line(aes(x = x, y = Tiempo, color = "Real"), size = 1) +
  geom_density(aes(x = Tiempo, color = "Modelo"), data = ppdist) +
  coord_cartesian(xlim = c(0, 40)) +
  facet_grid(Edad ~ Sexo) +
  theme_bw() +
  scale_color_manual("Tiempo de recuperación",
                     values = c("Modelo" = "#BF1363",
                                "Real" = "#39A6A3")) +
  scale_fill_manual("Tiempo de recuperación",
                    values = c("Modelo" = "#BF1363",
                               "Real" = "#39A6A3")) +
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
  scale_y_continuous(labels = scales::percent)
ggsave(paste0("Atipicos_edad_",simdist,".png"), width = 10, height = 10, dpi = 750)
