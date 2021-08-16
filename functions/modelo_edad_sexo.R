
compilar_modelo <- function(model = "models/Modelo_Gamma_Outliers_Edad.stan",
                            compiler_path_cxx = NULL, 
                            cxx_flags = "-O3 -march=native", 
                            use_stan_threads = T){
  
  message("Compilando modelo en C++")
  if (!is.null(compiler_path_cxx)){
    cpp_options <- list(cxx_flags = cxx_flags, 
                        cxx = compiler_path_cxx, 
                        stan_threads = use_stan_threads)
  } else {
    cpp_options <- list(cxx_flags = cxx_flags, stan_threads = use_stan_threads)
  }
  
  cmdstan_model(model, cpp_options = cpp_options)
  
}

inicializar_modelo <- function(data_cie10, edades_interpol = c(20,30,40,50,60)){
  
  #Datos para inicializar el modelo
  datos_init <- list(N = nrow(data_cie10), 
                     Tobs = as.vector(data_cie10$dias_acum), 
                     Edad = as.vector(data_cie10$EDAD), 
                     Sexo = as.vector(data_cie10$Sexo),
                     EdadesResultado = as.vector(edades_interpol), 
                     M = length(edades_interpol))
  
  #Parámetros iniciales para MCMC
  params_init <- function(chain_id = 1) {
    list(error1       = runif(nrow(data_cie10), 0, 3.5), 
         error2       = runif(nrow(data_cie10), 0, 1000), 
         beta         = (rnorm(2, 2.5, 1) %>% abs()) + 1,
         gamma        = rnorm(2, 0, 0.1),
         eta          = rnorm(2, 0, 0.1),
         nu           = rnorm(2, 100, 0.1) %>% abs(),
         theta        = runif(1,0,1)
    )}
  
  return(list(datos_init = datos_init, params_init = params_init))
}

ajusta_modelo <- function(modelo, init, 
                          method = "variational",
                          output_dir = "cmdstan",
                          stan_seed = 237235,
                          hmc_params = list(
                            chains = 4,
                            seed = 2457924,
                            iter_warmup = 250,
                            adapt_delta = 0.95,
                            iter_sampling = 250,
                            max_treedepth = 2^(11),
                            threads_per_chain = 4
                          ),
                          variational_params = list(
                            output_samples = 1000,
                            adapt_engaged = T,
                            adapt_iter = NULL,
                            threads = 2,
                            iter_variational = 100000,
                            sig_figs = 4
                          )){
  
  #Crear directorio de outputs
  if (!dir.exists(output_dir)){dir.create(output_dir)}
  
  #Ver que el cálculo de chains y threads per chain tenga sentido  
  tph <- parallel::detectCores()/hmc_params$chains
  if (hmc_params$threads_per_chain > tph & method == "HMC"){
    stop("Ajusta chains o threads per chain pues no tienes tantos núcleos.")
  } 
  
  #Este es el mejor método pero es más lento
  if (method == "HMC"){
    model_sample <- modelo$sample(data = init$datos_init, 
                                  chains = hmc_params$chains, 
                                  seed = stan_seed, 
                                  iter_warmup = hmc_params$iter_warmup,
                                  adapt_delta = hmc_params$adapt_delta, 
                                  iter_sampling = hmc_params$iter_sampling,
                                  init = init$params_init,
                                  output_dir = output_dir, 
                                  max_treedepth = hmc_params$max_treedepth,
                                  threads_per_chain = hmc_params$threads_per_chain)
  } else if (method == "variational"){
    model_sample <- modelo$variational(data = init$datos_init, 
                                       seed = stan_seed, 
                                       sig_figs = variational_params$sig_figs,
                                       iter= variational_params$iter_variational,
                                       init = init$params_init,
                                       adapt_iter = variational_params$adapt_iter,
                                       adapt_engaged = variational_params$adapt_engaged,
                                       output_samples	= variational_params$output_samples,
                                       threads = variational_params$threads,
                                       output_dir = output_dir)
  } else {
    message(paste0("Method ", method, " not found. Try 'HMC' or 'variational'"))
  }
  
  return(model_sample)
}

resultados_modelo <- function(model_sample, edades_interpol = c(20,30,40,50,60),
                              exclusivo = "Todos"){
  
  #Obtenemos la distribución posterior
  suppressWarnings({
    ppdist_0  <- model_sample$draws(variables="Tpred_0") %>% as_draws_df()
    ppdist_0  <- ppdist_0 %>% select(-starts_with("."))
  })
  colnames(ppdist_0) <- edades_interpol
  
  if (exclusivo == "Todos"){
    suppressWarnings({
      ppdist_1  <- model_sample$draws(variables="Tpred_1") %>% as_draws_df()
      ppdist_1  <- ppdist_1 %>% select(-starts_with("."))
    })
    colnames(ppdist_1) <- edades_interpol
    ppdist_0  <- ppdist_0 %>% mutate("Sexo" = "Hombre")
    ppdist_1  <- ppdist_1 %>% mutate("Sexo" = "Mujer")
    ppdist    <- rbind(ppdist_0, ppdist_1)
  } else {
    ppdist_0  <- ppdist_0 %>% mutate("Sexo" = !!exclusivo)
    ppdist    <- ppdist_0
  }
  
  ppdist    <- ppdist %>% 
    pivot_longer(cols =  as.character(edades_interpol), 
                 values_to = "Tiempo", names_to = "Edad")
  ppdist    <- ppdist %>% mutate(Edad = paste0("Edad = ", Edad)) 
  
  #Proporción estimada de outliers (1 - theta)
  prop_true <- model_sample$draws(variables="theta") %>% summarise_draws()
  
  return(list(results = ppdist, prop_validos = prop_true))
}

#Genera los csv de resumen
resume_resultados <- function(result, 
                              data_cie10,
                              quantiles = c(0.25, 0.5, 0.75),
                              qnames = c("Ligera","Moderada","Pesada"),
                              fname_results = paste0("results/",data_cie10$icd_code,".csv"),
                              fname_logs = paste0("logs/log_",data_cie10$icd_code,".csv")){
  
  #Obtenemos los cuantiles de interés
  resumen <- result$results %>% 
    group_by(Edad, Sexo) %>%
    summarise(`Duración de la ITT (días)` = quantile(Tiempo, quantiles), .groups = "keep")
  edades <- unique(result$results$Edad)
  sexos  <- unique(result$results$Sexo)
  resumen$`Carga de trabajo` <- rep(rep(qnames, length(sexos)), length(edades))
  resumen <- resumen %>% arrange(Sexo, `Carga de trabajo`, Edad)
  write_excel_csv(resumen, fname_results)
  
  #Obtenemos el log:
  sample_size  <- c("Muestra total (n)", scales::comma(nrow(data_cie10$data_cie10)))
  sex_size     <- data_cie10$data_cie10 %>% 
    group_by(sexo) %>% 
    summarise(n = scales::comma(n()))
  age_range    <- data_cie10$data_cie10 %>% 
    group_by(sexo) %>% 
    summarise(`Rango` = paste0("[",min(EDAD),", ", max(EDAD),"]"), .groups = "keep")
  itt_range    <- data_cie10$data_cie10 %>% 
    group_by(sexo) %>% 
    summarise(`Rango` = paste0("[",min(dias_acum),", ", max(dias_acum),"]"), .groups = "keep")
  prop_validos <- c("Proporción de datos válidos (no atípicos)", 
                    paste0("[",
                         scales::percent(as.numeric(result$prop_validos["q5"]),0.01),
                         ",", 
                         scales::percent(as.numeric(result$prop_validos["q95"],0.01)),
                         "]"))
  datos               <- rbind(sample_size, prop_validos, sex_size)
  colnames(datos)     <- c("Variable","Tamaño/rango de valores")
  age_range$sexo      <- paste0("Rango de edad en ",age_range$sexo)
  colnames(age_range) <- colnames(datos)
  datos               <- rbind(datos, age_range)
  
  itt_range$sexo      <- paste0("Rango registrado de duración de una ITT en ",itt_range$sexo)
  colnames(itt_range) <- colnames(datos)
  datos               <- rbind(datos, itt_range)
  write_excel_csv(datos, fname_logs)
  
}

grafica_resultados <- function(resultados_modelo, 
                               icd_code = "SIN CÓDIGO",
                               diagmed  = "SIN DIAGNÓSTICO",
                               title = paste0(icd_code,": ", diagmed),
                               fname = paste0("fitted_plots/", icd_code,".png"), 
                               edades_interpol = c(20,30,40,50,60)){
  
  ggplot(as_tibble(resultados_modelo$results)) +
    geom_density(aes(x = Tiempo, color = Edad)) +
    facet_wrap(~Sexo) +
    theme_classic() +
    scale_color_viridis(discrete = T) + 
    labs(
      x = "Duración de la incapacidad temporal en el trabajo", 
      y = "",
      title = title,
      subtitle = "ATIPICO: Ajuste del TIempo de Prescripción de Incapacidades Con Outliers",
      caption  = paste0("Instituto Mexicano del Seguro Social | Datos de la Dirección de Prestaciones Económicas | ", lubridate::today())
    ) +
    scale_y_continuous(labels = scales::percent)
  ggsave(fname, width = 12, height = 6, dpi = 750)
}
