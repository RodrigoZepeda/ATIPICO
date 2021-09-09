rm(list = ls())
library(readxl)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(scales)
library(viridis)
library(foreach)
library(loo)

setwd("~/Dropbox/GUIAS_VSTAN/GuiasITT_V3")

source("functions/read_file.R")
source("functions/modelo_edad_sexo.R")

set.seed(1322121)
stan_seed <- sample(1:10000000, 1)


#Este es el camino al compilador puedes NO ponerlo
#y correr compilar_modelo() y en automático detecta.
#Yo pongo mi path porque uso otro compilador no el default de mac
#puedes ponerle NULL y usa el automático
cpath <-  "/usr/local/opt/llvm/bin/clang++" #NULL
modelo_gamma_un_sexo          <- compilar_modelo(model = "models/Modelo_Gamma_Outliers_Log_Edad_sin_Sexo.stan", compiler_path_cxx = cpath)
modelo_gamma                  <- compilar_modelo(model = "models/Modelo_Gamma_Log_Lineal.stan", compiler_path_cxx = cpath)
modelo_sin_outliers           <- compilar_modelo(model = "models/Modelo_Gamma_sin_Outliers_Log_Edad.stan", compiler_path_cxx = cpath)
modelo_sin_outliers_sin_sexo  <- compilar_modelo(model = "models/Modelo_sin_outliers_Log_edad_sin_sexo.stan", compiler_path_cxx = cpath)
modelo_sin_edad               <- compilar_modelo(model = "models/Modelo_sin_edad.stan", compiler_path_cxx = cpath)
modelo_basico                 <- compilar_modelo(model = "models/Modelo_Gamma.stan", compiler_path_cxx = cpath)
modelo_basico_con_outliers    <- compilar_modelo(model = "models/Modelo_Gamma_Outliers.stan", compiler_path_cxx = cpath)
#modelo_gp_con_outliers        <- compilar_modelo(model = "models/Modelo_Gamma_Log_Edad_GP.stan", compiler_path_cxx = cpath)

#Leemos los diagnósticos que existen
diags <- read_excel("datasets/Diagnósticos_Guias_f.xlsx") %>% 
  mutate(`Diagnóstico exclusivo de mujer` = 
           (`Diagnóstico exclusivo de mujer` == "x")) %>%
  mutate(`Diagnóstico exclusivo de hombre` = 
           (`Diagnóstico exclusivo de hombre` == "x")) %>%
  mutate(`Diagnóstico exclusivo de mujer` = 
           replace_na(`Diagnóstico exclusivo de mujer`, FALSE)) %>%
  mutate(`Diagnóstico exclusivo de hombre` = 
           replace_na(`Diagnóstico exclusivo de hombre`, FALSE)) 

#Lectura de los archivos
fnames    <- list.files("datasets/bases_waldo_agosto_2021","*.csv", full.names = T)

#Edades a interpolar
edades_interpol <- c(19, seq(20,60, by = 1),61)

#Hace dos corridas: una para aquellos con menos de 80,000 observaciones
#que mi compu sí logra paralelizar en terminal (no en RStudio) y una
#segunda que va de uno por uno para aquellas con >80k observaciones
#porque mi compu muere su RAM
for (intento in c("paralelo_memoria_limitada","serie_memoria_completa")){

  #Si ya hay archivos previos que convergen
  files_convergen <- str_replace(list.files("logs"),"log_","")
  files_convergen <- paste0("datasets/bases_waldo_agosto_2021/b_", files_convergen)
  files_interesan <- paste0("datasets/bases_waldo_agosto_2021/b_",unique(diags$`CIE 10`),".csv")
  if (length(files_convergen) > 0){
    fnames <- fnames[which(fnames %in% files_interesan)]
    fnames <- fnames[which(!(fnames %in% files_convergen))]
  }
  #fnames <- fnames[-1]
  
  #Loopeamos por cada código
  n <- length(fnames)
  
  #Crear los directorios por si no están
  if (!dir.exists("results")){dir.create("results")}
  if (!dir.exists("logs")){dir.create("logs")}
  if (!dir.exists("fitted_plots")){dir.create("fitted_plots")}
  if(file.exists("log.out")){file.remove("log.out")}
  
  if (intento == "paralelo_memoria_limitada"){
    #Máximos intentos para decir que converge
    max_attempts    <- 10
    
    #Aquellos que no matan la memoria de mi compu
    mlimit <- 80000 #En 80k mi computadora aguanta súper bien
    cores  <- max(floor(parallel::detectCores()/2),4)
    cl     <- parallel::makeCluster(cores, outfile = "log.out")
    doParallel::registerDoParallel(cl)
    options(mc.cores = 2)
    `%scooby_do%` <- `%dopar%`
    
  } else {
    #Máximos intentos para decir que converge
    max_attempts <- 5
    
    mlimit <- Inf #Ya no es necesario limitar la computadora
    options(mc.cores = cores)
    `%scooby_do%` <- `%do%` 
  }
  
  #ISMAEL:
  #PARA CORRER EN PARALELO CAMBIA EL do POR dopar
  foreach (i = 1:n, .combine='c', .inorder=FALSE, 
           .packages=c("cmdstanr","scales","viridis","tidyverse","posterior")) %scooby_do% {
    
    #Selección del archivo         
    file_name <- fnames[i]
    
    #Convergencia
    convergencia <- FALSE
    
    #Lectura de la base de datos
    data_cie10 <- read_file(file_name, valid_list = unique(diags$`CIE 10`))
    message(data_cie10$icd_code)
    
    if (!is.null(data_cie10$data_cie10) & nrow(data_cie10$data_cie10) < mlimit){
      
      #Checamos que el diagnóstico no sea exclusivo de hombre o mujer
      row_diag <- diags %>% filter(`CIE 10` == data_cie10$icd_code)
      if (row_diag$`Diagnóstico exclusivo de mujer`){
        data_cie10$data_cie10 <- data_cie10$data_cie10 %>% filter(sexo == "Mujer")
        data_cie10$exclusivo  <- "Mujer"
        modelo_outliers       <- modelo_gamma_un_sexo
        modelo_no_outliers    <- modelo_sin_outliers_sin_sexo
      } else if (row_diag$`Diagnóstico exclusivo de hombre`){
        data_cie10$data_cie10 <- data_cie10$data_cie10 %>% filter(sexo == "Hombre")
        data_cie10$exclusivo  <- "Hombre"
        modelo_outliers       <- modelo_gamma_un_sexo
        modelo_no_outliers    <- modelo_sin_outliers_sin_sexo
      } else {
        data_cie10$exclusivo  <- "Todos"
        modelo_outliers       <- modelo_gamma
        modelo_no_outliers    <- modelo_sin_outliers
      }
      
      #Estimaciones con outliers (si la prop de outliers es muy baja no converge)
      attempts_outliers <- 1
      while(!convergencia & attempts_outliers <= max_attempts){
        try({
          init         <- inicializar_modelo(data_cie10$data_cie10, edades_interpol = edades_interpol)
          modelo       <- ajusta_modelo(modelo_outliers, init, icd_code = data_cie10$icd_code,
                                        stan_seed = stan_seed, method = "variational")
          convergencia <- (length(modelo$output_files()) > 0)
        })
        attempts_outliers <- attempts_outliers + 1
        if (convergencia){outliers     <- T; edad <- T; modelo_tipo <- "Edad_sexo_outliers"}
      }
  
      #Estimaciones x edad sin outliers (si la prop de outliers es muy baja converge)
      attempts_outliers <- 1
      while(!convergencia & attempts_outliers <= max_attempts){
        try({
          init         <- inicializar_modelo(data_cie10$data_cie10, edades_interpol = edades_interpol)
          modelo       <- ajusta_modelo(modelo_no_outliers, init, icd_code = data_cie10$icd_code,
                                        stan_seed = stan_seed, method = "variational")
          convergencia <- (length(modelo$output_files()) > 0)
        })
        attempts_outliers <- attempts_outliers + 1
        if (convergencia){outliers     <- F; edad <- T; modelo_tipo <- "Edad_sexo_sin_outliers"}
      }
      
      #Modelo sin edad con sexo
      attempts_outliers <- 1
      while(!convergencia & attempts_outliers <= max_attempts & data_cie10$exclusivo == "Todos"){
        try({
          init         <- inicializar_modelo(data_cie10$data_cie10, edades_interpol = edades_interpol)
          modelo       <- ajusta_modelo(modelo_basico_con_outliers, init, icd_code = data_cie10$icd_code, 
                                        stan_seed = stan_seed, method = "variational")
          convergencia <- (length(modelo$output_files()) > 0)
        })
        attempts_outliers <- attempts_outliers + 1
        if (convergencia){outliers     <- T; edad <- F; modelo_tipo <- "Sexo_outliers"}
      }
  
      #Modelo básico con outliers
      attempts_outliers <- 1
      while(!convergencia & attempts_outliers <= max_attempts){
        try({
          init         <- inicializar_modelo(data_cie10$data_cie10, edades_interpol = edades_interpol)
          modelo       <- ajusta_modelo(modelo_basico_con_outliers, init, icd_code = data_cie10$icd_code, 
                                        stan_seed = stan_seed, method = "variational")
          convergencia <- (length(modelo$output_files()) > 0)
        })
        attempts_outliers <- attempts_outliers + 1
        if (convergencia){outliers     <- T; edad <- F; modelo_tipo <- "Basico_outliers"}
      }
      
      # #Estimaciones sencillas sin edad ni sexo ni outliers
      attempts_outliers <- 1
      while(!convergencia & attempts_outliers <= max_attempts){
        try({
          init         <- inicializar_modelo(data_cie10$data_cie10, edades_interpol = edades_interpol)
          modelo       <- ajusta_modelo(modelo_basico, init, icd_code = data_cie10$icd_code, 
                                        stan_seed = stan_seed, method = "variational")
          convergencia <- (length(modelo$output_files()) > 0)
        })
        attempts_outliers <- attempts_outliers + 1
        if (convergencia){outliers     <- F; edad <- F;  modelo_tipo <- "Basico_sin_outliers"}
      }
      
      if (convergencia){
        message("Convergió!")
        #Obtener resultados si convergió
        result       <- resultados_modelo(modelo, edades_interpol = edades_interpol, 
                                          exclusivo = data_cie10$exclusivo, 
                                          outliers = outliers, edad = edad)
        save(result, file = paste0("results/Sims",data_cie10$icd_code,".rds"))
        
        #Excels
        message("Resumiendo resultados")
        resume_resultados(result, data_cie10, modelo_tipo = modelo_tipo)
        
        #Gráficas
        message("Graficando")
        grafica_resultados(result, data_cie10$icd_code, row_diag$Diagnóstico)
        message(file_name)
      } else {
        message(paste0(file_name, "-- no convergió"))
      }
    }
  }
  
  #SI HACES EN PARALELO MATA EL CLUSTER:
  if (intento == "paralelo_memoria_limitada"){
    parallel::stopCluster(cl)
    message("Terminé de correr en paralelo")
  } else {
    message("Terminé de correr todo")
  }
}
message("VERIFICA LOS 470 RESULTADOS MANUALMENTE PARA CHECAR QUE NO HAYA HABIDO ALGUNO QUE CON N MUY GRANDE HIZO EL MODELO BÁSICO LO QUE IMPLICARÍA UN MODELO MAL ESPECIFICADO.")
