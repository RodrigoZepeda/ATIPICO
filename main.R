rm(list = ls())
library(readxl)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(scales)
library(viridis)
library(foreach)

source("functions/read_file.R")
source("functions/modelo_edad_sexo.R")

set.seed(7345789)
stan_seed <- sample(1:10000000, 1)


#Este es el camino al compilador puedes NO ponerlo
#y correr compilar_modelo() y en automático detecta.
#Yo pongo mi path porque uso otro compilador no el default de mac
cpath <-  "/usr/local/opt/llvm/bin/clang++"
modelo_gamma         <- compilar_modelo(compiler_path_cxx = cpath)
modelo_gamma_un_sexo <- compilar_modelo(model = "models/Modelo_Gamma_Outliers_Edad_sin_Sexo.stan", compiler_path_cxx = cpath)

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
fnames    <- list.files("datasets/bases_output_enero_2021","*.csv", full.names = T)
file_name <- "datasets/bases_output_enero_2021/b_N760.csv"

#Edades a interpolar
edades_interpol <- seq(20, 60, by = 5)

#Lista de los que no convergen
convergen <- c()

#Loopeamos por cada código
n <- length(fnames)

cores <- max(parallel::detectCores() - 2,1)
cl    <- parallel::makeCluster(cores/2)
doParallel::registerDoParallel(cl)
options(mc.cores = 2)

foreach (i = 1:n, .combine='c', 
         .packages=c("cmdstanr","scales","viridis","tidyverse","posterior")) %dopar% {
  
  file_name <- fnames[i]
  message(file_name)
  
  #Lectura de la base de datos
  data_cie10 <- read_file(file_name, valid_list = unique(diags$`CIE 10`))
  if (!is.null(data_cie10$data_cie10)){
    
    #Checamos que el diagnóstico no sea exclusivo de hombre o mujer
    row_diag <- diags %>% filter(`CIE 10` == data_cie10$icd_code)
    if (row_diag$`Diagnóstico exclusivo de mujer`){
      data_cie10$data_cie10 <- data_cie10$data_cie10 %>% filter(sexo == "Mujer")
      data_cie10$exclusivo  <- "Mujer"
      modelo <- modelo_gamma_un_sexo
    } else if (row_diag$`Diagnóstico exclusivo de hombre`){
      data_cie10$data_cie10 <- data_cie10$data_cie10 %>% filter(sexo == "Hombre")
      data_cie10$exclusivo  <- "Hombre"
      modelo <- modelo_gamma_un_sexo
    } else {
      data_cie10$exclusivo  <- "Todos"
      modelo <- modelo_gamma
    }
    
    #Estimaciones
    try({
      init      <- inicializar_modelo(data_cie10$data_cie10, edades_interpol = edades_interpol)
      modelo    <- ajusta_modelo(modelo, init, stan_seed = stan_seed)
      result    <- resultados_modelo(modelo, edades_interpol = edades_interpol, exclusivo = data_cie10$exclusivo)
      convergen <- c(convergen, data_cie10$icd_code)
      #Excels
      if (!dir.exists("results")){dir.create("results")}
      if (!dir.exists("logs")){dir.create("logs")}
      resume_resultados(result, data_cie10)
      
      #Gráficas
      if (!dir.exists("fitted_plots")){dir.create("fitted_plots")}
      grafica_resultados(result, data_cie10$icd_code, row_diag$Diagnóstico)
    })
  }
  convergen
}

parallel::stopCluster(cl)

