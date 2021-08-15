#Función para leer la base de datos enviada por
#Prestaciones Económicas. 
#Optimizada para leer una base con ese formato.
read_file <- function(file_name, 
                      icd_code = str_replace(str_replace(file_name,".*b_",""),".csv",""),
                      valid_list = NULL){
  
  #Checamos que la enfermedad exista en la lista de validación
  data_cie10   <- NULL
  if (!is.null(valid_list) && (icd_code %in% valid_list)){
    data_cie10   <- read_csv(file_name, locale = locale(encoding = "UTF-8"),
                             col_types = cols(
                               EDAD                  = col_double(),
                               sexo                  = col_character(),
                               dias_acum             = col_double(),
                               CODCIE10              = col_character(),
                               NUM_AFILIACION        = col_character(),
                               NUM_FOLIO_INCAPACIDAD = col_character(),
                               CVE_DELEGACION_EXP    = col_double(),
                               CVE_DELEGACION_TRAM   = col_double(),
                               CVE_UNIDAD_EXPEDICION = col_character(),
                               CVE_UNIDAD_TRAMITE    = col_character(),
                               CVE_REG_PATRONAL      = col_character(),
                               colaps_fec_inicio     = col_date("%d/%m/%y"),
                               colaps_fec_termino    = col_date("%d/%m/%y"),
                               FEC_PROCESO           = col_date("%d/%m/%y"),
                               FEC_EXPEDICION        = col_date("%d/%m/%y"),
                               DIAGNOSTICO           = col_character()
                             ))
    data_cie10 <- data_cie10 %>% 
      mutate(Sexo = if_else(sexo == "Hombre", 0, 1)) %>%
      select(EDAD, Sexo, sexo, dias_acum)
  } 
  return(list(data_cie10 = data_cie10, icd_code = icd_code))
}
