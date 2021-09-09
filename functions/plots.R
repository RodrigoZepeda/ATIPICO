plot_hist <- function(datos, original_data, sexname = "Hombre", 
                      size = 4, ciename = "sincie",
                      breaks_edades = c(0,unique(datos$edad),Inf)){

  
  original_data$data_cie10 <- original_data$data_cie10 %>%
    mutate(edad_cat = cut(EDAD, breaks_edades, right = FALSE)) %>%
    filter(sexo == !!sexname)
  
  if (nrow(original_data$data_cie10) > 0){
    edad_count <- original_data$data_cie10 %>% 
      group_by(edad_cat) %>% tally()
    
    #Histograma de hombres
    plotsex <- ggplot(original_data$data_cie10) +
      geom_histogram(aes(x = fct_rev(edad_cat)), stat = "count", 
                     fill = "#134E39") +
      theme_classic() +
      coord_flip(clip = "off") +
      geom_text(aes(x = edad_cat, y = n, label = scales::comma(n, accuracy = 1)), hjust = -0.1,
                data = edad_count, color = "#404040", size = size) +
      theme(plot.margin = unit(c(1, 5, 1, 1), "lines")) +
      theme(axis.line.x  = element_blank(),
            axis.text.x  = element_blank(), 
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y  = element_text(color = "#404040", size = 18),
            axis.ticks.y = element_line(color = "#404040"),
            axis.line.y  = element_line(color = "#404040")) 
  } else {
    plotsex <- ggplot() + theme_void()
  }
    ggsave(paste0(sexname,"_",ciename,".pdf"), plotsex, width = 4, height =6)
}
