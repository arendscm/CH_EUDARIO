library(ggthemes)

##Themes
my_theme <- function(base_size=12, base_family="Roboto") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           legend.direction = "vertical",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing.y = unit(0.2, "cm"),
           legend.title = element_text(),
           plot.margin=unit(c(10,5,5,5),"mm"),
           panel.spacing.x = unit(1,"line"),
           strip.background=element_rect(colour=NA,fill=NA)
           #strip.text = element_text(face="bold")
   ))
}


##Percent function
percent <- function(x, digits = 0, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

theme_graphicalabstract <- function(base_size=12, base_family="Roboto") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           legend.direction = "vertical",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing.y = unit(0.2, "cm"),
           legend.title = element_text(),
           plot.margin=unit(c(10,5,5,5),"mm"),
           panel.spacing.x = unit(1,"line"),
           strip.background=element_rect(colour=NA,fill=NA)
           #strip.text = element_text(face="bold")
   ))
  
}
