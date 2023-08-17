#####  Dependencies   #####
library(base)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(readxl)

### read event data fram Sampleregistry
plotdata <- read_excel("data/external/Sample Registry.xlsx", sheet = "events")

### functions and themes
source("src/global_functions_themes.R") #functions and themes

##in months
plotdata%>%
  mutate(tEvent = ifelse(Event == "priorPARPi",-0.2, tEvent))%>%
  mutate(tEvent = ifelse(Event %in% c("Priorline1", "Priorline>1"),-0.1, tEvent))%>%
  mutate(tEvent = ifelse(Event == "Progression",(tEvent -0.15), tEvent))->plotdata2
  

palette = pal_npg("nrc")(9)

ggplot(plotdata2) +
  ## timelength
  geom_segment(
    data = subset(plotdata, Event == "timelength"),
    aes(
      x = 0, xend = tEvent,
      y = reorder(as.character(Patient.ID), tEvent),
      yend = reorder(as.character(Patient.ID), tEvent)
    ),
    size = 1, color = "lightgrey"
  ) +
  ## previous treatment (Priorline1, Priorline>1, priorPARPi)
  geom_point(
    data = subset(plotdata, Event %in% c("Priorline1", "Priorline>1", "priorPARPi")),
    aes(x = tEvent, y = Patient.ID, color = Event)
  ) +
  
  ## EUDARIO Treatment
  geom_bar(
    data = subset(plotdata, Event == "Maintanance"),
    aes(x = tEvent, y = Patient.ID, fill = "EUDARIO Maintenance"), 
    stat = "identity"
  ) +
  geom_bar(
    data = subset(plotdata, Event == "MaintananceStart"),
    aes(x = tEvent, y = Patient.ID, fill = "EUDARIO Chemotherapy"), 
    stat = "identity"
  ) +
  geom_bar(
    data = subset(plotdata, Event == "Chemotherapy"),
    aes(x = tEvent, y = Patient.ID, fill = "EUDARIO Chemotherapy"), 
    stat = "identity"
  ) +
  ## Additional events with different colors, shapes, and sizes
  geom_point(
    data = subset(plotdata, Event %in% c("c1d1_cf","c7d1_cf","ue_cf","eot_cf")),
    aes(x = tEvent, y = Patient.ID, color = type),
    shape = 15, size = 1
  ) +
  geom_point(
    data = subset(plotdata, Event %in% c("c1d1_wb","eot_wb")),
    aes(x = tEvent+0.2, y = Patient.ID, color = type),
    shape = 15, size = 1
  ) +
  geom_point(
    data = subset(plotdata, Event == "Progression"),
    aes(x = tEvent, y = Patient.ID, color = Event),
    shape = 8, size = 0.7
  ) +
  geom_point(
    data = subset(plotdata, Event == "Death"),
    aes(x = tEvent, y = Patient.ID, color = Event),
    shape = 18
  ) +
  geom_point(
    data = subset(plotdata, Event %in% c("PriorLines1","PriorLines>1","PriorPARPi")),
    aes(x = tEvent, y = Patient.ID, color = Event),
    shape = 15, size = 1
  ) +
  ylab("Patients") +
  theme(axis.text.y = element_text(size = 4))+
  xlab("Time in months") +
  scale_fill_manual(values=c("EUDARIO Chemotherapy"=alpha(palette[6], 1), "EUDARIO Maintenance" =alpha(palette[7], 0.8)),name = "Therapy")+
  scale_color_manual(
    values = c(
      "Whole blood" = palette[1], "Plasma" = palette[5], "Progression" = palette[4], "Death" = "black", "PriorLines>1" = palette[2], "PriorLines1"= palette[9],  "PriorPARPi" = palette[3]
    ),
    name = "Event",
    breaks = c("Whole blood", "Plasma", "Progression", "Death" ,"PriorLines>1", "PriorLines1", "PriorPARPi"),
    labels = c("Whole blood", "Plasma", "Progression", "Death" ,"> 1 Prior line", "1 prior line", "Prior PARPi treatment")
  ) +
  guides(
    color = guide_legend(title = "Prior Therapy", override.aes = list(shape = 16))
  ) +
  my_theme() +
  theme(axis.line.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())-> p.timeline

p.timeline

png("output/figures/p.timeline_final.png",width=12, height=8,units="in",res=500,type="cairo")
p.timeline
dev.off()



####for ech patient
plotdata$Patient.ID%>%unique()->pat

plotdata->plotdatabackup
for (current_patient in pat)
{print(current_patient)
  plotdatabackup%>%
    filter(Patient.ID == current_patient)->plotdata
  
  
  palette = pal_npg("nrc")(9)
  
  ggplot(plotdata) +
    ## timelength
    geom_segment(
      data = subset(plotdata, Event == "timelength"),
      aes(
        x = 0, xend = tEvent,
        y = reorder(as.character(Patient.ID), tEvent),
        yend = reorder(as.character(Patient.ID), tEvent)
      ),
      size = 1, color = "lightgrey"
    ) +
    ## previous treatment (Priorline1, Priorline>1, priorPARPi)
    geom_point(
      data = subset(plotdata, Event %in% c("Priorline1", "Priorline>1", "priorPARPi")),
      aes(x = tEvent, y = Patient.ID, color = Event)
    ) +
    
    ## EUDARIO Treatment
    geom_bar(
      data = subset(plotdata, Event == "Maintanance"),
      aes(x = tEvent, y = Patient.ID, fill = "EUDARIO Maintenance"), 
      stat = "identity"
    ) +
    geom_bar(
      data = subset(plotdata, Event == "MaintananceStart"),
      aes(x = tEvent, y = Patient.ID, fill = "EUDARIO Chemotherapy"), 
      stat = "identity"
    ) +
    geom_bar(
      data = subset(plotdata, Event == "Chemotherapy"),
      aes(x = tEvent, y = Patient.ID, fill = "EUDARIO Chemotherapy"), 
      stat = "identity"
    ) +
    ## Additional events with different colors, shapes, and sizes
    geom_point(
      data = subset(plotdata, Event %in% c("c1d1_cf","c1d1_wb","c7d1_cf","ue_cf","eot_cf","eot_wb")),
      aes(x = tEvent, y = Patient.ID, color = type),
      shape = 15, size = 1
    ) +
    geom_text(
      data = subset(plotdata, Event %in% c("c1d1_cf","c1d1_wb","c7d1_cf","ue_cf","eot_cf","eot_wb")),
      aes(x = tEvent, y = Patient.ID, label = Event),
      vjust = -0.5, color = "black"
    ) +
    geom_point(
      data = subset(plotdata, Event == "Progression"),
      aes(x = tEvent, y = Patient.ID, color = Event),
      shape = 8, size = 0.7
    ) +
    geom_point(
      data = subset(plotdata, Event == "Death"),
      aes(x = tEvent, y = Patient.ID, color = Event),
      shape = 18
    ) +
    geom_point(
      data = subset(plotdata, Event %in% c("PriorLines1","PriorLines>1","PriorPARPi")),
      aes(x = tEvent, y = Patient.ID, color = Event),
      shape = 15, size = 1
    ) +
    ylab("Patient") +
    theme(axis.text.y = element_text(size = 4))+
    xlab("Time in months") +
    scale_fill_manual(values=c("EUDARIO Chemotherapy"=alpha(palette[7], 1), "EUDARIO Maintenance" =alpha(palette[8], 0.8)),name = "Therapy")+
    scale_color_manual(
      values = c(
        "Whole blood" = palette[1], "Plasma" = palette[2], "Progression" = palette[4], "Death" = "black", "PriorLines>1" = palette[5], "PriorLines1"= palette[6],  "PriorPARPi" = palette[3]
      ),
      name = "Event",
      breaks = c("Whole blood", "Plasma", "Progression", "Death" ,"PriorLines>1", "PriorLines1", "PriorPARPi"),
      labels = c("Whole blood", "Plasma", "Progression", "Death" ,"PriorLines>1", "PriorLines1", "PriorPARPi")
    ) +
    guides(
      color = guide_legend(title = "Prior Therapy", override.aes = list(shape = 16))
    ) -> p.timeline
  
  png(paste0(current_patient,"_timeline.png"),
      width=10,
      height=6,
      units="in",
      res=500,
      type="cairo")
  ### plot image to file
  print(p.timeline)
  ### close file again
  dev.off()
}
