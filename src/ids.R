library(base)
library(readxl)


##Patient ID table that identifies Sample IDs with Patient ID and timepoints
ids <- read_excel("data/external/Sample Registry.xlsx", 
                         sheet = "Sample IDs") %>%
  mutate(visit_material = paste(Visite,Material,sep="_"))

intervals <- read_excel("data/external/Sample Registry.xlsx", 
                        sheet = "Intervals")
left_join(ids,intervals,by="Int.Patient.ID")->ids

rm(intervals)
