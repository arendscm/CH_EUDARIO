library(base)
library(dplyr)
library(forcats)
library(ggplot2)
library(ggsci)

source("src/create_comut.R")
source("src/global_functions_themes.R")

##load filtered mutation data
load('data/interim/seqdata_filtered.RData')

sub <- df.filtered.c1d1 %>% 
  mutate(Gene = ifelse(Gene=="U2AF1;U2AF1L5","U2AF1",Gene))%>%
  filter(TVAF > 0.01) %>% 
  filter(tag=="true") %>%
  group_by(Patient.ID) %>%
  mutate(nom = n())%>% 
  filter(nom > 1) %>%
  data.frame %>% 
  group_by(Patient.ID,Gene)%>%
  mutate(maxVAF = max(TVAF))%>%
  data.frame%>%
  dplyr::select(Sample,Gene,maxVAF)%>%
  mutate(TVAF = maxVAF) %>%
  unique

sub$Gene %>% table %>% data.frame -> order 
names(order)<- c("Gene","order")

createComut(sub,dis=0.01)%>%
  full_join(.,order,by=c("Gene_1"="Gene"))%>%
  full_join(.,order,by=c("Gene_2"="Gene"))%>%
  ggplot(., aes(y=forcats::fct_rev(reorder(Gene_1,-order.x)), x = reorder(Gene_2,-order.y))) +
  geom_point(aes(size = number_1_and_2, color = fraction_1_then_2), shape = 15, stat = "identity") +
  #geom_point(aes(size = number_1_and_2), shape = 0, color = "#374E55FF") +
  scale_color_gradient2(high =  "#4DBBD5FF", low = "#E64B35FF", mid = "grey",
                        midpoint = 50,
                        name= "Fraction of \nunambiguous\ncases\n(Gene 1 > Gene2)") +
  scale_size(range = c(1,10), name= "Number of\nco-occurences", breaks = c(2,5,10,20)) +
  #theme_cowplot() +
  my_theme()+
  theme(axis.text.x = element_text(angle = 90, vjust = .5,
                                   size = 10, hjust = 0),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "#f0f0f0"),
        axis.text.y = element_text(size = 10))+
  theme(axis.text.x.top = element_text(vjust = 0.5)) +
  theme(legend.key.height=unit(0.5,"cm"))+
  guides(    color = guide_colorbar(order = 1),
             fill = guide_legend(order = 0)) +
  scale_x_discrete(position = "top") +
  xlab(label= "Gene 2") +
  ylab(label="Gene 1") +
  labs(title = NULL) -> p.comut

png("ouput/figures/comut_cvd.png",width=7, height=5,units="in",res=500,type="cairo")
p.comut
dev.off()

###compare comutation burden of PPM1D mutated patients with others CH positive patients
df.filtered.c1d1 %>% 
  mutate(Gene = ifelse(Gene=="U2AF1;U2AF1L5","U2AF1",Gene))%>%
  filter(TVAF > 0.01) %>% 
  filter(tag=="true") %>%
  group_by(Patient.ID) %>%
  mutate(nom = n())%>% 
  data.frame %>%
  mutate(ppm1d = is.element(Patient.ID, df.filtered.c1d1 %>%  
                              filter(TVAF > 0.01) %>% 
                              filter(tag=="true") %>% 
                              filter(Gene=="PPM1D") %>% .$Patient.ID)) %>%
  dplyr::select(Patient.ID,ppm1d,nom)%>%
  unique%>%
  ggplot(aes(x=ppm1d,group=factor(nom),fill=factor(nom)))+
  geom_bar(stat="count",position="stack")+
  scale_fill_manual(values=scales::seq_gradient_pal(high = "#E64B35FF", low = "#4DBBD5FF",space="Lab")(seq(0,1,length.out=8)))+
  my_theme()
