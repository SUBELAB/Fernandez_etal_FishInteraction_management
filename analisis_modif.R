library(car) 
library(tidyverse)
library(MASS)
library(tweedie)
library(statmod)
library(cowplot)
library(qqplotr)
library(broom)
library(readr)
library(lme4)
library(MuMIn)
library(gridExtra)
library(emmeans)
library(glmmTMB)
library(DHARMa)
##Titulo:Evaluación del efecto del régimen de manejo sobre las interacciones ecológicas de peces en bosques dominados por macroalgas pardas en Chile
##Autores: Catalina Castillo-Sáez1, Catalina S. Ruz1, Ítalo Fernández-Cisternas1, Alejandro Pérez-Matus1,2, Mauricio F. Landaeta2,3,4. 
##fecha: Diciembre 2024

#INFO UVCs------
Values<-read.csv("UVC_values.csv") 
Values$date <- as.Date(Values$date)
Values$month <- as.numeric(format(Values$date, "%m"))
Values$year <- as.numeric(format(Values$date, "%y"))
Values <- Values %>%
  mutate(
    date = ymd(date),
    season = case_when(
      (month(date) == 12 & day(date) >= 21) | (month(date) %in% c(1, 2)) | (month(date) == 3 & day(date) <= 20) ~ "summer",
      (month(date) == 3 & day(date) >= 21) | (month(date) %in% c(4, 5)) | (month(date) == 6 & day(date) <= 20) ~ "fall",
      (month(date) == 6 & day(date) >= 21) | (month(date) %in% c(7, 8)) | (month(date) == 9 & day(date) <= 20) ~ "winter",
      (month(date) == 9 & day(date) >= 21) | (month(date) %in% c(10, 11)) | (month(date) == 12 & day(date) <= 20) ~ "spring"
    )
  )

Values2 <- Values %>% 
  group_by(site, management,date, year, season, method, transect) %>% 
  summarise(
    richness=sum(rich_std),
    abundance=sum(abun_std),
    biomass=sum(biom_std)) 
Values2$siteMan<- paste(Values2$site,Values2$management)
Values2$time_rf <- paste(Values2$year, Values2$season)

Values_summary <- Values2 %>% 
  group_by(site, management) %>% 
  summarise(total_richness = mean(richness, na.rm = TRUE),
            se_richness = sd(richness)/sqrt(n()),
            total_abundance= mean(abundance),
            se_abundance= sd(abundance)/sqrt(n()),
            total_biomass= mean(biomass),
            se_biomass= sd(biomass)/sqrt(n()))

Values_summary$management <-factor(Values_summary$management, 
                                   levels = c("MPA", "TURF", "OA")) 

table(Values2$site, Values2$management)
table(Values2$site)
table(Values2$time_rf,Values2$site)

#1. Richness----
ggplot(Values_summary, aes(x = management, fill = management)) +
  geom_bar(aes(y = total_richness), stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = total_richness - se_richness, 
                    ymax = total_richness + se_richness), 
                position = position_dodge(width = 0.9), width = 0.2) +
  theme_bw() +
  labs(title = "UVC-registered Richness",
       x = "Management regime",
       y = "Richness (Species/m²)",
       fill = "Management regime") +
  scale_fill_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1")) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 16),
        panel.grid = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Expande el eje y automáticamente
  facet_grid(~site)

#ggsave("uvc_richness_years.png",width = 10, height = 7)

#Analisis richness
hist(Values2$richness)
richnUVC_gauss_mix <- glmmTMB(richness ~ site + management + (1 | time_rf) ,
                              family = gaussian(),
                              data = Values2)

summary(richnUVC_gauss_mix)
Anova(richnUVC_gauss_mix, type = "III")
res <- simulateResiduals(richnUVC_gauss_mix)
plot(res)

richnUVC_gauss_mix <- glmmTMB(richness ~ site + management + (1 | time_rf) ,
                              family = gaussian(),
                              dispformula = ~ management,
                              data = Values2)

summary(richnUVC_gauss_mix)
Anova(richnUVC_gauss_mix, type = "III")
res <- simulateResiduals(richnUVC_gauss_mix)
plot(res)

#2. Abundance------
ggplot(Values_summary, aes(x = management, fill = management)) +
  geom_bar(aes(y = total_abundance), stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = total_abundance - se_abundance, 
                    ymax = total_abundance + se_abundance), 
                position = position_dodge(width = 0.9), width = 0.2) +
  theme_bw() +
  labs(title = "UVC-registered Abundance",
       x = "Management regime",
       y = "Abundance (Ind/m²)",
       fill = "Management regime") +
  scale_fill_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1")) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 16),
        panel.grid = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Expande el eje y automáticamente
  facet_grid(~site)

#ggsave("uvc_abundance_year.png",width = 10, height = 7)

##Analisis abundance
hist(Values2$abundance)
abunUVC_gauss<- glmmTMB(abundance ~ site+ management + (1 | time_rf),
                        family = Gamma(link = "log"),
                        data = Values2)

summary(abunUVC_gauss)
Anova(abunUVC_gauss)
res <- simulateResiduals(abunUVC_gauss)
plot(res)

#3. Biomass------
ggplot(Values_summary, aes(x = management, fill = management)) +
  geom_bar(aes(y = total_biomass), stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = total_biomass - se_biomass, 
                    ymax = total_biomass + se_biomass), 
                position = position_dodge(width = 0.9), width = 0.2) +
  theme_bw() +
  labs(title = "UVC-registered Biomass",
       x = "Management regime",
       y = "Biomass (Kg/m²)",
       fill = "Management regime") +
  scale_fill_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1")) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 16),
        panel.grid = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Expande el eje y automáticamente
  facet_grid(~site)

#ggsave("uvc_biomass.png",width = 10, height = 7)

#Analisis biomass
hist(Values2$biomass)
BiomUVC_mix <- glmmTMB(biomass ~ site+ management + (1 | time_rf),
                              family = Gamma(link = "log"),
                              data = Values2)

summary(BiomUVC_mix)
Anova(BiomUVC_mix, type="III")
res <- simulateResiduals(BiomUVC_mix)
plot(res)

#varianza
plotResiduals(res, Values2$management)
plotResiduals(res, Values2$site)
plotResiduals(res, Values2$siteMan)
plotResiduals(res, Values2$season)

Biom_UVC_mix <- glmmTMB(biomass ~ site + management + (1 | time_rf),
                        family = Gamma(link = "log"),
                        dispformula = ~ season ,
                        data = Values2)

summary(Biom_UVC_mix)
Anova(Biom_UVC_mix, type="III")
res <- simulateResiduals(Biom_UVC_mix)
plot(res)

# Comparaciones de Tukey # ComparacionessetClassUnion() de Tukey 
tukey_biom <- emmeans(Biom_UVC_mix, pairwise ~ site * management  , adjust = "tukey")
print(tukey_biom)

#INFO RUVS---------
#write.csv(Variables,"Alldata_final.csv")
Variables<-read.csv("Alldata_final.csv") %>% 
  group_by(site, management, video_id) %>% 
  summarise(rich= sum(richness),
            rich_std=sum(richness_std),
            abun=sum(abundance),
            abun_std=sum(abundance_std),
            biom=sum(biomass),
            biom_std=sum(biomass_std),
            int_str=sum(interaction_strength),
            detec=sum(detections),
            time=first(useful_time)) 

random_var <- read.csv("Alldata.csv") %>% 
  dplyr::select(year, season, video_id, date) %>%
  distinct(video_id, .keep_all = TRUE)  # Mantiene una sola fila por video_id

Variables <- left_join(Variables, random_var, by = "video_id", "year")
Variables$time_rf=paste(Variables$year,Variables$season)
Variables$siteMan= paste(Variables$site,Variables$management)

Variables_summary <- Variables %>% 
  group_by(management,site) %>%
  summarise(total_abundance = mean(abun_std, na.rm = TRUE),
            se_abundance = sd(abun_std)/sqrt(length(management)),
            total_biomass = mean(biom_std, na.rm = TRUE),
            se_biomass = sd(biom_std)/sqrt(length(management)),
            total_richness = mean(rich_std, na.rm = TRUE),
            se_richness = sd(rich_std)/sqrt(length(management)),
            total_int=mean(int_str),
            se_int = sd(int_str)/sqrt(length(management)))
Variables_summary$management<- factor(levels = c("MPA", "TURF", "OA"),Variables_summary$management)

#4. Interaction strength------------------
ggplot(Variables_summary, aes(x = management, fill = management)) +
  geom_bar(aes(y = total_int), stat = "identity", position = "dodge") +
  # Barras de error
  geom_errorbar(aes(ymin = total_int - se_int, 
                    ymax = total_int + se_int), 
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_bw() +
  labs(title = "Fish Interaction Strength",
       x = "Management regime",
       y = "Strength (n° interactions/2m²*hr)",
       fill = "Management regime") +
  facet_grid( ~site)+
  scale_fill_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1")) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 16),
        panel.grid = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

#ggsave("fish_strength.png",width = 10, height = 7)

#analisis fuerza interaccion
hist(Variables$detec)
Finteract_mixed <- glmmTMB(
  detec ~ site + management + offset(log(time)) + (1 | time_rf),
  family = nbinom2(link = "log"),
  data = Variables
)

summary(Finteract_mixed)
Anova(Finteract_mixed, type="III")
res <- simulateResiduals(Finteract_mixed)
plot(res)

# Comparaciones de Tukey 
tukey_strength <- emmeans(Finteract_mixed, pairwise ~ site + management  , adjust = "tukey")
print(tukey_strength)

#CORRELACIONES UVC-strength------
#promediar un valor por dia
Values_days <- Values2 %>% 
  group_by(site, management,date) %>% 
  summarise(rich_UVC = mean(richness, na.rm = TRUE),
            se_richness = sd(richness)/sqrt(length(management)),
            abun_UVC= mean(abundance),
            se_abundance= sd(abundance)/sqrt(length(management)),
            biom_UVC= mean(biomass),
            se_biomass= sd(biomass)/sqrt(length(management)))

#convertir a m²
Variables_days <- Variables %>% 
  mutate(rich_std = rich_std/2,
         abun_std= abun_std/2,
         biom_std= biom_std/2,
         int_str= int_str) %>% 
  group_by(site, management,date) %>% 
  summarise(rich_RUV = mean(rich_std, na.rm = TRUE),
            abun_RUV= mean(abun_std),
            biom_RUV= mean(biom_std),
            int_str= mean(int_str)) %>% 
  mutate(date=as.Date(date))

UVC_RUVS <- Variables_days %>% 
  inner_join(Values_days, by = c("site", "management", "date")) 

#1. Richness----
# Agregamos coordenadas para ubicar el texto en cada panel
cor_rich_labelsUVC <- UVC_RUVS %>%
  group_by(management) %>%
  summarise(
    R = cor(rich_UVC, int_str, use = "complete.obs"),
    x = max(rich_UVC, na.rm = TRUE) * 0.7,
    y = max(int_str, na.rm = TRUE) * 0.9
  )

richness_plot2 <- ggplot(UVC_RUVS, aes(x = rich_UVC, y = int_str, color = management, fill = management)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(title = "UVC Richness",
       x = "Richness (species/m²)",
       y = "Interaction strength (int/2m²*hr)") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 16),
        panel.grid = element_blank()) + 
  geom_text(data = cor_rich_labelsUVC,
            aes(x = x, y = y, label = paste0("R = ", round(R, 3))),
            inherit.aes = FALSE, size = 5, hjust = 0) +
  scale_color_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1")) +
  scale_fill_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1"))

#2. Abundance----
cor_abun_labelsUVC <- UVC_RUVS %>%
  group_by(management) %>%
  summarise(
    R = cor(abun_UVC, int_str, use = "complete.obs"),
    x = max(abun_UVC, na.rm = TRUE) * 0.7,
    y = max(int_str, na.rm = TRUE) * 0.9
  )

abundance_plot2 <- ggplot(UVC_RUVS, aes(x = abun_UVC, y = int_str, color = management, fill = management)) +
  geom_point() +  
  geom_smooth(method = "lm", se = F) +  # Línea de regresión
  labs(title = "UVC Abundance",
       x = "Abundance (ind/m²)",
       y = "Interaction strength (int/2m²*hr)") +
  theme_bw() +  
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 16),
        panel.grid = element_blank()) +
  geom_text(data = cor_abun_labelsUVC,
            aes(x = x, y = y, label = paste0("R = ", round(R, 3))),
            inherit.aes = FALSE, size = 5, hjust = 0) +
  scale_color_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1")) +
  scale_fill_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1")) 

#3. Biomass----

cor_biom_labelsUVC <- UVC_RUVS %>%
  group_by(management) %>%
  summarise(
    R = cor(biom_UVC, int_str, use = "complete.obs"),
    x = max(biom_UVC, na.rm = TRUE) * 0.7,
    y = max(int_str, na.rm = TRUE) * 0.9
  )

biomass_plot2 <- ggplot(UVC_RUVS, aes(x = biom_UVC, y = int_str, color = management, fill = management)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(title = "UVC Biomass",
       x = "Biomass (kg/m²)",
       y = "Interaction strength (int/2m²*hr)") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 16),
        panel.grid = element_blank()) +
  geom_text(data = cor_biom_labelsUVC,
            aes(x = x, y = y, label = paste0("R = ", round(R, 3))),
            inherit.aes = FALSE, size = 5, hjust = 0) +
  scale_color_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1")) +
  scale_fill_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1"))

grid.arrange(richness_plot2, abundance_plot2, biomass_plot2, ncol = 3)
