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
library(ggdist)




##Titulo:Evaluación del efecto del régimen de manejo sobre las interacciones ecológicas de peces en bosques dominados por macroalgas pardas en Chile
##Autores: Catalina Castillo-Sáez1, Catalina S. Ruz1, Ítalo Fernández-Cisternas1, Alejandro Pérez-Matus1,2, Mauricio F. Landaeta2,3,4. 
##fecha: Diciembre 2024

#write.csv(Variables,"Alldata_final.csv")

#INFO RUVS----------------------------

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
VariablesCHA <- Variables %>% filter(site == "Chañaral de Aceituno")
VariablesLC <- Variables %>% filter(site == "Las Cruces")

Variables_summary <- Variables %>% 
  group_by(management, site,year) %>%
  summarise(total_abundance = mean(abun_std, na.rm = TRUE),
            se_abundance = sd(abun_std)/sqrt(length(management)),
            total_biomass = mean(biom_std, na.rm = TRUE),
            se_biomass = sd(biom_std)/sqrt(length(management)),
            total_richness = mean(rich_std, na.rm = TRUE),
            se_richness = sd(rich_std)/sqrt(length(management)),
            total_int=mean(int_str),
            se_int = sd(int_str)/sqrt(length(management)))

Variables_summary$management<- factor(levels = c("MPA", "TURF", "OA"),Variables_summary$management)

#1. Richness-------
##plot richness
ggplot(Variables_summary, aes(x = management, fill = management)) +
  geom_bar(aes(y = total_richness), stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = total_richness - se_richness, 
                    ymax = total_richness + se_richness), 
                position = position_dodge(width = 0.9), width = 0.2) +
  theme_bw() +
  labs(title = "RUV-registered Richness",
       x = "Management regime",
       y = "Richness (Species/2m²*hr)",
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
  facet_grid(~year ~site)

ggsave("fish_richness_year.png",width = 10, height = 7)

##Analisis general
richness_normal_mix <- lmer(
  rich ~ site + management + offset(log(time)) + (1 | year) ,
  data = Variables  # Por defecto para gaussian
)
summary(richness_normal_mix)
Anova(richness_normal_mix, type="III")

#Ajustar modelo a Chañaral
hist(VariablesCHA$rich)
richness_normal_mix1 <- lmer(rich ~ management + offset(log(time)) + (1 | year) ,
                             data = VariablesCHA)
summary(richness_normal_mix1)
Anova(richness_normal_mix1, type="III")

r.squaredGLMM(richness_normal_mix1)*100

#Ajustar modelo a Las Cruces
hist(VariablesLC$rich)
richness_normal_mix2 <- lmer(rich ~ management + offset(log(time)) + (1 | year) ,
                             data = VariablesLC,
)
summary(richness_normal_mix2)
Anova(richness_normal_mix2, type="III")                            
                             
r.squaredGLMM(richness_normal_mix2)*100

# Analisis de residuales del modelo 
# Calcular residuales y valores ajustados
VariablesLC$residuals <- residuals(richness_normal_mix2, type = "pearson")  # Residuales tipo Pearson
VariablesLC$fitted <- fitted(richness_normal_mix2)  # Valores ajustados
# guarda la informacion requerida para el analisis de residuales

# Q-Q plot
qq.model <- ggplot(data = VariablesLC, aes(sample = residuals)) + 
  stat_qq_point() +
  stat_qq_line() +
  stat_qq_band(alpha = 0.3) +
  theme_bw() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16)
  )

# Residuales vs valores predichos
res.fit <- ggplot(data = VariablesLC, aes(x = fitted, y = residuals)) + 
  geom_point(aes(col = management), size = 2) +
  geom_hline(yintercept = 0, linetype = 2, size = 1.2) +
  theme_bw() +
  labs(x = "Fitted", y = "Residuals") +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16),
    legend.position = "none"
  )

plot_grid(qq.model,res.fit, nrow=1)

#test homocedasticidad
leveneTest(VariablesLC$residuals ~ VariablesLC$management * VariablesLC$site)##son homocedasticos


#2. Abundance -----------------------------------------------------

#Plot 
ggplot(Variables_summary, aes(x = management,fill = management)) +
  geom_bar(aes(y = total_abundance), stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = total_abundance - se_abundance, 
                    ymax = total_abundance + se_abundance), 
                position = position_dodge(width = 0.9), width = 0.2) +
  theme_bw() +
  labs(title = "RUV-registered Abundance",
       x = "Management regime",
       y = "Abundance (MaxN/2m²*hr)",
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
  facet_grid(~year ~site)

ggsave("fish_abundance_year.png",width = 10, height = 7)

#Distribucion abundance
hist(VariablesCHA$abun)

# Modelo mixto general con distribución binomial negativa
abundance_nb_mixed <- glmer.nb(
  abun ~ site + management + offset(log(time)) + (1 | year),
  data = Variables
)

# Ajustar el modelo mixto a Chañaral
abundance_nb_mixed1 <- glmer.nb(
  abun ~ management + offset(log(time)) + (1 | year),
  data = VariablesCHA)

summary(abundance_nb_mixed1)
Anova(abundance_nb_mixed1)
r.squaredGLMM(abundance_nb_mixed1)*100

# Ajustar el modelo mixto a Las Cruces
abundance_nb_mixed2 <- glmer.nb(
  abun ~ management + offset(log(time)) + (1 | year),
  data = VariablesLC)

summary(abundance_nb_mixed2)
Anova(abundance_nb_mixed2)
r.squaredGLMM(abundance_nb_mixed2)*100

# Comparaciones de Tukey para 'management'
tukey_abun <- emmeans(abundance_nb_mixed1, pairwise ~ management, adjust = "tukey")
print(tukey_abun)

# Analisis de residuales del modelo 
# Calcular residuales y valores ajustados
Variables$residuals <- residuals(abundance_nb_mixed, type = "pearson")  # Residuales tipo Pearson
Variables$fitted <- fitted(abundance_nb_mixed)  # Valores ajustados
# guarda la informacion requerida para el analisis de residuales

# Q-Q plot
qq.model <- ggplot(data = Variables, aes(sample = residuals)) + 
  stat_qq_point() +
  stat_qq_line() +
  stat_qq_band(alpha = 0.3) +
  theme_bw() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16)
  )

# Residuales vs valores predichos
res.fit <- ggplot(data = Variables, aes(x = fitted, y = residuals)) + 
  geom_point(aes(col = management), size = 2) +
  geom_hline(yintercept = 0, linetype = 2, size = 1.2) +
  theme_bw() +
  labs(x = "Fitted", y = "Residuals") +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16),
    legend.position = "none"
  )

plot_grid(qq.model,res.fit, nrow=1)


#3. Biomass -----------------------------------------
#plot biomass
ggplot(Variables_summary, aes(x = management, fill= management)) +
  geom_bar(aes(y = total_biomass), stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = total_biomass - se_biomass, 
                    ymax = total_biomass + se_biomass), 
                position = position_dodge(width = 0.9), width = 0.2) +
  theme_bw() +
  labs(title = "RUV-registered Biomass",
       x = "Management regime",
       y = "Biomass (Kg/2m²*hr)",
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
  facet_grid(~year ~site)

ggsave("fish_biomass_year.png",width = 10, height = 7)

hist(Variables$biom_std)
hist(Variables$biom)

##lognormal
# Transformar la variable de respuesta
Variables$log_biom <- log(Variables$biom)

# Ajustar el modelo lineal mixto con distribución log-normal a Chañaral
biomass_lognormal_mixto1 <- lmer(
  log_biom ~ management + offset(log(time)) + (1 | year),
  data = VariablesCHA)

summary(biomass_lognormal_mixto1)
Anova(biomass_lognormal_mixto1, type="III")
r.squaredGLMM(biomass_lognormal_mixto1)*100

# Ajustar el modelo lineal mixto con distribución log-normal a Las Cruces
biomass_lognormal_mixto2 <- lmer(
  log_biom ~ management + offset(log(time)) + (1 | year),
  data = VariablesLC)

summary(biomass_lognormal_mixto2)
Anova(biomass_lognormal_mixto2, type="III")
r.squaredGLMM(biomass_lognormal_mixto2)*100

# Analisis de residuales del modelo 
# Calcular residuales y valores ajustados
Variables$residuals <- residuals(biomass_lognormal_mixto, type = "pearson")  # Residuales tipo Pearson
Variables$fitted <- fitted(biomass_lognormal_mixto)  # Valores ajustados
# guarda la informacion requerida para el analisis de residuales

# Q-Q plot
qq.model <- ggplot(data = Variables, aes(sample = residuals)) + 
  stat_qq_point() +
  stat_qq_line() +
  stat_qq_band(alpha = 0.3) +
  theme_bw() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16)
  )

# Residuales vs valores predichos
res.fit <- ggplot(data = Variables, aes(x = fitted, y = residuals)) + 
  geom_point(aes(col = management), size = 2) +
  geom_hline(yintercept = 0, linetype = 2, size = 1.2) +
  theme_bw() +
  labs(x = "Fitted", y = "Residuals") +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16),
    legend.position = "none"
  )

plot_grid(qq.model,res.fit, nrow=1)

#4. Fuerza de interaccion------------------

#plot
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
  facet_grid(~year ~site)+
  scale_fill_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1")) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 16),
        panel.grid = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("fish_strength_year.png",width = 10, height = 7)

#Modelo general
Finteract_nb_mixed <- glmer.nb(detec ~ site + management + offset(log(time)) + (1 | year),
                                data = Variables)
summary(Finteract_nb_mixed)
Anova(Finteract_nb_mixed, type="III")

# Modelo mixto con distribución binomial negativa para Chañaral
Finteract_nb_mixed1 <- glmer.nb(detec ~ management + offset(log(time)) + (1 | year),
                               data = VariablesCHA)
                               
summary(Finteract_nb_mixed1)
Anova(Finteract_nb_mixed1, type="III")
r.squaredGLMM(Finteract_nb_mixed1)*100

# Modelo mixto con distribución binomial negativa para Las Cruces
Finteract_nb_mixed2 <- glmer.nb(detec ~ management + offset(log(time)) + (1 | year),
                                data = Variables %>% filter(site == "Las Cruces"))

summary(Finteract_nb_mixed2)
Anova(Finteract_nb_mixed2, type="III")
r.squaredGLMM(Finteract_nb_mixed2)*100

# Analisis de residuales del modelo 
# Calcular residuales y valores ajustados
Variables$residuals <- residuals(Finteract_nb_mixed, type = "pearson")  # Residuales tipo Pearson
Variables$fitted <- fitted(Finteract_nb_mixed)  # Valores ajustados
# guarda la informacion requerida para el analisis de residuales

# Q-Q plot
qq.model <- ggplot(data = Variables, aes(sample = residuals)) + 
  stat_qq_point() +
  stat_qq_line() +
  stat_qq_band(alpha = 0.3) +
  theme_bw() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16)
  )

# Residuales vs valores predichos
res.fit <- ggplot(data = Variables, aes(x = fitted, y = residuals)) + 
  geom_point(aes(col = management), size = 2) +
  geom_hline(yintercept = 0, linetype = 2, size = 1.2) +
  theme_bw() +
  labs(x = "Fitted", y = "Residuals") +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16),
    legend.position = "none"
  )

plot_grid(qq.model,res.fit, nrow=1)

#CORRELACIONES RUV-strength-----
Variables$management<- factor(levels = c("MPA", "TURF", "OA"),Variables$management)

##1. Richness----
cor_rich_labels <- Variables %>%
  group_by(management) %>%
  summarise(R = cor(rich_std, int_str, use = "complete.obs"))
# Agregamos coordenadas para ubicar el texto en cada panel
cor_labels <- Variables %>%
  group_by(management) %>%
  summarise(
    R = cor(rich_std, int_str, use = "complete.obs"),
    x = max(rich_std, na.rm = TRUE) * 0.7,
    y = max(int_str, na.rm = TRUE) * 0.9
  )

richness_plot1 <- ggplot(Variables, aes(x = rich_std, y = int_str, color = management, fill = management)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "RUVs richness",
       x = "Richness (species/2m²*hr)",
       y = "Interaction strength (int/2m²*hr)") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 16),
        panel.grid = element_blank()) +   facet_grid(~management) +
  geom_text(data = cor_labels,
            aes(x = x, y = y, label = paste0("R = ", round(R, 3))),
            inherit.aes = FALSE, size = 5, hjust = 0) +
  scale_color_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1")) +
  scale_fill_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1"))

##2. Abundance----
cor_abun_labels <- Variables %>%
  group_by(management) %>%
  summarise(
    R = cor(abun_std, int_str, use = "complete.obs"),
    x = max(abun_std, na.rm = TRUE) * 0.7,
    y = max(int_str, na.rm = TRUE) * 0.9
  )

abundance_plot1 <- ggplot(Variables, aes(x = abun_std, y = int_str, color = management, fill = management)) +
  geom_point() +  
  geom_smooth(method = "lm", se = TRUE) +  # Línea de regresión
  labs(title = "RUVs abundance",
       x = "Abundance (maxN/2m²*hr)",
       y = "Interaction strength (int/2m²*hr)") +
  theme_bw() +  
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 16),
        panel.grid = element_blank()) +
  facet_grid(~management) +
  geom_text(data = cor_abun_labels,
            aes(x = x, y = y, label = paste0("R = ", round(R, 3))),
            inherit.aes = FALSE, size = 5, hjust = 0) +
  scale_color_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1")) +
  scale_fill_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1"))
##3. Biomass----
cor_biom_labels <- Variables %>%
  group_by(management) %>%
  summarise(
    R = cor(biom_std, int_str, use = "complete.obs"),
    x = max(biom_std, na.rm = TRUE) * 0.7,
    y = max(int_str, na.rm = TRUE) * 0.9
  )

biomass_plot1 <- ggplot(Variables, aes(x = biom_std, y = int_str, color = management, fill = management)) +
  geom_point() +  
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "RUVs biomass",
    x = "Biomass (kg/2m²*hr)",
    y = "Interaction strength (int/2m²*hr)"
  ) +
  theme_bw() +  
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 16),
        panel.grid = element_blank()) +
  facet_grid(~management) +
  geom_text(data = cor_biom_labels,
            aes(x = x, y = y, label = paste0("R = ", round(R, 3))),
            inherit.aes = FALSE, size = 5, hjust = 0) +
  scale_color_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1")) +
  scale_fill_manual(values = c("MPA" = "springgreen4", "TURF" = "dodgerblue", "OA" = "brown1"))

grid.arrange(richness_plot1, abundance_plot1, biomass_plot1, ncol = 1)


#CORRELACIONES RUV-UVC-----
Values<-read.csv("UVC_values.csv")

#promediar un valor por dia
Values_days <- Values %>% 
  group_by(site, management,date, transect) %>%
  summarise(rich_std = sum(rich_std),
            abun_std = sum(abun_std),
            biom_std = sum(biom_std)) %>% 
  group_by(site, management,date) %>% 
  summarise(rich_UVC = mean(rich_std, na.rm = TRUE),
            se_richness = sd(rich_std)/sqrt(length(management)),
            abun_UVC= mean(abun_std),
            se_abundance= sd(abun_std)/sqrt(length(management)),
            biom_UVC= mean(biom_std),
            se_biomass= sd(biom_std)/sqrt(length(management)))
  
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
            int_str= mean(int_str))

UVC_RUVS <- Variables_days %>% 
  inner_join(Values_days, by = c("site", "management", "date")) 

#1. Richness------
richness_plot <- ggplot(UVC_RUVS, aes(x = rich_RUV, y = rich_UVC)) +
  geom_point() +  # Agregar puntos
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # Línea de regresión
  labs(title = "Richness",
       x = "RUVs Richness (species/m²*hr)",
       y = "UVCs Richness (species/m²)") +
  annotate("text", x = max(UVC_RUVS$rich_RUV) * 0.7, 
           y = max(UVC_RUVS$rich_UVC) * 0.9, 
           label = paste("R =", round(0.365, 3)), 
           size = 5, hjust = 0) +
  theme_linedraw() + # Aplicar un tema minimalista
  theme(plot.title = element_text(hjust = 0.5))  # Centrar el título

richness_correlation <- cor.test(UVC_RUVS$rich_RUV, UVC_RUVS$rich_UVC, use = "complete.obs")

#2. Abundance------
abundance_plot <- ggplot(UVC_RUVS, aes(x = abun_RUV, y = abun_UVC)) +
  geom_point() +  # Agregar puntos
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # Línea de regresión
  labs(title = "Abundance",
       x = "RUVs Abundance (maxN/m²*hr)",
       y = "UVCs Abundance (ind/m²)") +
  annotate("text", x = max(UVC_RUVS$abun_RUV) * 0.7, 
           y = max(UVC_RUVS$abun_UVC) * 0.9, 
           label = paste("R =", round(0.701, 3)), 
           size = 5, hjust = 0) +
  theme_linedraw() + # Aplicar un tema minimalista
  theme(plot.title = element_text(hjust = 0.5))  # Centrar el título

abundance_correlation <- cor.test(UVC_RUVS$abun_RUV, UVC_RUVS$abun_UVC, use = "complete.obs")

#3. Biomass------
biomass_plot <-ggplot(UVC_RUVS, aes(x = biom_RUV, y = biom_UVC)) +
  geom_point() +  # Agregar puntos
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # Línea de regresión
  labs(title = "Biomass",
       x = "RUVs Biomass (kg/m²*hr)",
       y = "UVCs Biomass (kg/m²)") +
  annotate("text", x = max(UVC_RUVS$biom_RUV) * 0.7, 
           y = max(UVC_RUVS$biom_UVC) * 0.9, 
           label = paste("R =", round(0.652, 3)), 
           size = 5, hjust = 0) +
  theme_linedraw() + # Aplicar un tema minimalista
  theme(plot.title = element_text(hjust = 0.5))  # Centrar el título

biomass_correlation <- cor.test(UVC_RUVS$biom_RUV, UVC_RUVS$biom_UVC, use = "complete.obs")

grid.arrange(richness_plot, abundance_plot, biomass_plot, ncol = 3)

#INFO UVCs------
Values_summary <- Values%>% 
  group_by(site, management, transect, date, year) %>%
  summarise(rich_std = sum(rich_std),
            abun_std = sum(abun_std),
            biom_std = sum(biom_std)) %>% 
  group_by(site, management, year) %>% 
  summarise(total_richness = mean(rich_std, na.rm = TRUE),
            se_richness = sd(rich_std)/sqrt(length(management)),
            total_abundance= mean(abun_std),
            se_abundance= sd(abun_std)/sqrt(n()),
            total_biomass= mean(biom_std),
            se_biomass= sd(biom_std)/sqrt(n()))

ValuesCHA <- Values %>% filter(site == "Chañaral de Aceituno")
ValuesLC <- Values %>% filter(site == "Las Cruces")


Values_summary$management <-factor(Values_summary$management, 
                                   levels = c("MPA", "TURF", "OA")) 

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
  facet_grid(~year~site)

ggsave("uvc_richness_years.png",width = 10, height = 7)

##Analisis 
hist(Values$richness)

# Ajustar el modelo mixto general
richnUVC_normal_mix <- lmer(
  richness ~ site + management + offset(log(meters)) + (1 | year) ,
  data = Values  # Por defecto para gaussian
)
summary(richnUVC_normal_mix)
Anova(richnUVC_normal_mix, type="III")

#Ajustar el modelo a Chañaral
richnUVC_normal_mix1 <- lmer(
  richness ~ management + offset(log(meters)) + (1 | year) ,
  data = ValuesCHA
)
summary(richnUVC_normal_mix1)
Anova(richnUVC_normal_mix1, type="III")

#Ajustar el modelo a Las Cruces

# Analisis de residuales del modelo 
# Calcular residuales y valores ajustados
Values$residuals <- residuals(richnUVC_normal_mix, type = "pearson")  # Residuales tipo Pearson
Values$fitted <- fitted(richnUVC_normal_mix)  # Valores ajustados
# guarda la informacion requerida para el analisis de residuales

# Q-Q plot
qq.model <- ggplot(data = Values, aes(sample = residuals)) + 
  stat_qq_point() +
  stat_qq_line() +
  stat_qq_band(alpha = 0.3) +
  theme_bw() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16)
  )

# Residuales vs valores predichos
res.fit <- ggplot(data = Values, aes(x = fitted, y = residuals)) + 
  geom_point(aes(col = management), size = 2) +
  geom_hline(yintercept = 0, linetype = 2, size = 1.2) +
  theme_bw() +
  labs(x = "Fitted", y = "Residuals") +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16),
    legend.position = "none"
  )

plot_grid(qq.model,res.fit, nrow=1)

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
  facet_grid(~year ~site)

ggsave("uvc_abundance_year.png",width = 10, height = 7)

#Analisis
#Distribucion abundance
hist(Values$abundance)

##Binomial negativa
abunUVC_nb <- glm.nb(abundance ~ site+ management + offset(log(meters)), 
                       data = Values #%>% filter(abun <50)
)
summary(abunUVC_nb)
Anova(abunUVC_nb)

##con factores aleatorios 
# Modelo mixto con distribución binomial negativa
abunUVC_nb_mixed <- glmer.nb(
  abundance ~ site + management + offset(log(meters)) + (1 | year),
  data = Values# %>% filter(abun <50)
)

summary(abunUVC_nb_mixed)
Anova(abunUVC_nb_mixed)
AIC(abunUVC_nb_mixed,abunUVC_nb)

# Analisis de residuales del modelo 
# Calcular residuales y valores ajustados
Values$residuals <- residuals(abunUVC_nb_mixed, type = "pearson")  # Residuales tipo Pearson
Values$fitted <- fitted(abunUVC_nb_mixed)  # Valores ajustados
# guarda la informacion requerida para el analisis de residuales

# Q-Q plot
qq.model <- ggplot(data = Values, aes(sample = residuals)) + 
  stat_qq_point() +
  stat_qq_line() +
  stat_qq_band(alpha = 0.3) +
  theme_bw() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16)
  )

# Residuales vs valores predichos
res.fit <- ggplot(data = Values, aes(x = fitted, y = residuals)) + 
  geom_point(aes(col = management), size = 2) +
  geom_hline(yintercept = 0, linetype = 2, size = 1.2) +
  theme_bw() +
  labs(x = "Fitted", y = "Residuals") +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16),
    legend.position = "none"
  )

plot_grid(qq.model,res.fit, nrow=1)

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
  facet_grid(~year ~site)

ggsave("uvc_biomass_year.png",width = 10, height = 7)

#Analisis
hist(Values$biomass)
#Gamma mixto
biomUVC_gamma_mixto<- glmer(
  biomass~ site + management + offset(log(meters)) + (1 | year),
  family = Gamma(link = "log"),
  data = Values
)

summary(biomUVC_gamma_mixto)
Anova(biomUVC_gamma_mixto,  type = "III")

hist(Values$log_biom)
##logvariables_management##lognormal mixto
# Transformar la variable de respuesta
Values$log_biom <- log(Values$biomass)

# Ajustar el modelo lineal mixto con distribución log-normal
biomUVC_lognormal_mixto <- lmer(
  log_biom ~ site + management + offset(log(meters)) + (1 | year),
  data = Values
)
summary(biomUVC_lognormal_mixto)
Anova(biomUVC_lognormal_mixto, type="III")

AIC(biomUVC_lognormal_mixto, biomUVC_gamma_mixto)

#CORRELACIONES UVC-strength------
#1. Richness----
richness_plot3 <-ggplot(UVC_RUVS, aes(x = rich_UVC, y = int_str)) +
  geom_point() +  # Agregar puntos
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # Línea de regresión
  labs(title = "UVCs richness",
       x = "Richness (species/m²)",
       y = "Interaction strength (interactions/2m²*hr)") +
  theme_bw()  +# Aplicar un tema minimalista
  theme(plot.title = element_text(hjust = 0.5)) +  # Centrar el título
  annotate("text", x = max(UVC_RUVS$rich_UVC) * 0.7, y = max(UVC_RUVS$int_str) * 0.9, 
           label = paste0("R = ", round(strength_richness_cor3$estimate, 3)), 
           size = 5, hjust = 0) # Centrar el título # Centrar el título  # Centrar el título

strength_richness_cor3 <- cor.test(UVC_RUVS$rich_UVC, UVC_RUVS$int_str, use = "complete.obs")

#2. Abundance----
abundance_plot3 <- ggplot(UVC_RUVS, aes(x = abun_UVC, y = int_str)) +
  geom_point() +  # Agregar puntos
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # Línea de regresión
  labs(title = "UVCs abundance",
       x = "Abundance (maxN/m²)",
       y = "Interaction strength (interactions/2m²*hr)") +
  theme_bw()  +# Aplicar un tema minimalista
  theme(plot.title = element_text(hjust = 0.5)) +  # Centrar el título
  annotate("text", x = max(UVC_RUVS$abun_UVC) * 0.7, y = max(UVC_RUVS$int_str) * 0.9, 
           label = paste0("R = ", round(strength_abundance_cor3$estimate, 3)), 
           size = 5, hjust = 0) # Centrar el título # Centrar el título  # Centrar el título

strength_abundance_cor3 <- cor.test(UVC_RUVS$abun_UVC, UVC_RUVS$int_str, use = "complete.obs")

#3. Biomass----
biomass_plot3 <- ggplot(UVC_RUVS, aes(x = biom_UVC, y = int_str)) +
  geom_point() +  # Agregar puntos
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # Línea de regresión
  labs(
    title = "UVCs biomass",
    x = "Biomass (kg/m²)",
    y = "Interaction strength (interactions/2m²*hr)"
  ) +
  theme_bw() +  # Aplicar un tema minimalista 
  theme(plot.title = element_text(hjust = 0.5)) +  # Centrar el título
  annotate("text", x = max(UVC_RUVS$biom_UVC) * 0.7, y = max(UVC_RUVS$int_str) * 0.9, 
           label = paste0("R = ", round(strength_biomass_cor3$estimate, 3)), 
           size = 5, hjust = 0) # Centrar el título # Centrar el título

strength_biomass_cor3 <- cor.test(UVC_RUVS$biom_UVC, UVC_RUVS$int_str, use = "complete.obs")

grid.arrange(richness_plot3, abundance_plot3, biomass_plot3, ncol = 3)
