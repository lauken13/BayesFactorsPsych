library(tidyverse)

final_results <- readRDS("results/final_case2_results.rds")

logo_cv <- final_results[[1]]

loo_cv <- final_results[[2]]


logo_cv_cases <- rbind(data.frame(logo_cv[[1]],case = "No fixed effect/random slopes"),
                      data.frame(logo_cv[[2]],case = "Fixed effect only"),
                      data.frame(logo_cv[[3]],case = "Fixed effect and random slope"))

logo_cv_cases %>%
  pivot_longer(elpd_diff.mod6:se_diff.mod3, 
               names_to = c("measure","model"),
               names_sep = "\\.") %>% 
  pivot_wider(names_from = "measure")%>%
  mutate(data = factor(data, levels = c("Fully Aggr.",    "2 Items","5 Items", "10 Items", "20 items","50 Items","Full Data"
  )),
  case = factor(case, levels = c("No fixed effect/random slopes","Fixed effect only",
                                 "Fixed effect and random slope")),
  model = fct_recode(model, "Model 3" = "mod3",
                     "Model 4" = "mod4", 
                     "Model 5" = "mod5"))%>%
  filter(model!="mod6")%>%
  ggplot(.,aes(x=data, y = elpd_diff, colour = model))+
  geom_hline(aes(yintercept = 0), colour = "orange")+
  geom_point(position = position_dodge(width = .5))+
  geom_errorbar(aes(ymin = elpd_diff - 2*se_diff,
                    ymax=elpd_diff +2*se_diff),
                position = position_dodge(width = .5))+
  facet_wrap(case~.)+
  coord_flip() +
  theme_bw()+
  scale_colour_viridis_d(end=.9)+
  ylab("ELPD difference, compared to model 6")+
  theme(axis.title.y = element_blank(),
        legend.title = element_blank())


loo_cv_cases <- rbind(data.frame(loo_cv[[1]],case = "No fixed effect/random slopes"),
                       data.frame(loo_cv[[2]],case = "Fixed effect only"),
                       data.frame(loo_cv[[3]],case = "Fixed effect and random slope"))

loo_cv_cases %>%
  pivot_longer(elpd_diff.mod6:se_diff.mod3, 
               names_to = c("measure","model"),
               names_sep = "\\.") %>% 
  pivot_wider(names_from = "measure")%>%
  mutate(data = factor(data, levels = c("Fully Aggr.",    "2 Items","5 Items", "10 Items", "20 items","50 Items","Full Data"
  )),
  case = factor(case, levels = c("No fixed effect/random slopes","Fixed effect only",
                                 "Fixed effect and random slope")),
  model = fct_recode(model, "Model 3" = "mod3",
                     "Model 4" = "mod4", 
                     "Model 5" = "mod5"))%>%
  filter(model!="mod6")%>%
  ggplot(.,aes(x=data, y = elpd_diff, colour = model))+
  geom_hline(aes(yintercept = 0), colour = "orange")+
  geom_point(position = position_dodge(width = .5))+
  geom_errorbar(aes(ymin = elpd_diff - 2*se_diff,
                    ymax=elpd_diff +2*se_diff),
                position = position_dodge(width = .5))+
  facet_wrap(case~.)+
  coord_flip() +
  theme_bw()+
  scale_colour_viridis_d(end=.9)+
  ylab("ELPD difference, compared to model 6")+
  theme(axis.title.y = element_blank(),
        legend.title = element_blank())
