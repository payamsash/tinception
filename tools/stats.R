library(readxl)
library(MatchIt)
library(openxlsx)
library(readr)
library(gtsummary)
library(ggplot2)
library(dplyr)

data <- read_csv("/Users/payamsadeghishabestari/tinception/material/with_qc/behavioural/tinception_master_ready_for_matching_qc.csv")
data <- na.omit(data)
data$group <- as.factor(data$group)
data$sex <- as.factor(data$sex)

# data <- data %>%
#   mutate(PTA4_avg = rowMeans(select(., RPTA4, LPTA4), na.rm = TRUE)) %>%
#   select(-RPTA4, -LPTA4)

# data <- data %>%
#   mutate(HF_avg = rowMeans(select(., RHF, LHF), na.rm = TRUE)) %>%
#   select(-RHF, -LHF)


data$site <- as.factor(data$site)

## matching
match_out <- matchit(
                    group ~ age + sex + PTA + site, 
                    data = data, 
                    method = "optimal",
                    distance = "glm",
#                    ratio=2
                    )

matched_data <- match.data(match_out) %>% 
    dplyr::filter(weights > 0)

# matched_data <- matched_data %>%
#   group_by(subclass) %>%
#   filter(all(distance < 0.8)) %>%
#   ungroup()

# write.xlsx(matched_data, "./material/with_qc/behavioural/tinception_matched_optimal.xlsx")
# write_lines(unique(matched_data$`subject ID`), "./material/with_qc/behavioural/tinception_matched_optimal.txt")

# write.xlsx(matched_data, "/Users/payamsadeghishabestari/Downloads/matched_optimal_glm.xlsx")

## plotting
quartz()
plot(summary(match_out))
plot(match_out, type = "jitter", interactive = FALSE)
plot(match_out, type = "density", interactive = FALSE,
                    which.xs = ~age + sex + PTA)

## stat summery
summary_tbl <- data %>%
    select(group, age, PTA, sex, site) %>%
    tbl_summary(
        by = group,  
        statistic = list(
            all_continuous() ~ "{mean} ± {sd}",  
            all_categorical() ~ "{n} ({p}%)"     
        ),
        digits = all_continuous() ~ 1  
    ) %>%
    add_p() %>%  
    bold_labels()  

# summary_tbl
