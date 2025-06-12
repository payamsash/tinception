library(readxl)
library(MatchIt)
library(openxlsx)
library(readr)
library(gtsummary)
library(ggplot2)

data <- read_csv("./material/no_talaska/behavioural/tinception_master_ready_for_matching.csv")
data$group <- as.factor(data$group)
data$sex <- as.factor(data$sex)
data$site <- as.factor(data$site)

## matching
match_out <- matchit(
                    group ~ age + sex + site + PTA, 
                    data = data, 
                    method = "nearest",
                    distance = "glm"
                    )

matched_data <- match.data(match_out) %>% 
    dplyr::filter(weights > 0)

write.xlsx(matched_data, "./material/no_talaska/behavioural/tinception_matched_optimal.xlsx")
write_lines(unique(matched_data$`subject ID`), "./material/no_talaska/behavioural/tinception_matched_optimal.txt")

## plotting
quartz()
plot(summary(match_out))
plot(match_out, type = "jitter", interactive = FALSE)
plot(match_out, type = "density", interactive = FALSE,
                    which.xs = ~age + site + sex)

## stat summery
summary_tbl <- matched_data %>%
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

summary_tbl