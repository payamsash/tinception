library(readxl)
library(MatchIt)
library(openxlsx)
library(readr)
library(gtsummary)

data <- read_csv("./data/tinception_master_ready_for_matching.csv")
data$group <- as.factor(data$group)
data$sex <- as.factor(data$sex)
data$site <- as.factor(data$site)

## matching
match_out <- matchit(
                    group ~ age + sex + site, 
                    data = data, 
                    method = "nearest",
                    distance = "glm"
                    )

matched_data <- match.data(match_out) %>% 
    dplyr::filter(weights > 0)

write.xlsx(matched_data, "./data/tinception_matched.xlsx")
write_lines(unique(matched_data$`subject ID`), "./data/matched_subjects.txt")

## plotting
quartz()
plot(summary(match_out))
plot(match_out, type = "jitter", interactive = FALSE)
plot(match_out, type = "density", interactive = FALSE,
                    which.xs = ~age + sex + site)

## stat summery
summary_tbl <- matched_data %>%
    select(group, age, sex, site) %>%
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