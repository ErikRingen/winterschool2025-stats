fit_model <- function(data, age_moderator = TRUE){
    if(age_moderator){
        model <- lm(cognitive_flexibility ~ multilingual * age, data = data)
    } else {
        model <- lm(cognitive_flexibility ~ multilingual, data = data)
    }
    return(model)
}
