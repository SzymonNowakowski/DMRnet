
#wykonaj proszę te przypisania

n<-nrow(Xtruncated)
Xtruncated_normalized<- apply(Xtruncated, 2, function(x) sqrt(n/sum(x^2))*x)

mod<-glmnet::glmnet(Xtruncated, ytruncated, alpha = 1, intercept = TRUE, nlambda = 20, family = "gaussian")
mod_normalized<-glmnet::glmnet(Xtruncated_normalized, ytruncated, alpha = 1, intercept = TRUE, nlambda = 20, family = "gaussian")


#a na koniec wypisz sobie kolejno te rzeczy (najlepiej używając CTRL+ENTER)
mod$beta
mod_normalized$beta
which(Xtruncated_normalized[,4] != Xtruncated_normalized[,5])
