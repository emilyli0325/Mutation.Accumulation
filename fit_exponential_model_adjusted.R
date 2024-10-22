### fit 

## AT repeats

B <- structure(list(RepeatNo = c(7.17, 12.58679906, 17.39053206, 22.10804158), MutationRate = c(1.432335362, 6.361081674, 15.69207433, 30.56510681)), .Names = c("RepeatNo", "MutationRate"), row.names = c(1L, 2L,3L, 4L), class = "data.frame")


attach(B)

exponential.model <- lm(log(MutationRate)~ RepeatNo)

exponential.model$coefficients

liner.model <- lm(MutationRate~ RepeatNo)

anova(liner.model, exponential.model)


