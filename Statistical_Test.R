library(MASS)
# chisq.test
dt <- cbind(c(2,8,3,14,8,2),c(10.03,2.30,10.03,10.03,2.30,2.30))
dt <- cbind(c(17,5),c(20.62,2.51))
dt <- cbind(c(17,10,5,5),c(16,47,39,38))
dt <- cbind(c(9,6,4,3),c(12.3,4.6,1.5,3.6))
chisq.test(dt)
dt <- cbind(c(16,47,39,38),c(22.4,55.7,21.6,40.3))

dt <- cbind(c(6,5,3,3),c(16,31,26,23))


# two sample Wilcoxon rank sum 
data <- read.csv("Micro_repeatNo_genomicLocation.csv", header=TRUE)
wilcox.test(exp~GenomeLocation, data=data[data$GenomeLocation== c("CDS", "intergenic"),], exact = FALSE)

wilcox.test(Repeat.No~DIFF, data=data[data$Motif.Size==1,], exact = FALSE)


wilcox.test(Repeat.No.Change2~DIFF, data=micro.ins.del, exact = FALSE)

## ins & del copy number difference
setwd("C:/Users/xli/Desktop/06062018/filterMethodsOptimize/rm5loci")
data <- read.csv("micro_mut.csv", header=TRUE)
> wilcox.test(Repeat.No.Change2~DIFF, data=data, alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  Repeat.No.Change2 by DIFF
W = 2812.5, p-value = 0.03173
alternative hypothesis: true location shift is greater than 0

## <15bp & >= 15bp


# binomial test
binom.test(x, n, p=0.5, alternative = c("two.sided", "less", "greater"), conf.level = 0.95)

binom.test(10, 37, p=0.33, alternative = c("two.sided", "less", "greater"), conf.level = 0.95)

binom.test(77, 140, p=0.3615, alternative = c("two.sided"), conf.level = 0.95)


### 95% CI & se 
library(mosaic)
boot_dsn <- do(1000)*rflip(19435,prob=16/19435)
histogram(~prop,data=boot_dsn,col="grey")
se <- sd(~prop, data=boot_dsn)
> se
[1] 0.0002054791
confint(boot_dsn, level = 0.95, method="quantile")



boot_dsn <- do(1000)*rflip(2641.5,prob=17/2641.5)
histogram(~prop,data=boot_dsn,col="grey")
se <- sd(~prop, data=boot_dsn)
> se
[1] 0.0002054791
confint(boot_dsn, level = 0.95, method="quantile")


sd(~prop, data=boot_dsn)
confint(boot_dsn, level = 0.95, method="quantile")

boot_dsn <- do(1000)*rflip(11580293,prob=6/11580293)

boot_dsn <- do(1000)*rflip(3883857,prob=1/3883857)
boot_dsn <- do(1000)*rflip(3883857,prob=4/3883857)
boot_dsn <- do(1000)*rflip(16898216,prob=1/16898216)
boot_dsn <- do(1000)*rflip(16898216,prob=0/16898216)
boot_dsn <- do(1000)*rflip(16898216,prob=9/16898216)



boot_dsn <- do(1000)*rflip(46626,prob=13/46626)
se <- sd(~prop, data=boot_dsn)
> se
[1] 7.869357e-05

boot_dsn <- do(1000)*rflip(43903,prob=69/43903)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(13106,prob=16/13106)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(10689,prob=23/10689)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(4975,prob=10/4975)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(2152,prob=9/2152)
se <- sd(~prop, data=boot_dsn)
se

boot_dsn <- do(1000)*rflip(121451,prob=142/121451)
se <- sd(~prop, data=boot_dsn)
se


boot_dsn <- do(1000)*rflip(21899,prob=9/21899)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(39923,prob=35/39923)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(31896,prob=41/31896)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(15229,prob=40/15229)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(6754,prob=16/6754)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(115701,prob=142/115701)
se <- sd(~prop, data=boot_dsn)
se


boot_dsn <- do(1000)*rflip(23405,prob=9/23405)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(12297,prob=21/12297)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(6053,prob=25/6053)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(1828,prob=15/1828)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(43583,prob=69/43583)
se <- sd(~prop, data=boot_dsn)
se


boot_dsn <- do(1000)*rflip(19435,prob=22/19435)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(48304,prob=46/48304)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(18719,prob=38/18719)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(34993,prob=35/34993)
se <- sd(~prop, data=boot_dsn)
se
boot_dsn <- do(1000)*rflip(121451,prob=142/121451)
se <- sd(~prop, data=boot_dsn)
se


boot_dsn <- do(1000)*rflip(123722,prob=145/123722)
se <- sd(~prop, data=boot_dsn)
se



binom.test(6, 17, p=0.5572, alternative = c("two.sided"), conf.level = 0.95)





