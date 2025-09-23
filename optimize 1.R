xbar=13.4
n=15
s=16.6
margin_of_error=qt(0.975,14)*s/sqrt(n)
ll=xbar-margin_of_error
up=xbar+margin_of_error
ci=c(ll,up)
ci

margin_of_error2=qt(0.975,28)*s/sqrt(29)
margin_of_error2

mu0=2.06
t=(xbar-mu0)/(s/sqrt(n))
t
pval=pt(t,df=n-1,lower.tail=FALSE)
pval


margin_of_error=qt(0.99,24)*s/sqrt(n)
ll=xbar-margin_of_error
up=xbar+margin_of_error
ci=c(ll,up)
ci


ggplot2
Possums=read.csv('https://www.openintro.org/data/csv/possum.csv')
Possums
plot(skull_w~head_l,data=Possums)
attach(Possums)
Possums_new=Possums[skull_w < 66,]
Possums_new
line=lm(skull_w~head_l,Possums_new)
ggplot(Possums_new,aes(head_l,skull_w))
summary(line)#skull_w=1.25455+0.59851head_l
residuals=residuals(line)
RSS=sum(residuals^2)
RSS



beta_0=1.24
beta_1=0.62
Possums_new$predicted_skull_w=beta_0+beta_1*Possums_new$head_l
Possums_new$risiduals=(Possums_new$skull_w-Possums_new$predicted_skull_w)
new_Rss=sum(Possums_new$risiduals^2)
new_Rss


plot(hindtibia~winglength,data=Drosophila)
l1=lm(hindtibia~winglength,data=Drosophila)
summary(l1)
abline(l1)
b1=0.16653
SE_b1=0.00564
quantile=qt(p=0.995,df=427)
quantile
CI=b1+c(-1,1)*quantile*SE_b1
CI


b0=0.6417099
SE_b0=0.1222707
qt=qt(p=0.975,df=28)
ci=b0+c(-1,1)*qt*SE_b0
ci


p=pt(q=1.58,df=28,lower.tail=FALSE)
p
pval=2*p
pval





residuals=residuals(l1)
residuals
ebar=mean(residuals)
ebar
Rse=0.03475
quantile=qt(p=0.995,df=427)
quantile
CI=ebar+c(-1,1)*quantile*Rse
CI


par(mfrow = c(1,2))
plot(per2016~per2012,data=election20_samp1)
plot(per2016~per2012,data=election20_samp2)

mod1=lm(per2016~per2012,data=election20_samp1)
summary(mod1)#per2016=10.85946+0.88221per2012
mod2=lm(per2016~per2012,data=election20_samp2)
summary(mod2)#per2016=18.8997+0.7631per2012

confint(mod1,level=0.90)#lw:0.7473272. up:1.017097

b1=mod2$coefficients[2]
SE_b1=0.1072
quantile=qt(p=0.95,df=18)
CI=b1+c(-1,1)*quantile*SE_b1
CI



mod3=lm(Price~Decor,data=ItalianRestaurants)
summary(mod3)
prediction=predict(mod3,newdata=data.frame(Decor=20),interval="confidence",level=0.94)
prediction



par(mfrow = c(1,2))
plot(quality~helpfulness,data=Rateprof_1)
plot(quality~easiness,data=Rateprof_1)


mod1=lm(quality~helpfulness,data=Rateprof_1)
mod2=lm(quality~easiness,data=Rateprof_1)
anova(mod1)
anova(mod2)

l1=lm(quality~clarity,data=Rateprof_1)
plot(quality~clarity,data=Rateprof_1)
abline(l1)
par(mfrow = c(2,2))
plot(l1)

rstandard(l1)
rlist <- rstandard(l1)
which(rlist > 4) 
which(rlist < (-4))

cklist <- cooks.distance(l1)
which(cklist > 1)



cherry=read.csv("https://www.openintro.org/data/csv/cherry.csv")
par(mfrow=c(2,2))
plot(volume~height,data=cherry)
plot(log(volume) ~ log(height), data=cherry)
plot(1/volume ~ I(1/height), data=cherry)
plot(sqrt(volume) ~ sqrt(height), data=cherry)

mod1=lm(volume~height,data=cherry)
par(mfrow = c(2,2))
plot(mod1)

mod2=lm(log(volume) ~ log(height), data=cherry)
par(mfrow = c(2,2))
plot(mod2)
summary(mod2)

beta_1hat=mod2$coefficients[2]
beta1=1.01^(beta_1hat)-1 
beta1

install.packages("MASS")
library(MASS)
boxcox(volume~height, data=cherry)

boxcox(volume ~ log(height), data=cherry)

xtabs(hindtibia~ type + sex, data=Drosophila)/xtabs( ~ type + sex, data=Drosophila)

mod3=lm(hindtibia~ type + sex,data=Drosophila)
summary(mod3)

mod3_summary=summary(mod3)
rsqur=mod3_summary$r.squared
pro=1-rsqur
pro

length=mod3$coefficients[1] + mod3$coefficients[2]*1 + mod3$coefficients[3]*1
length

intercept=coef(mod3)[1]
typelab=coef(mod3)[2]
sexm=coef(mod3)[3]

fitted_length <- c(intercept + (typelab* 0) + (sexm * 0), 
                 intercept + (typelab * 1) + (sexm * 0),
                 intercept + (typelab * 0) + (sexm* 1),
                 intercept + (typelab * 1) + (sexm * 1))
groups = c("Group 1", "Group 2", "Group 3", "Group 4")
data = data.frame(Group = groups, length= fitted_length)
ordered_data = data[order(data$length, decreasing = TRUE), ]
print(ordered_data)



plot(mercury~temp,data=vapor)
mod1=lm(mercury ~ 1 + temp,data=vapor)
summary(mod1)

mod2=lm(mercury~1+temp+I(temp^2),data=vapor)
summary(mod2)
mod3=lm(mercury~1+temp+I(temp^2)+I(temp^3),data=vapor)
summary(mod3)

modA=lm(wage ~ 1 + age + union,data=cps1985)
summary(modA)
modB=lm(wage ~ 1 + age + union + age:union,data=cps1985)
summary(modB)


install.packages("leaps")
library(leaps)
a1=regsubsets( hindtibia ~ type + sex + winglength + wingwidth,data=Drosophila_1_)
summary(a1)

sum_a1=summary(a1)
n = nrow(Drosophila_1_) 
np = 2:5
bic = n*log(sum_a1$rss/n) + log(n)*np
which.min(bic)

m0=lm( hindtibia ~ 1,data=Drosophila_1_)
mFull=lm(hindtibia ~ 1 + morph + sex + winglength + wingwidth,data=Drosophila_1_)
step(mFull, k = log(n))

step(m0, direction="forward",scope = list(lower = m0, upper = mFull),k= log(n))


m1.1 = lm(hindtibia ~ sex, data = Drosophila_1_)
step(mFull, direction="backward", scope = list(lower = m1.1, upper = mFull), k = log(n))

count <- xtabs( ~ type, data = Drosophila_1_)
count
Sum_wage <- xtabs(winglength ~ type, data = Drosophila_1_)
Sum_wage
Sum_wage/count

pnorm(1.15,lower.tail=FALSE)
qnorm(0.01)

#a:H0:mu=2.37. Ha:mu not equal to 2.37
social_media_time=c(0,2,4,6,8,10,12)
t.test(social_media_time,alternative=c("two.side"),mu=2.37,conf.level = 0.99)
#p-value=0.067>0.01 do not reject H0

#b:99 percent confidence interval:(-0.0542, 12.0542) it includes 2.37, so it agree with our conclusion
#c:p-value=0.03  
#d:p-value=0.03 
#In the context of this problem, a Type II error means that:
#•	We concluded that the mean amount of time per day spent on social media for U of M freshmen and sophomores is 
#not different from the global mean of 2.37 hours.
#•	However, in reality, the true mean amount of time might actually be different from 2.37 hours, and we failed to detect this difference.

#Implication: This would mean that our sample did not provide strong enough evidence to show a difference when one actually exists, 
#potentially leading us to overlook a genuine difference in social media usage between U of M students and the global average.



pairs(price~carat+depth,data=FairDiamonds)
par(mfrow=c(2,2))
mod1=lm(price~carat+depth,data=FairDiamonds)
plot(mod1)

mod2=lm(log(price)~log(carat)+depth,data=FairDiamonds)
par(mfrow=c(2,2))
plot(mod2)

summary(mod2)
1.01^mod2$coefficients[2]-1

exp(mod2$coefficients[3])-1

nrow(TitanicPartial_v2)
table(TitanicPartial_v2$Survival)


mod3=glm( Survival ~ 1 + as.factor(Pclass),family = binomial,data=TitanicPartial_v2)
summary(mod3)

predict(mod3,newdata = data.frame(as.factor(Pclass)=c(3)), type = 'response')

mod4=glm(Survival ~ 1 + as.factor(Pclass)+ Age,family = binomial,data=TitanicPartial_v2)
summary(mod4)

exp(mod4$coefficients[4])-1

predict(mod4,newdata = data.frame(Pclass=c(3),Age=c(20)), type = 'response')

mod3=glm(Survival~ 1 + as.factor(Pclass)+ Age + as.factor(Pclass): Age,family = binomial,data=TitanicPartial_v2)
summary(mod3)

exp(-0.0404)-1

predict(mod3, newdata = data.frame(Age = c(70,40,10), 
                                   Pclass = c(1,2,3)), type = "response")
mod2=glm(Survival ~ 1 + as.factor(Pclass)+ Age,family = binomial,data=TitanicPartial_v2)
summary(mod2)


odds1 <- exp(-14.379 -3.907*s )
odds2 <- exp(-14.379-3.907*(s+0.1) )
(odds_ratio21 <- odds2/odds1)

pchisq(6.0217,df=1,lower.tail=FALSE)




s1 <- arima.sim(model = list(ar = 0.7), n = 50, sd = 0.5)+2
s2 <- arima.sim(model = list(ar = 0.7), n = 1000, sd = 0.5)+2
par(mfrow = c(2, 1))
plot.ts(s1) 
plot.ts(s2) 
mean(s1)
mean(s2)

var(s1)
var(s2)
mean_s1 <- mean(s1)
mean_s2 <- mean(s2)
mean_s1  # Sample mean for n = 50
mean_s2  # Sample mean for n = 1000

var_s1 <- var(s1)
var_s2 <- var(s2)
var_s1  # Sample variance for n = 50
var_s2  # Sample variance for n = 1000

acf(s1, plot = FALSE)
acf(s2, plot = FALSE)
acf_s1$acf[5]  # Sample autocorrelation at lag 4 for n = 50
acf_s2$acf[5]  # Sample autocorrelation at lag 4 for n = 1000






plot.ts(GoogleStockVolume2020$Volume)

par(mfrow=c(1,2))
acf(GoogleStockVolume2020$Volume)
pacf(GoogleStockVolume2020$Volume)


mod1 = arima(GoogleStockVolume2020$Volume, order = c(1,0,0))
mod1

acf(mod1$residuals)

predict(mod1, n.ahead = 1)


y75=2219267.4 + 0.8 * (1949000 - 2219267.4) 
y75

residual=1695500-2004816
residual



plot.ts(SampleData_HW11$value )

xt <- diff(SampleData_HW11$value)
plot.ts(xt)

par(mfrow=c(1,2))
acf(xt)
pacf(xt)


(modx <- arima(xt, order = c(0,0,1)))


predict(modx, n.ahead = 1)
x1000=-0.3574443
y1001=x1000+466.7438
y1001

y1000=466.7438
y999=467.9156
x999=y1000-y999
x1000=0.466+0.7918*x999
y1001=x1000+y1000

y1001


install.packages("tinytex")
tinytex::install_tinytex()  # 安装 TinyTeX

tinytex::is_tinytex()
tinytex::tlmgr_install("some-missing-package")


install.packages("rmarkdown")
library(rmarkdown)

rmarkdown::render("your_file.Rmd", output_format = "pdf_document")


set.seed(123)
n_values=c(5, 10, 15, 100)
theta=4
reps=10000
simulate_coverage=function(n, theta, reps=10000) {
  captured.list=numeric(reps)  
  
  for (i in 1:reps) {
    X_sample=invf(n, theta) 
    X_bar=mean(X_sample)  
    S=sd(X_sample)  
    
    t_val=qt(0.975, df=n-1)  
    lower=X_bar - t_val * S / sqrt(n)  
    upper=X_bar + t_val * S / sqrt(n)  
    captured.list[i]=1 * (lower <= 3) * (3 <= upper)  
  }
  
  coverage_results=mean(captured.list)
}
data.frame(n = n_values, Coverage_Probability = coverage_results)




set.seed(123)
n_values=c(5, 10, 15, 100)
theta=4
reps=10000
simulate_coverage=function(n, theta, reps=10000) {
  captured.list=numeric(reps)
  
  for (i in 1:reps) {
    X_sample=invf(n, theta) 
    X_bar=mean(X_sample)  
    S=sd(X_sample)  
    
    t_val=qt(0.975, df=n-1)  
    lower=X_bar - t_val * S / sqrt(n)  
    upper=X_bar + t_val * S / sqrt(n)  
    captured.list[i]=1 * (lower <= 3) * (3 <= upper)  
  }
  
  return(mean(captured.list))
}

data.frame(n=n_values,coverage_pro=mean(captured.list))


install.packages("tinytex")
tinytex::install_tinytex()  

tinytex::reinstall_tinytex()

rmarkdown::render("your_file.Rmd", output_format = "pdf_document")

tinytex::tlmgr_install("latexmk")
latexmk --version

system("latexmk --version")

tinytex::tlmgr_install("missing_package_name")
getwd()
rmarkdown::render("Feng.Rmd", output_format = "pdf_document", clean = TRUE)


#1.1
theta=3
x=seq(0, theta, length.out=100)
f_x=3*x^2/theta^3
plot(x,f_x,xlab="x", ylab="f(x;3)",main="Density Function f(x;3)")

#1.2
#Calculate the CDF:F(x)=x^3 / theta^3
#Set F(x)=U , then X=theta*U^1/3
#To generate a realization of X, take a uniform random variable #U~Unif(0,1) and compute  X=theta*U^1/3

#1.3
invf=function(n, theta) {
  U=runif(n) 
  X=theta * U^(1/3) 
  return(X)
}

set.seed(123)  
samples=invf(10000, 4)

hist(samples, probability=TRUE,main="Histogram of Generated Samples (Theta = 4)",xlab="X", ylab="Density")

x_vals=seq(0, 4, length.out=100)
lines(x_vals, 3 * x_vals^2 / 4^3)

#1.4
#we calculate E(x)=3 when theta=4
mean(samples) #2.995
#So,the simulation-based estimate close to the formula

#1.5
set.seed(123)
n_values=c(5, 10, 15, 100)
theta=4
reps=10000
simulate_coverage=function(n, theta, reps=10000) {
  captured.list=numeric(reps)
  
  for (i in 1:reps) {
    X_sample=invf(n, theta) 
    X_bar=mean(X_sample)  
    S=sd(X_sample)  
    
    t_val=qt(0.975, df=n-1)  
    lower=X_bar - t_val * S / sqrt(n)  
    upper=X_bar + t_val * S / sqrt(n)  
    captured.list[i]=1 * (lower <= 3) * (3 <= upper)  
  }
  
  return(mean(captured.list))
}
coverage_results = sapply(n_values, function(n) simulate_coverage(n, theta, reps))

result_df = data.frame(n = n_values, coverage_prob = coverage_results)
print(result_df)
#result close to 0.95

#1.6
#Var(X)=E(X^2)−(E(X))^2 ,E(X^2)=3*theta^2/5 , (E(X))^2=(3*theta/4)^2
#then theta=4, Var(X)= 1.8
var(samples) #which is 0.596 not equal to the formula variance.


#2.1
generate_normal_box_muller=function(n, mu=37, sigma=2) {
  U1=runif(n/2)
  U2=runif(n/2)
  Z1=sqrt(-2 * log(U1)) * cos(2 * pi * U2)
  Z2=sqrt(-2 * log(U1)) * sin(2 * pi * U2)
  Z=c(Z1, Z2) 
  X=mu + sigma * Z
  return(X)
}
normal.method1=generate_normal_box_muller(100)


generate_normal_clt=function(n, mu=37, sigma=2) {
  U=matrix(runif(n * 12), nrow=n) 
  Z=rowSums(U) - 6 
  X=mu + sigma * Z
  return(X)
}

normal.method2=generate_normal_clt(100)

#2.2
simulate_coverage_normal=function(n, mu=37, sigma=2, reps=10000) {
  captured=numeric(reps)
  
  for (i in 1:reps) {
    X_sample=generate_normal_box_muller(n, mu, sigma)
    X_bar=mean(X_sample)
    S=sd(X_sample)
    t_val=qt(0.975, df=n-1)
    lower=X_bar - t_val * S / sqrt(n)
    upper=X_bar + t_val * S / sqrt(n)
    
    captured[i]=1 * (lower <= mu) * (mu <= upper)
  }
  
  coverage_results=mean(captured)
}

coverage_probability=simulate_coverage_normal(n=100)
print(coverage_probability)  # it close to 0.95

#2.3
simulate_coverage_gamma=function(n, alpha, reps=10000) {
  captured=numeric(reps)
  mu=2 * 2  
  sigma=sqrt(2 * 2^2) 
  
  for (i in 1:reps) {
    X_sample=rgamma(n, shape=2, scale=2)
    X_bar=mean(X_sample)
    S=sd(X_sample)
    
    t_val=qt(1 - alpha/2, df=n-1)
    lower=X_bar - t_val * S / sqrt(n)
    upper=X_bar + t_val * S / sqrt(n)
    
    captured[i]=1*(lower <= mu) * (mu <= upper)
  }
  
  coverage_results=mean(captured)
}

n_values=c(10, 50, 500, 5000)
alpha_values=c(0.01, 0.05)

coverage_results=outer(n_values, alpha_values, Vectorize(function(n, alpha) {simulate_coverage_gamma(n, alpha)}))

print(coverage_results) 

#Expected coverage close to 1-alpha

#2.4
# t= (xbar-mu0)/(sigma/sqrt(n))

#2.5
#If the null hypothes is true,it is the exact t-distribution.

#2.6
est.pval.dist.P2=function(n, mu0, mu, sigma, reps=10000) {
  pvals=numeric(reps)
  
  for (i in 1:reps) {
    X_sample=generate_normal_box_muller(n, mu, sigma)
    X_bar=mean(X_sample)
    S=sd(X_sample)
    
    T_stat=(X_bar - mu0) / (S / sqrt(n))
    pvals[i]=1 - pt(T_stat, df=n-1)
  }
  
  return(pvals)
}

p_values_H0=est.pval.dist.P2(n=30, mu0=37, mu=37, sigma=2)
hist(p_values_H0, main="P-value Distribution under H0", xlab="p-value")

p_values_HA=est.pval.dist.P2(n=30, mu0=37, mu=38, sigma=2)
hist(p_values_HA, main="P-value Distribution under HA", xlab="p-value")
#Under H0, the p-value should be uniformly distributed; under HA, the #p-value should be biased toward 0.


#2.7
extra.credit.pvals=function(n, mu0, mu, sigma, reps=10000) {
  pvals=numeric(reps)
  
  for (i in 1:reps) {
    X_sample=generate_normal_box_muller(n, mu, sigma)
    X_bar=mean(X_sample)
    S=sd(X_sample)
    
    t_val=qt(0.975, df=n-1)
    lower=X_bar - t_val * S / sqrt(n)
    upper=X_bar + t_val * S / sqrt(n)
    
    pvals[i]= 1*(mu0 < lower) 
  }
  
  mean(pvals)  
}

print(extra.credit.pvals(n=30, mu0=37, mu=37, sigma=2))  
# 0.02<0.05,mu0 is outside the confidence interval.

tinytex::tlmgr_install("graphicx")





# Simulate the power curve
power_curve <- function(n1 = 50, alpha = 0.01) {
  power_vals <- numeric(81)
  n2_vals <- 20:100  # Sample sizes for Y (n2)
  
  for (n2 in n2_vals) {
    rejections <- 0
    for (i in 1:1000) {  # Simulate 1000 experiments
      samples <- sim.two.gamma(n1, n2)
      X_bar <- mean(samples$X)
      Y_bar <- mean(samples$Y)
      var_X <- var(samples$X)
      var_Y <- var(samples$Y)
      
      T <- (X_bar - Y_bar) / sqrt(var_X / n1 + var_Y / n2)  # Compute test statistic
      df <- n1 + n2 - 2  # Degrees of freedom
      
      p_value <- 2 * (1 - pt(abs(T), df))  # Two-tailed p-value
      if (p_value < alpha) {
        rejections <- rejections + 1  # Reject null hypothesis
      }
    }
    power_vals[n2 - 19] <- rejections / 1000  # Power for this n2
  }
  
  plot(n2_vals, power_vals, type = "b", col = "blue", main = "Power Curve", xlab = "n2", ylab = "Power")
  return(power_vals)
}


X=rgamma(n1, shape = mu1, scale = 1)  
Y=rgamma(n2, shape = mu2 , scale =2 )

T_stat= t_test_statistic(X, Y)
T_stat




# 安装必要包（如果未安装）
install.packages(c("cluster", "factoextra", "ggplot2", "scales"))

library(cluster)
library(factoextra)
library(ggplot2)
library(scales)

# 设置随机种子
set.seed(42)

# 生成模拟数据（类似make_blobs）
library(MASS)
centers <- matrix(c(1,1,1, -1,-1,-1, 1,-1,1, -1,1,-1, 0,0,0), ncol=3, byrow=TRUE)
X <- NULL
for (i in 1:5) {
  X <- rbind(X, mvrnorm(n = 60, mu = centers[i,], Sigma = diag(3)))
}

# 标准化数据
X_scaled <- scale(X)

# 初始化结果存储
k_range <- 2:8
sse <- numeric()
silhouette_scores <- numeric()

# 循环计算不同K值的聚类结果
for (k in k_range) {
  kmeans_result <- kmeans(X_scaled, centers = k, nstart = 20)
  sse[k - 1] <- kmeans_result$tot.withinss
  
  sil <- silhouette(kmeans_result$cluster, dist(X_scaled))
  silhouette_scores[k - 1] <- mean(sil[, 3])
}

# 构造数据框用于绘图
df <- data.frame(
  K = k_range,
  SSE = sse,
  Silhouette = silhouette_scores
)

# 作图：肘部法则与轮廓系数
par(mfrow = c(1, 2))  # 两图并排

# 肘部法则图
plot(df$K, df$SSE, type = "b", pch = 19, col = "red",
     xlab = "Number of clusters (K)", ylab = "Sum of Squared Errors (SSE)",
     main = "Elbow Method for Optimal K")
abline(v = 5, col = "blue", lty = 2)

# 轮廓系数图
plot(df$K, df$Silhouette, type = "b", pch = 19, col = "blue",
     xlab = "Number of clusters (K)", ylab = "Silhouette Score",
     main = "Silhouette Analysis for Optimal K")
abline(v = 5, col = "red", lty = 2)

# 恢复默认图形设置
par(mfrow = c(1, 1))




# 安装并加载绘图包
install.packages("plotly")
library(plotly)

# 使用前三个主成分
pca_result <- prcomp(scale(USArrests))
pca_scores <- as.data.frame(pca_result$x[, 1:3])

# 添加聚类结果（例如4类）
set.seed(123)
kmeans_result <- kmeans(pca_scores, centers = 4)
pca_scores$cluster <- as.factor(kmeans_result$cluster)

# 绘制三维散点图
plot_ly(pca_scores, x = ~PC1, y = ~PC2, z = ~PC3,
        color = ~cluster, colors = "Set2",
        type = "scatter3d", mode = "markers") %>%
  layout(title = "图4-1：因子得分散点图（三维）")


# 加载包
library(cluster)
library(factoextra)

# 计算距离矩阵
dist_matrix <- dist(scale(USArrests))

# K-means 聚类
kmeans_result <- kmeans(scale(USArrests), centers = 4, nstart = 20)

# 轮廓图绘制
fviz_silhouette(silhouette(kmeans_result$cluster, dist_matrix)) +
  ggtitle("图4-2：聚类轮廓图") +
  theme_minimal()













