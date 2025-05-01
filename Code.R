library(tidyverse)
library(e1071)

####### a)

finch.dat = read_csv("zebrafinches.csv")

t = t.test(finch.dat$further, alternative = "less")$statistic
skew = skewness(finch.dat$further)
pdf = dnorm(t)

error = (skew/sqrt(length(finch.dat$further)))*(((2*t^2) + 1)/6)*(pdf)


###### b)

calcError = function(t){
  error = (skew/sqrt(25))*(((2*t^2) + 1)/6)*(dnorm(t))
}

graph.dat = tibble(t = seq(-10, 10, length.out = 1000),
                   error = calcError(t))

ggplot(data = graph.dat) +
  geom_line(aes(x = t, y = error)) + 
  theme_bw()


###### c)

t = unname(t.test(finch.dat$further, alternative = "less")$statistic)
skew = skewness(finch.dat$further)
pdf = dnorm(t)

n.desired1 = ((skew/(6*0.10*0.05))*(((2*t^2)+1)*pdf))^2

t = qnorm(0.05)
skew = skewness(finch.dat$further)
pdf = dnorm(t)

n.desired2 = ((skew/(6*0.10*0.05))*(((2*t^2)+1)*pdf))^2


##### PART 2
##### a)


resample.t = function(data){
  t = c()
  for (i in 1:10000){
    set.seed(7272+i)
    samps = sample(data, length(data), replace = TRUE)
    x.bar = mean(samps)
    t[i] = x.bar/(sd(data)/sqrt(length(data)))
  }
  t.shifted = mean(t) - t
  return(list(shift = t.shifted, t = t))
}

t.stats = tibble(further.shifted = resample.t(finch.dat$further)$shift,
                 further.t = resample.t(finch.dat$further)$t,
                 closer.shifted = resample.t(finch.dat$closer)$shift,
                 closer.t = resample.t(finch.dat$closer)$t,
                 difference.shifted = resample.t(finch.dat$diff)$shift,
                 difference.t = resample.t(finch.dat$diff)$t)


###### b)

p.calc = function(data, data2, alt){
    if (alt == "less"){
      t.star = t.test(data2, alternative = "less")$statistic
      for (i in 1:10000){
        p = length(which(data <= t.star))/length(data)
      }
    } else if (alt == "greater"){
      t.star = t.test(data2, alternative = "greater")$statistic
      for (i in 1:10000){
        p = length(which(data >= t.star))/length(data)
      }
    } else if (alt == "two.sided"){
      t.star = t.test(data2, alternative = "two.sided")$statistic
      for (i in 1:10000){
        p = length(which(data >= t.star))/length(data)
      }
    }
  
  return(p)
}

p.stats = tibble(further = p.calc(t.stats$further.shifted, finch.dat$further, "less"),
                 closer = p.calc(t.stats$closer.shifted, finch.dat$closer, "greater"),
                 difference = p.calc(t.stats$difference.shifted, finch.dat$diff, "two.sided"))


###### c), d)

quantiles = tibble(further = quantile(t.stats$further.t, 0.05),
                   closer = quantile(t.stats$closer.t, 0.05),
                   difference = quantile(t.stats$difference.t, 0.05))


conf = tibble(further = quantile(t.stats$further.t, c(0.025, 0.975)),
                   closer = quantile(t.stats$closer.t, c(0.025, 0.975)),
                   difference = quantile(t.stats$difference.t, c(0.025, 0.975)))
###### PART 3
###### a)

rand = tibble(x = rep(NA, 1000))

x.shift = finch.dat$further

for(i in 1:1000){
  curr.rand = x.shift*sample(c(-1, 1), 25, replace = TRUE)
  rand$x[i] = mean(curr.rand)
}

delta = abs(mean(finch.dat$further))

low = -delta
high = delta

p = length(which(rand$x <= low | rand$x >= high))



mu.a = mean(finch.dat$further)
mu0 = mu.a
p = 1

x.shift = finch.dat$further - mu.a  # center ONCE outside the loop

while (p > 0.05) {
  mu0 = mu0 - 0.0001
  rands = numeric(1000)
  for (i in 1:1000) {
    curr.rand = x.shift * sample(c(-1, 1), 25, replace = TRUE) + mu0
    rands[i] = mean(curr.rand)
  }
  p = mean(rands >= mu.a)
}

lower = mu0


mu.a = mean(finch.dat$further)
mu0 = lower
p = 0

x.shift = finch.dat$further - mu.a  # center ONCE outside the loop

while (p > 0.05 | p == 0) {
  mu0 = mu0 + 0.0001
  rands = numeric(1000)
  x.shift = finch.dat$further - mu.a
  for (i in 1:1000) {
    curr.rand = mu0 + x.shift * sample(c(-1, 1), 25, replace = TRUE)
    rands[i] = mean(curr.rand)
  }
  p = mean(rands <= mu.a)
}

upper = mu0


randomizer = function(data, bound){
  if (bound == 'lower'){
    mu.a = mean(data)
    mu0 = mu.a
    p = 1
    
    x.shift = data - mu.a  # center ONCE outside the loop
    
    while (p > 0.05) {
      mu0 = mu0 - 0.0001
      rands = numeric(1000)
      for (i in 1:1000) {
        curr.rand = x.shift * sample(c(-1, 1), 25, replace = TRUE) + mu0
        rands[i] = mean(curr.rand)
      }
      p = mean(rands >= mu.a)
    }
    
    lower = mu0
    return(lower)
  }else if (bound == 'upper'){
    mu.a = mean(data)
    mu0 = mu.a
    p = 0
    
    x.shift = data - mu.a  # center ONCE outside the loop
    
    while (p > 0.05 | p == 0) {
      mu0 = mu0 + 0.0001
      rands = numeric(1000)
      x.shift = finch.dat$further - mu.a
      for (i in 1:1000) {
        curr.rand = mu0 + x.shift * sample(c(-1, 1), 25, replace = TRUE)
        rands[i] = mean(curr.rand)
      }
      p = mean(rands <= mu.a)
    }
    
    upper = mu0
    return(upper)
  }
}

randomized = tibble(data = c("Further", "Closer", "Difference"),
                    lower = c(lower, randomizer(finch.dat$closer, 'lower'),
                              randomizer(finch.dat$diff, 'lower')),
                    mean = c(mean(finch.dat$further), mean(finch.dat$closer),
                             mean(finch.dat$diff)),
                    upper = c(upper, randomizer(finch.dat$closer, 'upper'),
                              randomizer(finch.dat$diff, 'upper')))


