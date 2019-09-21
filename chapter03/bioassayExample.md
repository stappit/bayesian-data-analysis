Bioassay Experiment
================
Corrie
5/29/2019

``` r
library(MASS)
library(ggplot2)
library(mvtnorm)
library(dplyr)
library(cowplot)
```

This example is a nonconjugate model for a bioassay experiment. It is a
two-parameter example from the class of generalized linear models. Here,
we use a simulation grid-approximation approach to get the posterior
distribution.

## The data

To test how toxic a drug is, it is often given to animals at various
doses and then observed how many have adverse outcomes. Often, the
response is simply a dichothomous outcome: animal alive or dead. An
example of such data:

``` r
d <- data.frame(log.dose=c(-0.86, -0.3, -0.05, 0.73),
                n.animals=c(5,5,5,5),
                n.deaths=c(0,1,3,5))
d %>% knitr::kable()
```

| log.dose | n.animals | n.deaths |
| -------: | --------: | -------: |
|   \-0.86 |         5 |        0 |
|   \-0.30 |         5 |        1 |
|   \-0.05 |         5 |        3 |
|     0.73 |         5 |        5 |

``` r
d <- d %>%
  mutate(prop=n.deaths / n.animals) 
d %>%
  ggplot(aes(x=log.dose, y=prop)) +
  geom_point(col="#E41A1C") +
  theme_minimal() +
  labs(y="Proportion of deaths",
       x="Dose (log g/ml")
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## The model

It is reasonable to model the outcomes of the five animals *within each
group* as exchangeable and independent. The data points
![y\_i](https://latex.codecogs.com/png.latex?y_i "y_i")=`n.deaths` are
then binomially distributed:   
![y\_i | \\theta\_i \\sim \\text{Bin}(n\_i,
\\theta)](https://latex.codecogs.com/png.latex?y_i%20%7C%20%5Ctheta_i%20%5Csim%20%5Ctext%7BBin%7D%28n_i%2C%20%5Ctheta%29
"y_i | \\theta_i \\sim \\text{Bin}(n_i, \\theta)")  
where ![n\_i](https://latex.codecogs.com/png.latex?n_i "n_i") is the
number of animals (in this example it is 5 for each group).

We model the response
![\\theta](https://latex.codecogs.com/png.latex?%5Ctheta "\\theta") by
the dose using a linear model together with a logit-link:   
![\\text{logit}(\\theta\_i) = \\alpha + \\beta
x\_i](https://latex.codecogs.com/png.latex?%5Ctext%7Blogit%7D%28%5Ctheta_i%29%20%20%3D%20%5Calpha%20%2B%20%5Cbeta%20x_i
"\\text{logit}(\\theta_i)  = \\alpha + \\beta x_i")  
This is called a logistic regression model.

## A frequentist approach

To get a rough estimate around where we should plot our grid, we first
compute the maximum likelihood estimate using the standard logistic
regression
tools.

``` r
d.notagg <- data.frame(log.dose=c(rep(-0.86, 5), rep(-0.3,5), rep(-0.05,5), rep(0.73,5) ),
                death=c(rep(0,5), 1, rep(0,4), rep(1, 3), 0, 0, rep(1, 5)) )

fit <- glm(death ~ 1 + log.dose,
          data=d.notagg,
          family="binomial")

summary(fit)
```

``` 

Call:
glm(formula = death ~ 1 + log.dose, family = "binomial", data = d.notagg)

Deviance Residuals: 
     Min        1Q    Median        3Q       Max  
-1.37756  -0.64102  -0.07708   0.05473   1.83495  

Coefficients:
            Estimate Std. Error z value Pr(>|z|)
(Intercept)   0.8466     1.0191   0.831    0.406
log.dose      7.7488     4.8727   1.590    0.112

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 27.526  on 19  degrees of freedom
Residual deviance: 11.789  on 18  degrees of freedom
AIC: 15.789

Number of Fisher Scoring iterations: 7
```

The estimate is ![(\\hat{\\alpha}, \\hat{\\beta}) =
(0.85, 7.75)](https://latex.codecogs.com/png.latex?%28%5Chat%7B%5Calpha%7D%2C%20%5Chat%7B%5Cbeta%7D%29%20%3D%20%280.85%2C%207.75%29
"(\\hat{\\alpha}, \\hat{\\beta}) = (0.85, 7.75)") with standard errors
of 1.0 and 4.9 for
![\\alpha](https://latex.codecogs.com/png.latex?%5Calpha "\\alpha") and
![\\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\\beta"),
respectively.

# Approximating the posterior using a grid

First, we define some functions.

``` r
logit <- function(x) log(x / (1-x) )
invlogit <- function(x) exp(x) / (1 + exp(x))
```

We compute the log posterior using the log likelihood. This helps to
avoid numerical problems. Simplifying some of the expressions, we get
the following function for the log likelihood:

``` r
log.lkhd <- function(alpha, beta) {
  lin <- alpha + beta * d$log.dose
  y <- d$n.deaths; n <- d$n.animals
  sum( y*(lin - log( 1 + exp(lin ) )) + ( n-y )*(-log(1 + exp(lin))) )
}
```

Next, we define a prior function. Since we use a uniform prior, we just
define a constant function:

``` r
prior <- function(alpha, beta) {
  1
}
```

Next, we define our grid and a function to compute the posterior:

``` r
grid_size <- 100
alpha_seq <- seq(-5, 10, length.out=grid_size)
beta_seq <- seq(-10, 40, length.out=grid_size)

alpha_width <- alpha_seq[2] - alpha_seq[1]
beta_width <- beta_seq[2] - beta_seq[1]

post.grid <- expand.grid(alpha = alpha_seq, 
                         beta = beta_seq ) 

posterior.grid <- function(grid, prior_fun=prior) {
  grid %>%
  rowwise %>%
  mutate(loglkhd = log.lkhd(alpha, beta),
         prior = prior_fun(alpha, beta)) %>%
  mutate(log.post = loglkhd + log(prior) ) %>%
  ungroup() %>%
  mutate(log.postm = log.post - max(log.post),
         un.post = exp(log.postm),
         # normalize the posterior
         post = un.post / sum(un.post),
         prior = prior / sum(prior) ) %>%
    select(-log.postm, -un.post)
}
```

Now, we compute the posterior:

``` r
post.grid <- posterior.grid(post.grid)
```

We can plot the posterior density as contour lines. To get the right
contour lines, we use the mode and multiply it with 0.05, 0.1, 0.15, …,
0.95.

``` r
mode <- max(post.grid$post)
breaks <- seq(0.05, 0.95, by=0.1) * mode

unf_post_plot <- post.grid %>%
  ggplot(aes(x=alpha, y=beta, z=post)) + 
  stat_contour(breaks=breaks, col="#377EB8") +
  ylim(-10, 40) +
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8, 10), 
                     limits = c(-5, 10)) +
  theme_minimal() +
  labs(title="Posterior density", subtitle = "with uniform prior")
unf_post_plot
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## Sampling from the posterior

To sample from the posterior, we take the following steps:

1.  Compute the marginal posterior distribution of
    ![\\alpha](https://latex.codecogs.com/png.latex?%5Calpha "\\alpha")
    by numerically summing over
    ![\\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\\beta"):

<!-- end list -->

``` r
marg.alpha <- post.grid %>%
  group_by(alpha) %>%
  summarise(post = sum(post)) %>%
  pull(post)
```

2.  For ![s = 1,
    ..., 1000](https://latex.codecogs.com/png.latex?s%20%3D%201%2C%20...%2C%201000
    "s = 1, ..., 1000") (or more if you want to have more samples),

<!-- end list -->

1)  Draw samples from ![p(\\alpha |
    y)](https://latex.codecogs.com/png.latex?p%28%5Calpha%20%7C%20y%29
    "p(\\alpha | y)"):

<!-- end list -->

``` r
N <- 1000
alpha.sample <- sample(seq(-5, 10, length.out = 100), N, replace=T, prob=marg.alpha)
```

2)  Draw ![\\beta](https://latex.codecogs.com/png.latex?%5Cbeta
    "\\beta") from the discrete conditional distribution ![p(\\beta |
    \\alpha,
    y)](https://latex.codecogs.com/png.latex?p%28%5Cbeta%20%7C%20%5Calpha%2C%20y%29
    "p(\\beta | \\alpha, y)") given the just-sampled value of
    ![\\alpha](https://latex.codecogs.com/png.latex?%5Calpha "\\alpha"):

<!-- end list -->

``` r
beta.sample <- c()
for(i in 1:N) {
  cond.beta <- post.grid %>%
    filter(alpha == alpha.sample[i]) %>%
    mutate(post = post / sum(post)) %>%
    pull(post)
  bsamp <- sample(seq(-10, 40, length.out = 100), 1, prob=cond.beta)
  beta.sample[i] <- bsamp
}
```

3)  For each of the sampled
    ![\\alpha](https://latex.codecogs.com/png.latex?%5Calpha "\\alpha")
    and ![\\beta](https://latex.codecogs.com/png.latex?%5Cbeta
    "\\beta"), add a uniform random jitter, centered at zero with a
    width equal to the spacing of the sampling grid, This gives the
    simulation draws a continuous distribution:

<!-- end list -->

``` r
# add random jitter
alpha.sample <- alpha.sample + runif(N, min = 0 - alpha_width/2,
                                        max = 0 + alpha_width/2)
beta.sample <- beta.sample + runif(N, min = 0 - beta_width/2, 
                                      max = 0 + beta_width/2 )
```

The whole thing as a function:

``` r
extract.sample <- function(density.grid, N=1000, prior=FALSE) {
  if (prior) {
    density.grid <- density.grid %>%
      select(alpha, beta, post=prior)
  } else {
    density.grid <- density.grid %>%
      select(alpha,  beta, post)
  }
  marg.alpha <- density.grid %>%
  group_by(alpha) %>%
  summarise(post = sum(post)) %>%
  pull(post)
  
  alpha.sample <- sample(seq(-5, 10, length.out = 100), N, replace=T, prob=marg.alpha)
  beta.sample <- c()
  for(i in 1:N) {
    cond.beta <- density.grid %>%
      filter(alpha == alpha.sample[i]) %>%
      mutate(post = post / sum(post)) %>%
      pull(post)
    bsamp <- sample(seq(-10, 40, length.out = 100), 1, prob=cond.beta)
    beta.sample[i] <- bsamp
  }
  alpha.sample <- alpha.sample + runif(N, min = 0 - alpha_width/2,
                                        max = 0 + alpha_width/2)
  beta.sample <- beta.sample + runif(N, min = 0 - beta_width/2, 
                                      max = 0 + beta_width/2 )
  
  data.frame(alpha = alpha.sample,
           beta = beta.sample)
}
```

We can now plot the posterior sample:

``` r
post.sample <- extract.sample(post.grid)

unf_post_sample <- post.sample %>%
  ggplot(aes(x=alpha, y=beta)) +
  geom_point(size=0.5) +
  ylim(-10, 40) +
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8, 10), 
                     limits = c(-5, 10)) +
  theme_minimal() +
  labs(title="Posterior sample",
       subtitle="with uniform prior")
unf_post_sample 
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

We can plot the samples as the resulting logistic model together with
the
data:

``` r
plot_samples <- function(sample, n=100, title="Posterior samples", subtitle="") {
  log.dose <- seq(-1.5, 1.5, length.out = 100)

  sample %>%
    sample_n(size=n) %>%
    mutate(id=1:n) %>%
    purrr::pmap_df(~tibble(log.dose=log.dose, id=..3,
                          prop=invlogit(..1 + ..2*log.dose))) %>%
    ggplot(aes(x=log.dose, y=prop)) +
    geom_line(aes(group=id), alpha=0.2, col="#377EB8") +
    geom_hline(yintercept=0.5, linetype='dashed', col="grey") +
    geom_point(data=d, col="#E41A1C") +
    theme_minimal() +
    labs(x="Dose (log g/ml)", y="Proportion of deaths",
         title=title, subtitle = subtitle)
}
plot_samples(post.sample)
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Note that we can also sample and visualize our prior distribution:

``` r
prior.sample <- extract.sample(post.grid, prior=TRUE)

prior.sample %>%
  ggplot(aes(x=alpha, y=beta)) +
  geom_point(size=0.5) +
  ylim(-10, 40) +
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8, 10), 
                     limits = c(-5, 10)) +
  theme_minimal() +
  labs(title="Prior sample")
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

The prior is uniform and thus not very spectacular. For the logistic
model, this then looks as follows:

``` r
plot_samples(prior.sample, n=200, title="Prior samples")
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

## Posterior Distribution of the LD50

We can use the posterior sample to compute the LD50 - the dose level at
which probability of death is 50%. In our logistic model, a 50% survival
rate means   
![\\begin{align\*}
\\text{LD50}: && E(\\frac{y\_i}{n\_i}) = \\text{logit}^{-1}(\\alpha +
\\beta x\_i) = 0.5
\\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%5Ctext%7BLD50%7D%3A%20%26%26%20E%28%5Cfrac%7By_i%7D%7Bn_i%7D%29%20%3D%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28%5Calpha%20%2B%20%5Cbeta%20x_i%29%20%3D%200.5%0A%5Cend%7Balign%2A%7D
"\\begin{align*}
\\text{LD50}: && E(\\frac{y_i}{n_i}) = \\text{logit}^{-1}(\\alpha + \\beta x_i) = 0.5
\\end{align*}")  
Thus ![\\alpha + \\beta x\_i = \\text{logit}(0.5)
= 0](https://latex.codecogs.com/png.latex?%5Calpha%20%2B%20%5Cbeta%20x_i%20%3D%20%5Ctext%7Blogit%7D%280.5%29%20%3D%200
"\\alpha + \\beta x_i = \\text{logit}(0.5) = 0") and the LD50 is ![x\_i
= -\\alpha /
\\beta](https://latex.codecogs.com/png.latex?x_i%20%3D%20-%5Calpha%20%2F%20%5Cbeta
"x_i = -\\alpha / \\beta").

*Attention:* In this example, LD50 is a meaningless concept if ![\\beta
\\leq 0](https://latex.codecogs.com/png.latex?%5Cbeta%20%5Cleq%200
"\\beta \\leq 0"), in which case increasing the dose does not cause the
probability of death to increase.

We report:

1)  The posterior probability that ![\\beta
    \> 0](https://latex.codecogs.com/png.latex?%5Cbeta%20%3E%200
    "\\beta \> 0"), that is, that the drug is harmful:

<!-- end list -->

``` r
mean(post.sample$beta > 0)
```

    [1] 1

From this, we can conclude that the posterior probability of ![\\beta
\> 0](https://latex.codecogs.com/png.latex?%5Cbeta%20%3E%200
"\\beta \> 0") is roughly estimated to exceed 0.999.

2)  The posterior distribution for the LD50 conditional on ![\\beta
    \> 0](https://latex.codecogs.com/png.latex?%5Cbeta%20%3E%200
    "\\beta \> 0"). All draws had positive values of
    ![\\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\\beta"), so
    the distribution is given by the whole sample:

<!-- end list -->

``` r
LD50_samps <- post.sample %>%
  mutate( LD50 = - alpha / beta) 
LD50.mean <- LD50_samps %>%
  summarise(mean = mean(LD50))

unif.LD50.plot <- LD50_samps %>%
  ggplot(aes(x=LD50)) + 
  geom_histogram(bins=50, 
                 fill="#377EB8", col="white") +
  scale_y_continuous(labels = NULL, name="") +
  xlim(-0.6, 0.7) +
  geom_vline(data=LD50.mean, aes(xintercept=mean), col="#E41A1C") +
  theme_minimal() +
  labs(title="Posterior distribution for the LD50",
       subtitle="with uniform prior")
unif.LD50.plot
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

We can incorporate the LD50 data in the plot of the logistic model:

``` r
plot_samples(post.sample) +
  geom_point(data=LD50_samps[1:100,], 
             aes(x=LD50, y=0.5), alpha=0.3, size=0.5)
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

# A different prior

We want to replace the uniform prior density by a joint normal prior
distribution on ![(\\alpha,
\\beta)](https://latex.codecogs.com/png.latex?%28%5Calpha%2C%20%5Cbeta%29
"(\\alpha, \\beta)") with ![\\alpha \\sim
\\text{Normal}(0, 2^2)](https://latex.codecogs.com/png.latex?%5Calpha%20%5Csim%20%5Ctext%7BNormal%7D%280%2C%202%5E2%29
"\\alpha \\sim \\text{Normal}(0, 2^2)"), ![\\beta \\sim
\\text{Normal}(10, 10^2)](https://latex.codecogs.com/png.latex?%5Cbeta%20%5Csim%20%5Ctext%7BNormal%7D%2810%2C%2010%5E2%29
"\\beta \\sim \\text{Normal}(10, 10^2)"), and ![\\text{corr}(\\alpha,
\\beta)=0.5](https://latex.codecogs.com/png.latex?%5Ctext%7Bcorr%7D%28%5Calpha%2C%20%5Cbeta%29%3D0.5
"\\text{corr}(\\alpha, \\beta)=0.5").

``` r
mvn_prior <- function(alpha, beta) {
  rho <- matrix(c(2^2, 2*10*0.5, 2*10*0.5, 10^2), ncol=2)
  dmvnorm(c(alpha, beta), 
           mean=c(0, 10),
           sigma=rho)
}

prior.density <- function(grid, prior_fun){
  grid %>%
    rowwise %>%
    mutate(prior = prior_fun(alpha, beta)) %>%
    ungroup() %>%
    mutate(prior = prior / sum(prior)) 
}

mvn.prior.grid <- prior.density(post.grid, prior_fun=mvn_prior)
```

Our prior density then looks as follows:

``` r
mode <- max(mvn.prior.grid$prior)
breaks <- seq(0.05, 0.95, by=0.1) * mode

mvn.prior.grid %>%
  ggplot(aes(x=alpha, y=beta, z=prior)) + 
  stat_contour(breaks=breaks, col="#E41A1C") +
  ylim(-10, 40) +
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8, 10), 
                     limits = c(-5, 10)) +
  theme_minimal() +
  labs(title="Prior density")
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

We can visualize the samples again:

``` r
prior.sample <- extract.sample(mvn.prior.grid, prior=TRUE)

prior_points <- prior.sample %>%
  ggplot(aes(x=alpha, y=beta)) +
  geom_point(size=0.5) +
  ylim(-10, 40) +
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8, 10), 
                     limits = c(-5, 10)) +
  theme_minimal() +
  labs(title="Prior sample", subtitle="with multivariate normal prior")

prior_model <- plot_samples(prior.sample, n=200, 
                            title="Prior sample", subtitle="with multivariate normal prior")
plot_grid(prior_points, prior_model)
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

We can see that the prior still allows a wide range of different models
but different to the uniform prior, it is much more restricted to a
certain range that is already very close to the observed data.

## The new posterior density

We now use this prior to compute our posterior denstiy.

``` r
mvn.post.grid <- posterior.grid(post.grid, prior_fun=mvn_prior)
```

We can compare our new posterior density with the old posterior density
(obtained with a uniform prior)

``` r
mode <- max(mvn.post.grid$post)
breaks <- seq(0.05, 0.95, by=0.1) * mode

mvn_post_plot <- mvn.post.grid %>%
  ggplot(aes(x=alpha, y=beta, z=post)) + 
  stat_contour(breaks=breaks, col="#377EB8") +
  ylim(-10, 40) +
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8, 10), 
                     limits = c(-5, 10)) +
  scale_color_brewer(palette = "Set1", name="",
                     label=c("post"="posterior", "prior"), direction = -1) +
  theme_minimal() +
  labs(title="Posterior density", subtitle = "with multivariate normal prior")

plot_grid(unf_post_plot, mvn_post_plot)
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

We can see that our new prior is slightly more regularizing than the
uniform prior: The new posterior density is a bit tighter than the old
posterior.

We can also compare the maximum a posteriori estimates: With the uniform
prior, we have as MAP estimate:

``` r
post.grid[which.max(post.grid$post),]
```

    # A tibble: 1 x 6
      alpha  beta loglkhd  prior log.post    post
      <dbl> <dbl>   <dbl>  <dbl>    <dbl>   <dbl>
    1 0.909  8.18   -5.90 0.0001    -5.90 0.00322

and with the multivariate normal posterior, we get:

``` r
mvn.post.grid[which.max(mvn.post.grid$post),]
```

    # A tibble: 1 x 6
      alpha  beta loglkhd    prior log.post    post
      <dbl> <dbl>   <dbl>    <dbl>    <dbl>   <dbl>
    1 0.758  8.18   -5.92 0.000613    -10.8 0.00403

While the beta value is the same (at least in our grid approximation),
the new posterior alpha is a bit closer to zero.

## Sampling from the new posterior

``` r
mvn.post.sample <- extract.sample(mvn.post.grid)
mvn.prior.sample <- extract.sample(mvn.prior.grid, prior = T)

post_plot <- mvn.post.sample %>%
  ggplot(aes(x=alpha, y=beta)) +
  geom_point(size=0.5) +
  ylim(-10, 40) +
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8, 10), 
                     limits = c(-5, 10)) +
  theme_minimal() +
  labs(title="Posterior sample", subtitle = "with multivariate normal prior")

prior_plot <- mvn.prior.sample %>%
  ggplot(aes(x=alpha, y=beta)) +
  geom_point(size=0.5) +
  ylim(-10, 40) +
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8, 10), 
                     limits = c(-5, 10)) +
  theme_minimal() +
  labs(title="Prior sample", subtitle = "with multivariate normal prior")

plot_grid(post_plot, prior_plot)
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

We again use the sample to visualize the posterior for the logistic
model:

``` r
plot_samples(mvn.post.sample) 
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

We can compare this again to the prior samples which also shows how the
posterior is a compromise between the data and the prior.

``` r
plot_samples(mvn.prior.sample) +
  labs(title="Prior samples")
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

## New LD50

``` r
mvn.LD50_samps <- mvn.post.sample %>%
  mutate( LD50 = - alpha / beta) 

mvn.LD50.mean <- mvn.LD50_samps %>%
  summarise(mean = mean(LD50))

mvn.LD50.plot <- mvn.LD50_samps %>%
  ggplot(aes(x=LD50)) + 
  geom_histogram(bins=50, 
                 fill="#377EB8", col="white") +
  scale_y_continuous(labels = NULL, name="") +
  xlim(-0.6, 0.7) + 
  geom_vline(data=mvn.LD50.mean, aes(xintercept=mean), col="#E41A1C") +
  theme_minimal() +
  labs(title="Posterior distribution for the LD50", 
       subtitle="with multivariate normal prior")

plot_grid(mvn.LD50.plot, unif.LD50.plot)
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

# Normal approximation of the posterior

The fourth chapter *Asymptotics and Non-Bayesian Approaches* explains
how it is possible to approximate the posterior distribution, using the
mode and a normal distribution. We will now compute the normal
approximation and compare it to the exact posterior obtained by the
uniform prior.

Since we assume a uniform prior density for ![(\\alpha,
\\beta)](https://latex.codecogs.com/png.latex?%28%5Calpha%2C%20%5Cbeta%29
"(\\alpha, \\beta)"), the posterior mode is the same as the maximum
likelihood estimate. So we get the mode by computing the MLE:

``` r
# lkhd function
bioassayfun <- function(w, df) {
  z <- w[1] + w[2]*df$log.dose
  -sum(df$n.deaths*(z) - df$n.animals*log1p(exp(z)))
}

#' Optimize
w0 <- c(0,0)
optim_res <- optim(w0, bioassayfun, gr = NULL, d, hessian = T)
# w is the mode
w <- optim_res$par
# this computes the inverse of the hessian at the mode
S <- solve(optim_res$hessian)

#' Multivariate normal probability density function
dmvnorm <- function(x, mu, sig)
  exp(-0.5*(length(x)*log(2*pi) + log(det(sig)) + (x-mu) %*% solve(sig, x - mu)))

#' Evaluate likelihood at points (alpha, beta) 
ab_grid <- expand.grid(alpha = alpha_seq, 
                         beta = beta_seq ) 

ab_grid$lkhd <- apply(ab_grid, 1, dmvnorm, w, S)




#' Create a plot of the posterior density
norm_post_plot <- ggplot(data = ab_grid, aes(x = alpha, y = beta, z=lkhd)) +
  stat_contour( col="#377EB8") +
  ylim(-10, 40) +
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8, 10), 
                     limits = c(-5, 10)) +
  theme_minimal() +
  labs(x = 'alpha', y = 'beta',
       title="Posterior Density",
       subtitle="using normal Approximation" ) 

norm_post_plot
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

The posterior density is very similar to the one we obtained before, but
it is missing the slight skew in the upper corner.

Similarly for the posterior sample:

``` r
# sample from the multivariate model
norm_sample <- MASS::mvrnorm(N, w, S) %>%
  data.frame() %>%
  rename(alpha=X1, beta=X2)

norm_post_sample <- norm_sample %>%
  ggplot() +
  geom_point(aes(alpha, beta), size=0.5) +
  ylim(-10, 40) +
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8, 10), 
                     limits = c(-5, 10)) +
  theme_minimal() +
  labs(title="Posterior sample",
       subtitle="using normal approximation")
norm_post_sample
```

    Warning: Removed 1 rows containing missing values (geom_point).

![](bioassayExample_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

This also affects the resulting LD50 and the probability that ![\\beta
\> 0](https://latex.codecogs.com/png.latex?%5Cbeta%20%3E%200
"\\beta \> 0"):

``` r
mean(norm_sample$beta > 0)
```

    [1] 0.938

Before, this probability was ~1.

``` r
LD50_norm_samps <- norm_sample %>%
  filter( beta > 0 ) %>%
  mutate( LD50 = - alpha/beta )

LD50_norm_mean <- LD50_norm_samps %>%
  summarise(mean = mean(LD50))

norm.LD50.plot <- LD50_norm_samps %>%
  ggplot(aes(x=LD50)) + 
  geom_histogram(bins=50, 
                 fill="#377EB8", col="white") +
  scale_y_continuous(labels = NULL, name="") +
  xlim(-1, 1) +
  geom_vline(data=LD50_norm_mean, aes(xintercept=mean), col="#E41A1C") +
  theme_minimal() +
  labs(title="Posterior distribution for the LD50",
       subtitle="using normal approximation")

norm.LD50.plot
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

## Comparison

A direct comparison of the normal approximation with the exact posterior
makes the differences clearer:

``` r
#' Combine the plots
plot_grid(unf_post_plot, unf_post_sample, 
             unif.LD50.plot + xlim(-1, 1), norm_post_plot, 
             norm_post_sample, norm.LD50.plot, ncol = 3)
```

![](bioassayExample_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->
