---
title: "Chapter 2 - Exercise 12"
author: "Aaron McDaid - aaron.mcdaid@gmail.com"
date: "2 May 2019"
output:
    html_document:
        keep_md: true
---

<div style="display:none">
  $$
  \newcommand{\EE}[1]{\mathbb{E}\mathopen{}\left[ #1 \right]\mathclose{}}
  \newcommand{\Var}[1]{\mathrm{Var}\mathopen{}\left[ #1 \right]\mathclose{}}
  $$
</div>



## Chapter 2 Exercise 9

_12. Jeffreys’ prior distributions: suppose y|θ ∼ Poisson(θ). Find Jeffreys’ prior density for θ, and then find α and β for which the Gamma(α, β) density is a close match to Jeffreys’ density._

pdf for Poisson:

$$ f(y ; \theta) = \frac{\theta^y e^{-\theta}}{y!} $$

The Jeffreys prior density is proportional to:

\begin{align}
    p(\theta) \propto \sqrt{I(\theta)} = \sqrt{\EE{\left( \frac{d}{d\theta} \log f(y; \theta) \right)^2 \middle| \theta}}
\end{align}

Breaking this into steps, we start with the density $f(y; \theta)$

$$ \frac{\theta^y e^{-\theta}}{y!} $$

apply the logarithm

\begin{align}
\log(\cdot) & = y \log(\theta) - \theta \log(e) - \log (y!)
\\          & = y \log(\theta) - \theta - \log (y!)
\end{align}

then the derivative

$$
\frac{d}{d\theta}(\cdot) = \frac{y}{\theta} - 1
$$

Next, we need the expectation of the square, i.e.

$$ \EE{\left( \frac{y}{\theta} - 1 \right)^2 \middle| \theta} $$

Using the fact that $\Var{X} = \EE{X^2} - \EE{X}^2$, we can compute the expectation of the square as $\Var{X} + \EE{X}^2$, i.e:

\begin{align}
\EE{\left( \frac{y}{\theta} - 1 \right)^2 \middle| \theta} & = \Var{ \frac{y}{\theta} - 1  \middle| \theta} + \EE{\frac{y}{\theta} - 1 \middle| \theta}^2
\\                                         & = \frac1{\theta^2} \Var{ y \middle| \theta} + \left(\frac1{\theta}\EE{y \middle| \theta} - 1\right)^2
\end{align}

Then, because $y$ has a Poisson distribution, the mean and variance of $y$ are simply $\Var{y|\theta} = \EE{y|\theta} = \theta$,
\begin{align}
\EE{\left( \frac{y}{\theta} - 1 \right)^2 \middle| \theta} & = \frac1{\theta^2} \theta + \left(\frac1{\theta} \theta - 1\right)^2
\\                                                         & = \frac1{\theta} + \left(1 - 1\right)^2
\\                                                         & = \frac1{\theta}
\end{align}

Final step in the Jeffreys formula above is to apply the square root, giving us the Jeffreys' prior for Poisson:

$$ p(\theta) \propto \sqrt{I(\theta)} = \sqrt{\frac1{\theta}} = \theta^{-\frac12} $$


To complete the exercise, we are also asked to find a parameterization of the Gamma distribution which is a close match to this:

$$ \theta \sim Gamma(\alpha, \beta) $$

$$ f(\theta | \alpha, \beta) = \frac{\beta^\alpha}{\Gamma(\alpha)} \theta^{\alpha-1} e^{-\beta\theta} $$

This can be achieved by setting $\alpha=\frac12$ and $\beta=0$

\begin{align}
f(\theta | \alpha, \beta) & = \frac{0^\frac12}{\Gamma\left(\frac12\right)} \theta^{\frac12-1} e^{-0\theta}
\\                          = \frac{0^\frac12}{\Gamma\left(\frac12\right)} \theta^{-\frac12}
\end{align}

As required, that is proportional to $\theta^{-\frac12}$. But it's an [improper prior](https://en.wikipedia.org/wiki/Prior_probability#Examples).
