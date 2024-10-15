# mermboost

## Installation Instructions

-   Currently only available via github:

    ``` r
    devtools::install_github("boost-R/mboost")
    library("mboost")
    ```

Simulations as well as data applications show that users can benefit from using \textit{mermboost} as unbiased estimates as well as an enabled estimation of random effects variation is included.

Technically, it is calling a newly defined family (\textit{merMod}), in which random components are treated as nuisance parameters and are repeatedly estimated by an \texttt{glmer()} object from \textit{lme4}. Simultaneously, fixed effects are estimated in a common way by the established \textit{mboost} package. Since a \textit{lme4} model is included \textit{mermboost} is restricted to distributions of base-\texttt{R}-families with the negative binomial distribution in addition. From a user's perspective, four functions are provided, namely:

| Function       | Description                 | Class               |
|----------------|-----------------------------|---------------------|
| `find_ccc()`   | Check for cluster constants |                     |
| `glmerboost()` | For estimating a GLMM       | *glmermboost* class |
| `mermboost()`  | For estimating a GAMM       | *mermboost* class   |
| `mer_cvrisk()` | For individual-sensitive    | *mer_cv* class      |
|                | Cross-validation            |                     |

The function \texttt{find\_ccc()} takes the arguments of a dataframe and a name of a cluster-defining variable. It gives logical values corresponding to all further variables, whether or not these are constant across each cluster as this is causing the described bias. Still, by providing an estimation of random effects' covariance $\boldsymbol{Q}$ \textit{mermboost} is a good choice for boosting mixed models even when there are no cluster constant covariates in the data set.

The boosting functions are wrappers of the corresponding \textit{mboost} functions. Hence, \texttt{glmerboost()} is a wrapper of \texttt{glmboost()} and \texttt{mermboost()} is a wrapper of \texttt{mboost()}. After initialising a \textit{lme4} object, both start the corresponding \textit{mboost} model with the \textit{merMod} family. For the GLMM (\texttt{glmerboost}) the intercept is estimated with the \textit{lme4} object and hecen, needs a correction of $$\beta_0 = \bar{\tilde{\boldsymbol{x}}}_p \boldsymbol{\beta}_{-0},$$ where $\bar{\tilde{\boldsymbol{x}}}_p$ denotes the means of the picked slope covariates and $\hat{\boldsymbol{\beta}}_{-0}$ are corresponding estimated slope coefficients. This corrected intercept is stored as an extra element in the \textit{glmer} class object as \texttt{model\$intercept}.

To call the boosting functions some modifications are needed in comparison to \textit{mboost}. The random effects are not specified within the base-learner \texttt{brandom()} anymore but in the way they are specified in \textit{lme4} (and other mixed model packages) by "\texttt{... + (random formula | id)}". By this, the add-on package differs visually from the \textit{mboost} package even though this leads to a mixed formula format. Thus, to estimate a specified model formulation of $\boldsymbol{\eta} = \beta_0 + \beta_1 \boldsymbol{x}_1 + \beta_2 \boldsymbol{x}_2 + \boldsymbol{\gamma}_0$ for a binomial distributed response variable $\boldsymbol{y}$ with a cluster specifying \textit{id} the differences are as follows:

### mboost

``` r
mboost(y ~ bols(x1) + bols(x2) + brandom(id),
        family = Binomial(type = "glm"),
        data = df)
```

### mermboost

``` r
glmerboost(y ~ x1 + x2 + (1 | id),,
            family = binomial(),
            data = df)
```

To clarify this even further, assume a GAMM with an additional random slope defined as $\boldsymbol{\eta} = f_1 (\boldsymbol{x}_1) + f_2 (\boldsymbol{x}_2) + \boldsymbol{\gamma}_0 + \boldsymbol{\gamma}_1 \boldsymbol{z}_2$, where $\boldsymbol{x}_2 = \boldsymbol{z}_2$. The corresponding models are to be called as follows without consideration of possible adjustments within \texttt{bbs()}:

### mboost

``` r
mboost(y ~ bbs(x1) + bbs(x2) + brandom(id) + brandom(id, by = x2),
        family = Binomial(type = "glm"), data = df)
```

### mermboost

``` r
mermboost(y ~ bbs(x1) + bbs(x2) + (1 + x2 | id), 
          family = binomial(), data = df)
```

Some methods of functions got adjusted to, where \texttt{predict()} is considering the nuisance estimations and therefore the random effects as well. With the argument \texttt{RE} it can be controlled whether to include random effects (\texttt{TRUE}) or not (\texttt{FALSE}).

To perform a cross-validation with a cluster-wise test-train split \texttt{mer\_cvrisk()} takes the \textit{mermboost} object and the number of folds as an argument. The folds are then generated automatically within the function considering the cluster structure by its identifying variable, which is stored in the \textit{mermboost} object. However, this might lead to unwanted splits. This is why alternatively pre-specified folds can be given as an argument as well by a matrix containing 0's and 1's indicating the folds. The cross-validation can can further be conducted by parallel computing, where the number of cores are to be specified by the \texttt{core} argument.
