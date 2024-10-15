# mermboost

## Installation Instructions

-   Currently only available via github:

    ``` r
    devtools::install_github("larsknieper/mermboost")
    library("mermboost")
    ```

## Tutorial

Simulations as well as data applications show that users can benefit from using *mermboost* as unbiased estimates as well as an enabled estimation of random effects variation is included.

Technically, it is calling a newly defined family (*merMod*), in which random components are treated as nuisance parameters and are repeatedly estimated by an *glmer()* object from *lme4*. Simultaneously, fixed effects are estimated in a common way by the established *mboost* package. Since a *lme4* model is included *mermboost* is restricted to distributions of base-*R*-families with the negative binomial distribution in addition. From a user's perspective, four functions are provided, namely:

| Function       | Description                 | Class               |
|----------------|-----------------------------|---------------------|
| `find_ccc()`   | Check for cluster constants |                     |
| `glmerboost()` | For estimating a GLMM       | *glmermboost* class |
| `mermboost()`  | For estimating a GAMM       | *mermboost* class   |
| `mer_cvrisk()` | For individual-sensitive    | *mer_cv* class      |
|                | Cross-validation            |                     |

The function *find\\\_ccc()* takes the arguments of a dataframe and a name of a cluster-defining variable. It gives logical values corresponding to all further variables, whether or not these are constant across each cluster as this is causing the described bias. Still, by providing an estimation of random effects' covariance $\boldsymbol{Q}$ *mermboost* is a good choice for boosting mixed models even when there are no cluster constant covariates in the data set.

The boosting functions are wrappers of the corresponding *mboost* functions. Hence, *glmerboost()* is a wrapper of *glmboost()* and *mermboost()* is a wrapper of *mboost()*. After initialising a *lme4* object, both start the corresponding *mboost* model with the *merMod* family. For the GLMM (*glmerboost*) the intercept is estimated with the *lme4* object and hence, needs a correction of $$
\beta_0 = \bar{\tilde{\boldsymbol{x}}}_p \boldsymbol{\beta}_{-0},
$$ where $\bar{\tilde{\boldsymbol{x}}}_p$ denotes the means of the picked slope covariates and $\hat{\boldsymbol{\beta}}_{-0}$ are corresponding estimated slope coefficients. This corrected intercept is stored as an extra element in the *glmer* class object as *model\\\$intercept*.

To call the boosting functions some modifications are needed in comparison to *mboost*. The random effects are not specified within the base-learner *brandom()* anymore but in the way they are specified in *lme4* (and other mixed model packages) by "*\... + (random formula \| id)*". By this, the add-on package differs visually from the *mboost* package even though this leads to a mixed formula format. Thus, to estimate a specified model formulation of $\boldsymbol{\eta} = \beta_0 + \beta_1 \boldsymbol{x}_1 + \beta_2 \boldsymbol{x}_2 + \boldsymbol{\gamma}_0$ for a binomial distributed response variable $\boldsymbol{y}$ with a cluster specifying *id* the differences are as follows:

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

To clarify this even further, assume a GAMM with an additional random slope defined as $\boldsymbol{\eta} = f_1 (\boldsymbol{x}_1) + f_2 (\boldsymbol{x}_2) + \boldsymbol{\gamma}_0 + \boldsymbol{\gamma}_1 \boldsymbol{z}_2$, where $\boldsymbol{x}_2 = \boldsymbol{z}_2$. The corresponding models are to be called as follows without consideration of possible adjustments within *bbs()*:

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

Some methods of functions got adjusted to, where *predict()* is considering the nuisance estimations and therefore the random effects as well. With the argument *RE* it can be controlled whether to include random effects (*TRUE*) or not (*FALSE*).

To perform a cross-validation with a cluster-wise test-train split *mer\\\_cvrisk()* takes the *mermboost* object and the number of folds as an argument. The folds are then generated automatically within the function considering the cluster structure by its identifying variable, which is stored in the *mermboost* object. However, this might lead to unwanted splits. This is why alternatively pre-specified folds can be given as an argument as well by a matrix containing 0's and 1's indicating the folds. The cross-validation can can further be conducted by parallel computing, where the number of cores are to be specified by the *core* argument.
