# BDA project: Parameter fitting for damped harmonic oscillator using Stan

This repo contains the course project done for *CS-E5710 Bayesian Data Analysis* course in Aalto university fall 2022. The idea of the project was to invert system parameters of a damped harmonic oscillator using bayesian workflow and `Stan` probabilistic programming language. The work was conducted by the authors listed in the contributors of the repository. The `.rmd` file has all the code while the `.pdf` files have the actual report and the presentation. The `.Rda` files have the fit ocjects, so that the code can be run without running the `Stan` models as they are computationally heavy. 

If you want to run the `Stan` models, change the headers of the `knitr` code chunks from
`{cmdstan, output.var = "model_rk45_...", eval=FALSE}`
to
`{cmdstan, output.var = "model_rk45_...", eval=TRUE}`
and the `run_stan` varibale to `TRUE`. It's best not to try to knit the document in this setting as it doesn't work well, for compiling the pdf-document use the ready `.Rda` fit objects.
