#Rscript -e 'library(rmarkdown); rmarkdown::render("./Summary_stats.Rmd", "github_document")'
#Rscript -e 'library(rmarkdown); rmarkdown::render("./transmission_setup.Rmd", "github_document")'
#Rscript -e 'library(rmarkdown); rmarkdown::render("./transmission_models.Rmd", "github_document")'
#Rscript -e 'library(rmarkdown); rmarkdown::render("./WFABC.Rmd", "github_document")'
#python Kimura_diffusion.py
#julia discrete_Ne.jl
Rscript -e 'library(rmarkdown); rmarkdown::render("./intrahost_model.Rmd", "github_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("./mutational_model.Rmd", "github_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("./supplemental.Rmd", "github_document")'
