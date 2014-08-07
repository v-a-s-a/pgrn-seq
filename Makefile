## testing out GNU make to compile knitr docs
R_OPTS=--vanilla
lab_notebook/analysis-description.pdf: lab_notebook/analysis-description.Rmd
	cd lab_notebook; R ${R_OPTS} -e 'getwd();library(knitr);knit("analysis-description.Rmd")'
	cd lab_notebook; pandoc analysis-description.md -o analysis-description.pdf

clean:
	cd lab_notebook; rm analysis-description.md
