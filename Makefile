## testing out GNU make to compile knitr docs
R_OPTS=--vanilla
lab_notebook/analysis-description.pdf: lab_notebook/analysis-description.Rmd
	cd lab_notebook; R ${R_OPTS} -e 'library(knitr);knit("analysis-description.Rmd")'
	cd lab_notebook; pandoc analysis-description.md -o analysis-description.pdf

analyses/phenotype_preprocessing.pdf: analyses/phenotype_preprocessing.Rmd
	cd analyses/; R ${R_OPTS} -e 'library(knitr);knit("phenotype_preprocessing.Rmd")'
	cd analyses/; pandoc phenotype_preprocessing.md -o phenotype_preprocessing.pdf 

analyses/single_marker_analysis.pdf: analyses/single_marker_analysis.Rmd
	cd analyses/; R ${R_OPTS} -e 'library(knitr);knit("single_marker_analysis.Rmd")'
	cd analyses/; pandoc single_marker_analysis.md -o single_marker_analysis.pdf 

clean:
	cd lab_notebook; rm analysis-description.md
	cd analyses/; rm phenotype_preprocessing.md single_marker_analysis.md
