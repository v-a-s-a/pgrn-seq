## testing out GNU make to compile knitr docs
R_OPTS=--vanilla

all: reports/analysis-description.pdf reports/phenotype_preprocessing.pdf reports/single_marker_analysis.pdf

reports/analysis-description.pdf: lab_notebook/analysis-description.Rmd
	cd lab_notebook; R ${R_OPTS} -e 'library(knitr);knit("analysis-description.Rmd")'
	cd lab_notebook; pandoc analysis-description.md -o ../reports/analysis-description.pdf

reports/phenotype_preprocessing.pdf: analyses/phenotype_preprocessing.Rmd
	cd analyses/; R ${R_OPTS} -e 'library(knitr);knit("phenotype_preprocessing.Rmd")'
	cd analyses/; pandoc phenotype_preprocessing.md -o ../reports/phenotype_preprocessing.pdf 

reports/single_marker_analysis.pdf: analyses/single_marker_analysis.Rmd
	cd analyses/; R ${R_OPTS} -e 'library(knitr);knit("single_marker_analysis.Rmd")'
	cd analyses/; pandoc single_marker_analysis.md -o ../reports/single_marker_analysis.pdf 

reports/gene_based_analysis.pdf: analyses/gene_based_analysis.Rmd
	cd analyses/;  R ${R_OPTS} -e 'library(knitr);knit("gene_based_analysis.Rmd")'
	cd analyses/; pandoc gene_based_analysis.md -o ../reports/gene_based_analysis.pdf

clean:
	cd lab_notebook; rm analysis-description.md
	cd analyses/; rm phenotype_preprocessing.md single_marker_analysis.md gene_based_analysis.md
