## testing out GNU make to compile knitr docs
R_OPTS=--vanilla

all: reports/analysis-description.pdf reports/phenotype_preprocessing.pdf reports/single_marker_analysis.pdf

reports/analysis-description.pdf: lab_notebook/analysis-description.Rmd
	cd lab_notebook; R ${R_OPTS} -e 'library(knitr);knit("analysis-description.Rmd")'
	cd lab_notebook; pandoc analysis-description.md -o ../reports/analysis-description.pdf

reports/phenotype_preprocessing.pdf: analyses/phenotype_preprocessing.Rmd
	cd analyses/; R ${R_OPTS} -e 'library(knitr);knit("phenotype_preprocessing.Rmd")'
	cd analyses/; pandoc phenotype_preprocessing.md -o ../reports/phenotype_preprocessing.pdf 

reports/anc_single_marker.pdf: analyses/anc_single_marker.Rmd
	cd analyses/; R ${R_OPTS} -e 'library(knitr);knit("anc_single_marker.Rmd")'
	cd analyses/; pandoc anc_single_marker.md -o ../reports/anc_single_marker.pdf 

reports/loganc_single_marker.pdf: analyses/loganc_single_marker.Rmd
	cd analyses/; R ${R_OPTS} -e 'library(knitr);knit("loganc_single_marker.Rmd")'
	cd analyses/; pandoc loganc_single_marker.md -o ../reports/loganc_single_marker.pdf 

reports/invnormanc_single_marker.pdf: analyses/invnormanc_single_marker.Rmd
	cd analyses/; R ${R_OPTS} -e 'library(knitr);knit("invnormanc_single_marker.Rmd")'
	cd analyses/; pandoc invnormanc_single_marker.md -o ../reports/invnormanc_single_marker.pdf 

## I refactored the single marker analysis into separate reports for each phenotype transform
#reports/single_marker_analysis.pdf: analyses/single_marker_analysis.Rmd
#	cd analyses/; R ${R_OPTS} -e 'library(knitr);knit("single_marker_analysis.Rmd")'
#	cd analyses/; pandoc single_marker_analysis.md -o ../reports/single_marker_analysis.pdf 

reports/anc_skat-o.pdf: analyses/anc_skat-o.Rmd
	cd analyses/;  R ${R_OPTS} -e 'library(knitr);knit("anc_skat-o.Rmd")'
	cd analyses/; pandoc anc_skat-o.md -o ../reports/anc_skat-o.pdf

reports/loganc_skat-o.pdf: analyses/loganc_skat-o.Rmd
	cd analyses/;  R ${R_OPTS} -e 'library(knitr);knit("loganc_skat-o.Rmd")'
	cd analyses/; pandoc loganc_skat-o.md -o ../reports/loganc_skat-o.pdf

reports/invnormanc_skat-o.pdf: analyses/invnormanc_skat-o.Rmd
	cd analyses/;  R ${R_OPTS} -e 'library(knitr);knit("invnormanc_skat-o.Rmd")'
	cd analyses/; pandoc invnormanc_skat-o.md -o ../reports/invnormanc_skat-o.pdf

## refactored gene-based analysis into separate reports
#reports/gene_based_analysis.pdf: analyses/gene_based_analysis.Rmd
#	cd analyses/;  R ${R_OPTS} -e 'library(knitr);knit("gene_based_analysis.Rmd")'
#	cd analyses/; pandoc gene_based_analysis.md -o ../reports/gene_based_analysis.pdf

clean:
	cd lab_notebook; rm analysis-description.md
	cd analyses/; rm phenotype_preprocessing.md anc_single_marker.md gene_based_analysis.md
	cd analyses/; rm loganc_single_marker.md invnorm_single_marker.md
