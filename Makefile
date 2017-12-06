all: ./results/Figures/Figure3D.pdf ./results/Figures/Figure1A.pdf ./results/Figures/Figure3C.pdf
.PHONY : all
#results/results.table.tsv:Kimura_diffusion.ipnb data/processed/secondary/Intrahost_initially_present.csv


#./data/processed/secondary/simulated_fits.csv :  

#./data/processed/secondary/one_per_person.csv ./data/processed/seconday/removed_data.csv : ./data/processed/secondary/Intrahost_initially_present.csv ./results/discrete_Ne.ipynb


#./results/Figures/Figure2D.pdf ./results/Figures/Figure2C.pdf ./Figures/Figure2B.pdf ./results/Figures/Figure2A.pdf ./Figures/Supplemental_Figure6.pdf : ./results/intrahost_model.Rmd ./data/processed/secondary/one_per_person.csv ./data/processed/secondary/Intrahost_initially_present.csv ./data/processed/secondary/simulated_fits.csv ./data/processed/seconday/removed_data.csv


./results/Figures/Figure3D.pdf  ./results/Figures/Figure3F.pdf ./results/Figures/Supplemental_Figure7D.pdf:Figure3DF
	# Just another empty recipe.

.INTERMEDIATE: Figure3DF
Figure3DF: ./data/processed/secondary/transmission_pairs_freq.poly.donor.csv ./data/processed/secondary/no_cut_transmission_pairs_freq.poly.donor.csv ./results/transmission_models.Rmd
	cd ./results; Rscript -e 'library(rmarkdown); rmarkdown::render("transmission_models.Rmd", "github_document")'

./results/Figures/Figure3C.pdf ./results/Figures/Figure3E.pdf ./results/Supplemental_Figure7C.pdf: Figure3CE_sup7
	# Just another empty recipe

.INTERMEDIATE: Figure3CE_sup7
Figure3CE_sup7: ./data/reference/all_meta.sequence_success.csv ./data/processed/secondary/qual.snv.csv ./data/processed/secondary/transmission_pairs.csv ./data/processed/secondary/trans_freq.csv ./data/processed/secondary/transmission_pairs_freq.poly.donor.csv ./data/processed/secondary/no_cut_trans_freq.csv ./data/processed/secondary/no_cut_transmission_pairs_freq.poly.donor.csv ./data/processed/secondary/community_pairs.freq.csv ./data/processed/secondary/community_pairs_freq.poly.donor.csv ./results/transmission_setup.Rmd
	cd ./results; Rscript -e 'library(rmarkdown); rmarkdown::render("transmission_setup.Rmd", "github_document")' 

./results/Figures/Figure1A.pdf ./results/Figures/Figure1B.pdf ./results/Figures/Figure1C.pdf ./results/Figures/Figure1D2.pdf\
./results/Figures/Figure3A.pdf ./results/Figures/Figure3B.pdf\
./results/Figures/Supplemental_Figure4.pdf: Figure1_3AB
	#Empty recipe to propogate newness

.INTERMEDIATE: Figure1_3AB
Figure1_3AB:./data/reference/all_meta.sequence_success.csv ./data/processed/secondary/meta_one.sequence.success.csv ./data/processed/secondary/transmission_pairs.csv ./data/processed/secondary/qual.snv.csv ./data/processed/secondary/possible.pairs.dist.csv ./results/Summary_stats.Rmd
	cd ./results; Rscript -e 'library(rmarkdown); rmarkdown::render("Summary_stats.Rmd", "github_document")'

./data/processed/secondary/community_pairs.freq.csv ./data/processed/secondary/community_pairs_freq.poly.donor.csv ./data/processed/secondary/no_cut_transmission_pairs_freq.poly.donor.csv ./data/processed/secondary/no_cut_trans_freq.csv ./data/processed/secondary/transmission_pairs_freq.poly.donor.csv ./data/processed/secondary/trans_freq.csv : transmission.setup.intermediate
	#Empty recipe to propagate "newness" from the intermediate to final targets
.INTERMEDIATE: transmission.setup.intermediate
transmission.setup.intermediate: ./data/processed/secondary/qual.snv.csv ./data/processed/secondary/no_freq_cut.qual.snv.csv ./data/processed/secondary/transmission_pairs.csv
	cd ./scripts/R; Rscript settingup_transmission.R


./data/processed/secondary/Intrahost_initially_present.csv ./data/processed/secondary/Intrahost_pairs.csv ./data/processed/secondary/Intrahost_all.csv: intrahost.setup.intermediate
	# Empty recipe to propagate "newness" from the intermediate to final targets
.INTERMEDIATE: intrahost.setup.intermediate
intrahost.setup.intermediate: ./data/processed/secondary/qual.snv.csv ./data/reference/all_meta.sequence_success.csv ./scripts/R/Intrahost_setup.R
	cd ./scripts/R; Rscript Intrahost_setup.R

./data/processed/secondary/possible.pairs.dist.csv ./data/processed/secondary/transmission_pairs.csv: distance.intermediate
	# Empty recipe to propagate "newness" from the intermediate to final targets

.INTERMEDIATE: distance.intermediate
distance.intermediate: ./data/processed/secondary/qual.snv.csv ./scripts/R/L1_norm.R
	cd ./scripts/R; Rscript L1_norm.R


./data/processed/secondary/qual.snv.csv ./data/processed/secondary/no_freq_cut.qual.snv.csv : qual.intermediate
	# Empty recipe to propagate "newness" from the intermediate to final targets
	
.INTERMEDIATE: qual.intermediate
qual.intermediate: ./data/processed/secondary/average_coverages.csv ./scripts/R/processing_snv.R
	cd ./scripts/R; Rscript processing_snv.R

./data/processed/secondary/average_coverages.csv: ./scripts/R/processing_coverage.R
	cd ./scripts/R; Rscript processing_coverage.R

