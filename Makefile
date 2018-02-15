# --------------------- Secondary Figures ---------------------

Figures : ./results/Figures/Figure1A.pdf ./results/Figures/Figure1B.pdf \
./results/Figures/Figure1C.pdf ./results/Figures/Figure1D.pdf \
./results/Figures/Figure2A.pdf ./results/Figures/Figure2B.pdf \
./results/Figures/Figure2C.pdf ./results/Figures/Figure2D.pdf \
./results/Figures/Figure2E.pdf ./results/Figures/Figure3A.pdf \
./results/Figures/Figure3B.pdf ./results/Figures/Figure3C.pdf \
./results/Figures/Figure3D.pdf ./results/Figures/Figure4A.pdf ./results/Figures/Figure4B.pdf
	echo "Figures made"
####################### Figure 1 ######################

./results/Figures/Figure1A.pdf ./results/Figures/Figure1B.pdf \
./results/Figures/Figure1C.pdf ./results/Figures/Figure1D.pdf : Figure1.intermediate
	#Empty recipe to propagate "newness" from the intermediate to final targets

.INTERMEDIATE: Figure1.intermediate
Figure1.intermediate: ./data/reference/all_meta.sequence_success.csv ./data/processed/secondary/qual.snv.csv \
./data/reference/segs.csv ./scripts/secondary/Figures/Figure1.R
	Rscript ./scripts/secondary/Figures/Figure1.R

####################### Figure 2 ######################

./results/Figures/Figure2A.pdf ./results/Figures/Figure2B.pdf \
./results/Figures/Figure2C.pdf ./results/Figures/Figure2D.pdf ./results/Figures/Figure2E.pdf : Figure2.intermediate
        #Empty recipe to propagate "newness" from the intermediate to final targets
.INTERMEDIATE: Figure2.intermediate
Figure2.intermediate: ./data/processed/secondary/qual.snv.csv ./data/processed/secondary/antigenic_isnv.csv \
./data/processed/secondary/global_freq_antigenic.tsv ./data/processed/secondary/minor_nonsynom.csv \
./data/reference/all_meta.sequence_success.csv ./data/processed/secondary/Intrahost_all.csv ./scripts/secondary/Figures/Figure2.R
	Rscript ./scripts/secondary/Figures/Figure2.R

####################### Figure 3 ######################

./results/Figures/Figure3A.pdf ./results/Figures/Figure3B.pdf ./results/Figures/Figure3C.pdf ./results/Figures/Figure3D.pdf : Figure3.intermediate
        #Empty recipe to propagate "newness" from the intermediate to final targets
.INTERMEDIATE: Figure3.intermediate

Figure3.intermediate: ./data/processed/secondary/possible.pairs.dist.csv ./data/processed/secondary/transmission_pairs.csv \
./data/processed/secondary/trans_freq.csv ./data/reference/all_meta.sequence_success.csv \
./data/processed/secondary/transmission_pairs_freq.poly.donor.csv ./data/reference/accuracy_stringent.csv \
./scripts/secondary/Figures/Figure3.R
	Rscript ./scripts/secondary/Figures/Figure3.R

####################### Figure 4 ######################

./results/Figures/Figure4A.pdf ./results/Figures/Figure4B.pdf : Figure4.intermediate
        #Empty recipe to propagate "newness" from the intermediate to final targets
.INTERMEDIATE: Figure4.intermediate
Figure4.intermediate: ./data/processed/secondary/no_cut_trans_freq.csv \
./data/processed/secondary/no_cut_transmission_pairs_freq.poly.donor.csv ./scripts/secondary/Figures/Figure4.R
	Rscript ./scripts/secondary/Figures/Figure4.R


# --------------------- Secondary Processing ---------------------


./data/processed/secondary/community_pairs.freq.csv ./data/processed/secondary/community_pairs_freq.poly.donor.csv \
./data/processed/secondary/no_cut_transmission_pairs_freq.poly.donor.csv ./data/processed/secondary/no_cut_trans_freq.csv \
./data/processed/secondary/transmission_pairs_freq.poly.donor.csv ./data/processed/secondary/trans_freq.csv : transmission.setup.intermediate
	#Empty recipe to propagate "newness" from the intermediate to final targets
.INTERMEDIATE: transmission.setup.intermediate
transmission.setup.intermediate: ./data/processed/secondary/qual.snv.csv ./data/processed/secondary/no_freq_cut.qual.snv.csv \
./data/processed/secondary/transmission_pairs.csv
	Rscript ./scripts/secondary_analsysis/processing/settingup_transmission.R


./data/processed/secondary/Intrahost_initially_present.csv ./data/processed/secondary/Intrahost_pairs.csv \
./data/processed/secondary/Intrahost_all.csv: intrahost.setup.intermediate
	# Empty recipe to propagate "newness" from the intermediate to final targets
.INTERMEDIATE: intrahost.setup.intermediate
intrahost.setup.intermediate: ./data/processed/secondary/qual.snv.csv ./data/reference/all_meta.sequence_success.csv ./scripts/R/Intrahost_setup.R
	Rscript ./scripts/secondary_analsysis/processing/Intrahost_setup.R

./data/processed/secondary/possible.pairs.dist.csv ./data/processed/secondary/transmission_pairs.csv: distance.intermediate
	# Empty recipe to propagate "newness" from the intermediate to final targets

.INTERMEDIATE: distance.intermediate
distance.intermediate: ./data/processed/secondary/qual.snv.csv ./scripts/R/L1_norm.R
	Rscript ./scripts/secondary_analsysis/processing/L1_norm.R


./data/processed/secondary/qual.snv.csv ./data/processed/secondary/no_freq_cut.qual.snv.csv : qual.intermediate
	# Empty recipe to propagate "newness" from the intermediate to final targets
	
.INTERMEDIATE: qual.intermediate
qual.intermediate: ./data/processed/secondary/average_coverages.csv ./scripts/R/processing_snv.R
	Rscript ./scripts/secondary_analsysis/processing/processing_snv.R

./data/processed/secondary/average_coverages.csv: ./scripts/R/processing_coverage.R
	Rscript ./scripts/secondary_analsysis/processing/processing_coverage.R

