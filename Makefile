.PHONY: figures secondary_analysis
# --------------------- Figures ---------------------

figures : ./results/Figures/Figure1A.pdf ./results/Figures/Figure1B.pdf \
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
./data/reference/segs.csv 
	Rscript ./scripts/secondary/Figures/Figure1.R

####################### Figure 2 ######################

./results/Figures/Figure2A.pdf ./results/Figures/Figure2B.pdf \
./results/Figures/Figure2C.pdf ./results/Figures/Figure2D.pdf ./results/Figures/Figure2E.pdf : Figure2.intermediate
        #Empty recipe to propagate "newness" from the intermediate to final targets
.INTERMEDIATE: Figure2.intermediate
Figure2.intermediate: ./data/processed/secondary/qual.snv.csv ./data/processed/secondary/antigenic_isnv.csv \
./data/processed/secondary/global_freq_antigenic.tsv ./data/processed/secondary/minor_nonsynom.csv \
./data/reference/all_meta.sequence_success.csv ./data/processed/secondary/Intrahost_all.csv 
	Rscript ./scripts/secondary/Figures/Figure2.R

####################### Figure 3 ######################

./results/Figures/Figure3A.pdf ./results/Figures/Figure3B.pdf ./results/Figures/Figure3C.pdf ./results/Figures/Figure3D.pdf : Figure3.intermediate
        #Empty recipe to propagate "newness" from the intermediate to final targets
.INTERMEDIATE: Figure3.intermediate

Figure3.intermediate: ./data/processed/secondary/possible.pairs.dist.csv ./data/processed/secondary/transmission_pairs.csv \
./data/processed/secondary/trans_freq.csv ./data/reference/all_meta.sequence_success.csv \
./data/processed/secondary/transmission_pairs_freq.poly.donor.csv ./data/reference/accuracy_stringent.csv 
	Rscript ./scripts/secondary/Figures/Figure3.R

####################### Figure 4 ######################

./results/Figures/Figure4A.pdf ./results/Figures/Figure4B.pdf : Figure4.intermediate
        #Empty recipe to propagate "newness" from the intermediate to final targets
.INTERMEDIATE: Figure4.intermediate
Figure4.intermediate: ./data/processed/secondary/no_cut_trans_freq.csv \
./data/processed/secondary/no_cut_transmission_pairs_freq.poly.donor.csv 
	Rscript ./scripts/secondary/Figures/Figure4.R


# --------------------- Secondary Processing ---------------------

secondary_analysis: ./data/processed/secondary/community_pairs.freq.csv ./data/processed/secondary/community_pairs_freq.poly.donor.csv \
./data/processed/secondary/no_cut_transmission_pairs_freq.poly.donor.csv ./data/processed/secondary/no_cut_trans_freq.csv \
./data/processed/secondary/transmission_pairs_freq.poly.donor.csv ./data/processed/secondary/trans_freq.csv \
./data/processed/secondary/Intrahost_initially_present.csv ./data/processed/secondary/Intrahost_pairs.csv \
./data/processed/secondary/Intrahost_all.csv \
./data/processed/secondary/possible.pairs.dist.csv ./data/processed/secondary/transmission_pairs.csv \
./data/processed/secondary/qual.snv.csv ./data/processed/secondary/no_freq_cut.qual.snv.csv \
./data/processed/secondary/average_coverages.csv
	echo "Done with secondary analysis"

./data/processed/secondary/community_pairs.freq.csv ./data/processed/secondary/community_pairs_freq.poly.donor.csv \
./data/processed/secondary/no_cut_transmission_pairs_freq.poly.donor.csv ./data/processed/secondary/no_cut_trans_freq.csv \
./data/processed/secondary/transmission_pairs_freq.poly.donor.csv ./data/processed/secondary/trans_freq.csv : transmission.setup.intermediate
	#Empty recipe to propagate "newness" from the intermediate to final targets
.INTERMEDIATE: transmission.setup.intermediate
transmission.setup.intermediate: ./data/processed/secondary/qual.snv.csv ./data/processed/secondary/no_freq_cut.qual.snv.csv \
./data/processed/secondary/transmission_pairs.csv
	Rscript ./scripts/secondary_analysis/processing/settingup_transmission.R


./data/processed/secondary/Intrahost_initially_present.csv ./data/processed/secondary/Intrahost_pairs.csv \
./data/processed/secondary/Intrahost_all.csv: intrahost.setup.intermediate
	# Empty recipe to propagate "newness" from the intermediate to final targets
.INTERMEDIATE: intrahost.setup.intermediate
intrahost.setup.intermediate: ./data/processed/secondary/qual.snv.csv ./data/reference/all_meta.sequence_success.csv 
	Rscript ./scripts/secondary_analysis/processing/Intrahost_setup.R

./data/processed/secondary/possible.pairs.dist.csv ./data/processed/secondary/transmission_pairs.csv: distance.intermediate
	# Empty recipe to propagate "newness" from the intermediate to final targets

.INTERMEDIATE: distance.intermediate
distance.intermediate: ./data/processed/secondary/qual.snv.csv 
	Rscript ./scripts/secondary_analysis/processing/L1_norm.R


./data/processed/secondary/qual.snv.csv ./data/processed/secondary/no_freq_cut.qual.snv.csv ./data/reference/all_meta.sequence_success.csv: qual.intermediate
	# Empty recipe to propagate "newness" from the intermediate to final targets
	
.INTERMEDIATE: qual.intermediate
qual.intermediate: ./data/processed/secondary/average_coverages.csv \
./data/processed/HK_1/all.variants.csv ./data/processed/HK_2/all.variants.csv ./data/processed/HK_6/all.variants.csv \
./data/processed/HK_7/all.variants.csv ./data/processed/HK_8/all.variants.csv ./data/processed/cali09/all.variants.csv \
./data/processed/cali09_2/all.variants.csv ./data/processed/victoria/all.variants.csv ./data/processed/victoria_2/all.variants.csv \
./data/processed/perth/all.variants.csv ./data/processed/perth_2/all.variants.csv ./data/reference/all_meta.csv
	Rscript ./scripts/secondary_analysis/processing/processing_snv.R


# Empty targes for the all.varaint files - These are made by the variant_calling_pipeline and not by the make file
./data/processed/HK_1/all.variants.csv: ;
./data/processed/HK_2/all.variants.csv: ;
./data/processed/HK_6/all.variants.csv: ;
./data/processed/HK_7/all.variants.csv: ;
./data/processed/HK_8/all.variants.csv: ;
./data/processed/cali09/all.variants.csv: ;
./data/processed/cali09_2/all.variants.csv: ;
./data/processed/victoria/all.variants.csv: ;
./data/processed/victoria_2/all.variants.csv: ;
./data/processed/perth/all.variants.csv: ;
./data/processed/perth_2/all.variants.csv: ;
./data/reference/all_meta.csv: ;



./data/processed/secondary/average_coverages.csv: ./data/processed/HK_1/all.coverage.csv ./data/processed/HK_2/all.coverage.csv \
./data/processed/HK_6/all.coverage.csv ./data/processed/HK_7/all.coverage.csv ./data/processed/HK_8/all.coverage.csv \
./data/processed/cali09/all.coverage.csv ./data/processed/cali09_2/all.coverage.csv ./data/processed/victoria/all.coverage.csv \
./data/processed/victoria_2/all.coverage.csv ./data/processed/perth/all.coverage.csv ./data/processed/perth_2/all.coverage.csv
	Rscript ./scripts/secondary_analysis/processing/processing_coverage.R

