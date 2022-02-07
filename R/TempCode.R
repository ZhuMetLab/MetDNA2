

# creat demo data for substructure annotation 200619 ----------------------
#
# load('I:/00_projects/03_MetDNA2/00_data/20200520_zhulib_annotation_with_msfinder_lib/merged_info_zhulib/cpd_info_zhulib.RData')
# load('I:/00_projects/03_MetDNA2/00_data/20200520_zhulib_annotation_with_msfinder_lib/merged_info_zhulib/zhumetlib.rda')
# load('I:/00_projects/03_MetDNA2/00_data/20200520_zhulib_annotation_with_msfinder_lib/lib_class_substructure_200527.RData')
#
# ms2_spec <- ImmsTools::GetMsMs(database = zhumetlib,
#                                polarity = 'neg',
#                                id = 'L0067',
#                                instrument = 'Agilent',
#                                ce = '20')
# mz_dev <- 0.0005
# mz_precursor <- 505.9884
# chem_class <- 'Purine ribonucleoside triphosphates'
#
# demo_data_substructure <- list(cpd_name = 'ATP_20v_neg',
#                                mz_precursor = 505.9884,
#                                chem_class = 'Purine ribonucleoside triphosphates',
#                                ms2_spec = ms2_spec)
#
# dir.create("./extdata", showWarnings = FALSE, recursive = TRUE)
# save(demo_data_substructure,
#      file = './extdata/demo_data_substructure.RData',
#      version = 2)


# lib_class_substructure --------------------------------------------------
# # 200619
# load('I:/00_projects/03_MetDNA2/00_data/20200520_zhulib_annotation_with_msfinder_lib/lib_class_substructure_200527.RData')
# devtools::use_data(lib_class_substructure)

# lib_adduct_nl ----------------------------------------------------------------
# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/lib_adduct_nl_200714.RData')
#
# devtools::use_data(lib_adduct_nl)

# # lib_formula ------------------------------------------------------------------
# lib_formula <- readr::read_tsv('I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/MsfinderFormulaDB-VS10.efd')
# lib_formula <- lib_formula %>%
#   na_if(., 'N/A') %>%
#   type_convert()
# colnames(lib_formula) <- c('exact_mass', 'formula', 'id_pubchem',
#                            'records', 'hmdb_included', 'knapsack_included',
#                            'chebi_included', 'drugbank_included', 'smpdb_included',
#                            'ymdb_included', 't3db_included',
#                            'foodb_included', 'nanpdb_included', 'stoff_included',
#                            'bmdb_included', 'lmsb_included', 'urine_included',
#                            'saliva_included', 'feces_included', 'ecmdb_included',
#                            'csf_included', 'serum_included', 'pubchem_included',
#                            'plantcyc_included', 'unpd_included', 'mine_included')
#
# devtools::use_data(lib_formula, overwrite = TRUE)

# cpd_zhulib -------------------------------------------------------------------
# load('I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/demo_formula/cpd_zhulib_200811.RData')
# devtools::use_data(cpd_zhulib)



################################################################################
# 2020 update MetDNA libs ------------------------------------------------------
#
# library(devtools)
#
# load('H:/00_projects/03_MetDNA2/06_files_for_package/library/20200922_zhumetlib/zhuMetlib_200923.RData')
# use_data(zhuMetlib, version = 2, overwrite = TRUE)
#
# load('H:/00_projects/03_MetDNA2/06_files_for_package/adduct_table/200923_adduct_for_annotation/lib_adduct_annot_200923.RData')
# use_data(lib_adduct_annot, version = 2, overwrite = TRUE)


################################################################################
# 20200928 update MetDNA ZhuMetLib ---------------------------------------------
# library(devtools)
#
# load('H:/00_projects/03_MetDNA2/06_files_for_package/library/20200928/zhuMetlib_200928.RData')
# use_data(zhuMetlib, version = 2, overwrite = TRUE)


################################################################################
# 20200929 use inHouse_cpd_md and kegg_cpd_md ----------------------------------
# load('D:/01_r_package/MetDNA_V1.21/data/inHouse.compound.md.rda')
# load('D:/01_r_package/MetDNA_V1.21/data/kegg.compound.md.rda')
# load('D:/01_r_package/MetDNA_V1.21/data/kegg.compound.rda')
#
# md_inHouse_cpd <- inHouse.compound.md
# use_data(md_inHouse_cpd, overwrite = TRUE, version = 2)
#
# md_kegg_cpd <- kegg.compound.md
# use_data(md_kegg_cpd, overwrite = TRUE, version = 2)
#
# cpd_kegg <- kegg.compound
# use_data(cpd_kegg, overwrite = TRUE, version = 2)

# # replace md_inHouse_cpd with new version
# load('H:/00_projects/03_MetDNA2/06_files_for_package/md/200920/md_zhumetlib_200920.RData')
# md_inHouse_cpd <- md_zhumetlib
# use_data(md_inHouse_cpd, overwrite = TRUE, version = 2)



################################################################################
# 20200929 use kegg metabolite and
# load('D:/01_r_package/MetDNA_V1.21/data/kegg.compound.rda')



################################################################################
# 20201007 update cpd_zhuMetLib, cpd_emrn, md_zhuMetLib, md_emrn ---------------
# library(devtools)
# load('H:/00_projects/03_MetDNA2/06_files_for_package/adduct_table/201007/lib_adduct_nl_201007.RData')
# use_data(lib_adduct_nl, overwrite = TRUE, version = 2)
#
# load('D:/01_r_package/MetIMMS102/data/lib_rt_exp.rda')
# temp_inchikey1 <- lib_rt_exp$inchikey %>% ImmsTools::SplitInchiKey()
# lib_rt_exp <- lib_rt_exp %>% mutate(inchikey1 = temp_inchikey1$inchikey1)
# lib_rt <- list('zhumetlib_exp' = lib_rt_exp)
# use_data(lib_rt, overwrite = TRUE, version = 2)
#
# load('D:/01_r_package/MetIMMS102/data/ccs_lib_zhumetlibPredCCS.rda')
# load('D:/01_r_package/MetIMMS102/data/lib_ccs_exp.rda')
# lib_ccs <- list('zhumetlib_exp' = lib_ccs_exp,
#                 'zhumetlib_pred' = ccs_lib_zhumetlibPredCCS)
# use_data(lib_ccs, overwrite = TRUE, version = 2)


################################################################################
# 20201009 update zhuMetLib to remove inchikey1 redunancy ---------------
# library(devtools)
# load('H:/00_projects/03_MetDNA2/06_files_for_package/library/20201009/zhuMetlib_201009.RData')
# use_data(zhuMetlib, overwrite = TRUE, version = 2)


################################################################################
# 20201010 update zhuMetLib to remove inchikey1 redunancy ---------------
# library(devtools)
# load('H:/00_projects/03_MetDNA2/06_files_for_package/library/20201010/zhuMetlib_201010.RData')
# use_data(zhuMetlib, overwrite = TRUE, version = 2)


################################################################################
# 20201012 update kegg cpds ---------------------------------------------
## load('D:/01_r_package/MetDNA_V1.21/data/inHouse.compound.md.rda')
# library(devtools)
# load('H:/00_projects/03_MetDNA2/06_files_for_package/md/201012/md_zhumetlib_201012.RData')
# use_data(md_zhumetlib, overwrite = TRUE, version = 2)
#
# load('D:/01_r_package/MetDNA_V1.21/data/kegg.compound.md.rda')
# load('H:/00_projects/03_MetDNA2/06_files_for_package/md/201007/md_emrn_201007.RData')
#
# md_mrn_emrn <- list('version1' = kegg.compound.md,
#                     'version2' = md_emrn)
# use_data(md_mrn_emrn, overwrite = TRUE, version = 2)
#
# load('D:/01_r_package/MetDNA_V1.21/data/kegg.compound.rda')
# load('H:/00_projects/03_MetDNA2/06_files_for_package/library/20201007/cpd_emrn_201007.RData')
# use_data(cpd_emrn, overwrite = TRUE, version = 2)
#
# load('D:/01_r_package/MetDNA_V1.21/data/kegg.rpair2.rda')
# load('H:/00_projects/03_MetDNA2/06_files_for_package/mrn/200929/emrn_rpair_200929.RData')
# rpair_emrn <- emrn_rpair
# rpair_mrn <- kegg.rpair2
# reaction_pair_network <- list('version1' = rpair_mrn,
#                               'version2' = rpair_emrn)
#
# use_data(reaction_pair_network, overwrite = TRUE, version = 2)


################################################################################
# 20201014 update cpd_emrn -----------------------------------------------------
# load('H:/00_projects/03_MetDNA2/06_files_for_package/library/20201014/cpd_emrn_201014.RData')
# library(devtools)
# use_data(cpd_emrn, overwrite = TRUE, version = 2)


################################################################################
# 20201103 update lib_rt, rt_reference -----------------------------------------
# load('/home/zhouzw/Data_processing/20201102_metdna2_rt_lib_update/RecalibrateRT/rt_ref_201103.RData')
# library(devtools)
# use_data(rt_ref, overwrite = TRUE, version = 2)
#
# load('/home/zhouzw/Data_processing/20201102_metdna2_rt_lib_update/RecalibrateRT/lib_rt_201103.RData')
# use_data(lib_rt, overwrite = TRUE, version = 2)
#
# load('/home/zhouzw/Data_processing/20201102_update_libraries/cpd_emrn_201102.RData')
# use_data(cpd_emrn, overwrite = TRUE, version = 2)
#
# load('/home/zhouzw/Data_processing/20201102_update_libraries/lib_adduct_nl_201102.RData')
# use_data(lib_adduct_nl, overwrite = TRUE, version = 2)
#
# load('/home/zhouzw/Data_processing/20201102_update_libraries/lib_formula_201103.rda')
# use_data(lib_formula, overwrite = TRUE, version = 2)


################################################################################
# # 20201106 update RT prediction model ------------------------------------------
# load('/home/zhouzw/Data_processing/20201102_metdna2_rt_lib_update/RT_prediction_model/201106/rf_reg_Amide12min_201106.RData')
# rf_reg_amide12min <- rf.reg
# load('/home/zhouzw/Data_processing/20201102_metdna2_rt_lib_update/RT_prediction_model/201106/rf_reg_Amide23min_201106.RData')
# rf_reg_amide23min <- rf.reg
#
# rt_model_pretrained <- list('Amide12min' = rf_reg_amide12min,
#                             'Amide23min' = rf_reg_amide23min)
#
# library(devtools)
# use_data(rt_model_pretrained, version = 2, overwrite = TRUE)

## only keep used MDs of emrn --------------------------------------------------
# load('/home/zhouzw/Data_processing/20201102_update_libraries/md_emrn_201106.RData')
# md_mrn_emrn$version2 <- md_emrn
# library(devtools)
# use_data(md_mrn_emrn, overwrite = TRUE, version = 2)

################################################################################
# # 20201118 update lib_adduct_nl (remove some uncommon adduct and nl in credentrial)
#
# library(devtools)
# load('/home/zhouzw/Data_processing/20201118_update_lib_adduct_nl/lib_adduct_nl_201118.RData')
# use_data(lib_adduct_nl, overwrite = TRUE, version = 2)

################################################################################
# # 20201202 update zhuMetLib and md_zhumetlib ---------------------------------
# library(devtools)
#
# load('/home/zhouzw/Data_processing/20201202_update_cpd_zhumetlib_and_md_zhumetlib/zhuMetlib_201201.RData')
# load('/home/zhouzw/Data_processing/20201202_update_cpd_zhumetlib_and_md_zhumetlib/md_zhumetlib_201202.RData')
#
# use_data(zhuMetlib, overwrite = TRUE, version = 2)
# use_data(md_zhumetlib, overwrite = TRUE, version = 2)


################################################################################
# # 20201212 update ZhuMetLib --------------------------------------------------
#
# library(devtools)
# load('/home/zhouzw/Data_processing/20201212_update_cpd_zhumetlib/zhuMetlib_201212.RData')
# use_data(zhuMetlib, overwrite = TRUE, version = 2)


################################################################################
# # 20210105 add fiehnHilicLib -------------------------------------------------
#
# library(devtools)
# load('/home/zhouzw/Data_processing/20210105_add_fiehn_hilic_lib/fiehnHilicLib_210105.RData')
# use_data(fiehnHilicLib, overwrite = TRUE, version = 2)
# load('/home/zhouzw/Data_processing/20210105_add_fiehn_hilic_lib/md_fiehnHilicLib_210105.RData')
# md_fiehnHilicLib <- md_fiehnHilicLib %>% tibble::column_to_rownames(var = 'name')
# use_data(md_fiehnHilicLib, overwrite = TRUE, version = 2)


################################################################################
# # 20210223 update extended MRN -----------------------------------------------
#
# update MRN and EMRN ----------------------------------------------------------
# load('/home/zhouzw/Data_processing/20210224_update_emrn_md_formula_lib/kegg.rpair2.rda')
# load('/home/zhouzw/Data_processing/20210224_update_emrn_md_formula_lib/emrn_rpair_210223.RData')
# rpair_emrn <- emrn_rpair
# rpair_mrn <- kegg.rpair2
# reaction_pair_network <- list('version1' = rpair_mrn,
#                               'version2' = rpair_emrn)
# usethis::use_data(reaction_pair_network, overwrite = TRUE, version = 2)
#
#
# update MD --------------------------------------------------------------------
# load('/home/zhouzw/Data_processing/20210224_update_emrn_md_formula_lib/kegg.compound.md.rda')
# load('/home/zhouzw/Data_processing/20210224_update_emrn_md_formula_lib/md_emrn_final_210224.RData')
# md_mrn_emrn <- list('version1' = kegg.compound.md,
#                     'version2' = md_emrn)
# usethis::use_data(md_mrn_emrn, overwrite = TRUE, version = 2)
#
#
# update compound library ------------------------------------------------------
# load('/home/zhouzw/Data_processing/20210224_update_emrn_md_formula_lib/cpd_emrn_210223.RData')
# usethis::use_data(cpd_emrn, overwrite = TRUE, version = 2)
#
#
# update formula library -------------------------------------------------------
# load('/home/zhouzw/Data_processing/20210224_update_emrn_md_formula_lib/lib_formula_210223.RData')
# usethis::use_data(lib_formula, overwrite = TRUE, version = 2)
# # update the corrected lib_formula
# load('/home/zhouzw/Data_processing/20210224_update_emrn_md_formula_lib/lib_formula_210308.RData')
# usethis::use_data(lib_formula, overwrite = TRUE, version = 2)

################################################################################
# # 202010309 add lib_pathway --------------------------------------------------
# load('/home/zhouzw/Data_processing/20210309_update_lib_pathway/list_pathway_210303.RData')
# usethis::use_data(lib_pathway, overwrite = TRUE, version = 2)

################################################################################
# # 20210310 update zhuMetlib, ZhuMetlib_obitrap, lib_rt, md_zhumetlib ---------
#
# load('/home/zhouzw/Data_processing/20210310_update_zhumetlib_zhuObitrapLib_rtLib_md/lib_rt_210310.RData')
# usethis::use_data(lib_rt, overwrite = TRUE, version = 2)
#
# load('/home/zhouzw/Data_processing/20210310_update_zhumetlib_zhuObitrapLib_rtLib_md/md_zhumetlib_210310.RData')
# usethis::use_data(md_zhumetlib, overwrite = TRUE, version = 2)
#
# load('/home/zhouzw/Data_processing/20210310_update_zhumetlib_zhuObitrapLib_rtLib_md/zhuMetlib_210310.RData')
# usethis::use_data(zhuMetlib, overwrite = TRUE, version = 2)
#
# load('/home/zhouzw/Data_processing/20210310_update_zhumetlib_zhuObitrapLib_rtLib_md/zhuMetlib_obitrap_210310.RData')
# usethis::use_data(zhuMetlib_obitrap, overwrite = TRUE, version = 2)

################################################################################
# # 20210312 update lib_adduct_nl ----------------------------------------------
#
# load('/home/zhouzw/Data_processing/20210312_update_lib_adduct_nl/lib_adduct_nl_210312.RData')
# usethis::use_data(lib_adduct_nl, overwrite = TRUE, version = 2)

################################################################################
# # 20210315 update lib_ccs ----------------------------------------------------
#
# load('/home/zhouzw/Data_processing/20210315_update_lib_ccs/lib_ccs_210315.RData')
# lib_ccs$emrnlib_pred$adduct[lib_ccs$emrnlib_pred$adduct == '[M+H-H2O]+'] <- '[M-H2O+H]+'
# usethis::use_data(lib_ccs, overwrite = TRUE, version = 2)

################################################################################
# # 20210326 add mrn_v0.3 for evaluation ---------------------------------------
#
# load('/home/zhouzw/Data_processing/20210326_metdna2_mrn_0.3/MetDNA2/data/cpd_emrn.rda')
# cpd_emrn_v03 <- cpd_emrn
# usethis::use_data(cpd_emrn_v03, overwrite = TRUE, version = 2)
#
# load('/home/zhouzw/Data_processing/20210326_metdna2_mrn_0.3/MetDNA2/data/md_mrn_emrn.rda')
# md_mrn_emrn_v03 <- md_mrn_emrn
# usethis::use_data(md_mrn_emrn_v03, overwrite = TRUE, version = 2)
#
# load('/home/zhouzw/Data_processing/20210326_metdna2_mrn_0.3/MetDNA2/data/reaction_pair_network.rda')
# reaction_pair_network_v03 <- reaction_pair_network
# usethis::use_data(reaction_pair_network_v03, overwrite = TRUE, version = 2)

################################################################################
# # 20210329 update mrn_v0.5 ---------------------------------------------------
#
# load('/home/zhouzw/Data_processing/20210224_update_emrn_md_formula_lib/kegg.rpair2.rda')
# load('/home/zhouzw/Data_processing/20210329_update_metdna2_mrn_v0.5/emrn_rpair_210329.RData')
# rpair_emrn <- emrn_rpair
# rpair_mrn <- kegg.rpair2
# reaction_pair_network <- list('version1' = rpair_mrn,
#                               'version2' = rpair_emrn)
# usethis::use_data(reaction_pair_network, overwrite = TRUE, version = 2)

################################################################################
# # 20210330 update ZhuMetlib_obitrap v3 ---------------------------------------
#
# load('/home/zhouzw/Data_processing/20210330_update_zhumetlib_orbitrap_v3/zhuMetlib_obitrap_210330.RData')
# usethis::use_data(zhuMetlib_obitrap, overwrite = TRUE, version = 2)

################################################################################
# # 20210403 update zhuMetlib and zhuMetlib_obitrap (compound information) -----
#
# load('/home/zhouzw/Data_processing/20210403_update_zhuMetlib_zhuMetlib_orbitrap/zhuMetlib_20210403.RData')
# usethis::use_data(zhuMetlib, overwrite = TRUE, version = 2)
#
# load('/home/zhouzw/Data_processing/20210403_update_zhuMetlib_zhuMetlib_orbitrap/zhuMetlib_obitrap_20210403.RData')
# usethis::use_data(zhuMetlib_obitrap, overwrite = TRUE, version = 2)

################################################################################
# # 20210413 change zhuMetLib_obitrap to zhuMetLib_orbitrap --------------------
#
# load('/home/zhouzw/Data_processing/20210403_update_zhuMetlib_zhuMetlib_orbitrap/zhuMetlib_obitrap_20210403.RData')
# zhuMetlib_orbitrap <- zhuMetlib_obitrap
# usethis::use_data(zhuMetlib_orbitrap, overwrite = TRUE, version = 2)

################################################################################
# # 20210515 update md_mrn_emrn$version1 ---------------------------------------
# load('H:/00_projects/03_MetDNA2/00_data/20210514_modify_MetDNA1/02_new_RData/md_kegg_metdna1_210515.RData')
# load('/home/zhouzw/Data_processing/20210516_update_md_kegg_metdna1/md_kegg_metdna1_210515.RData')
# data("md_mrn_emrn")
# md_mrn_emrn$version1 <- md_kegg_metdna1
# usethis::use_data(md_mrn_emrn, overwrite = TRUE)


################################################################################
# # 20210615 update lib_adduct_nl ----------------------------------------------
# adduct_nl_pos <- readr::read_csv('/home/zhouzw/Data_processing/20210615_update_adduct_nl_lib/lib_adduct_nl_pos_210615_manual.csv')
# adduct_nl_neg <- readr::read_csv('/home/zhouzw/Data_processing/20210615_update_adduct_nl_lib/lib_adduct_nl_neg_210615_manual.csv')
# lib_adduct_nl$positive <- adduct_nl_pos
# lib_adduct_nl$negative <- adduct_nl_neg
#
# usethis::use_data(lib_adduct_nl, overwrite = TRUE)
#
# load('/home/zhouzw/Data_processing/20210615_update_adduct_nl_lib/zhuMetlib_210615.RData')
# zhuMetlib$meta$compound
#
# load('/home/zhouzw/Data_processing/20210615_update_adduct_nl_lib/zhuMetlib_obitrap_210615.RData')
# zhuMetlib_orbitrap <- zhuMetlib_obitrap
#
# usethis::use_data(zhuMetlib, overwrite = TRUE)
# usethis::use_data(zhuMetlib_orbitrap, overwrite = TRUE)

################################################################################
# # 20210616 update lib_adduct_nl, zhuMetlib, zhuMetlib_orbitrap ---------------
#
# adduct_nl_pos <- readxl::read_xlsx('/home/zhouzw/Data_processing/20210615_update_adduct_nl_lib/lib_adduct_nl_pos_210616_manual.xlsx')
# adduct_nl_neg <- readxl::read_xlsx('/home/zhouzw/Data_processing/20210615_update_adduct_nl_lib/lib_adduct_nl_neg_210616_manual.xlsx')
# lib_adduct_nl$positive <- adduct_nl_pos
# lib_adduct_nl$negative <- adduct_nl_neg
# usethis::use_data(lib_adduct_nl, overwrite = TRUE)
#
# load('/home/zhouzw/Data_processing/20210615_update_adduct_nl_lib/zhuMetlib_20210616.RData')
# load('/home/zhouzw/Data_processing/20210615_update_adduct_nl_lib/zhuMetlib_orbitrap_20210616.RData')
# usethis::use_data(zhuMetlib, overwrite = TRUE)
# usethis::use_data(zhuMetlib_orbitrap, overwrite = TRUE)

################################################################################
# # 20210621 update zhuMetlib and zhuMetlib_orbitrap ---------------------------
# load('/home/zhouzw/Data_processing/20210615_update_adduct_nl_lib/zhuMetlib_20210621.RData')
# load('/home/zhouzw/Data_processing/20210615_update_adduct_nl_lib/zhuMetlib_orbitrap_20210621.RData')
# usethis::use_data(zhuMetlib, overwrite = TRUE)
# usethis::use_data(zhuMetlib_orbitrap, overwrite = TRUE)

################################################################################
# # 20210629 update cpd_emrn (add name for some in-silico metabolites in KEGG)
# load('/home/zhouzw/Data_processing/20210629_update_emrn/cpd_emrn_210629.RData')
# usethis::use_data(cpd_emrn, overwrite = TRUE)

################################################################################
# # 20210720 add 200STD compound related files ---------------------------------
# load('/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/files_for_200std/cpd_200stdExd_210720.RData')
# usethis::use_data(cpd_200stdExd, overwrite = TRUE)
# load('/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/files_for_200std/md_200std_210720.RData')
# usethis::use_data(md_200std, overwrite = TRUE)
# load('/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/files_for_200std/reaction_pair_network_200std_210720.RData')
# usethis::use_data(reaction_pair_network_200std, overwrite = TRUE)

################################################################################
# # 20210823 add 46STD compound related files ----------------------------------
# load('/home/zhouzw/Data_processing/20210823_modify_MetDNA2_for_46std_validation/cpd_46stdExd_210823.RData')
# usethis::use_data(cpd_46stdExd, overwrite = TRUE)
# load('/home/zhouzw/Data_processing/20210823_modify_MetDNA2_for_46std_validation/md_46std_210823.RData')
# usethis::use_data(md_46std, overwrite = TRUE)
# load('/home/zhouzw/Data_processing/20210823_modify_MetDNA2_for_46std_validation/reaction_pair_network_46std_210823.RData')
# usethis::use_data(reaction_pair_network_46std, overwrite = TRUE)

################################################################################
# # 20210928 add metlin RT library ---------------------------------------------
# load('/home/zhouzw/Data_processing/20210928_metlin_rt/lib_rt_210928.RData')
# usethis::use_data(lib_rt, overwrite = TRUE)
#
# load('/home/zhouzw/Data_processing/20210928_metlin_rt/rt_ref_210928.RData')
# usethis::use_data(rt_ref, overwrite = TRUE)

################################################################################
# # 20211028 add zhulab RP library ---------------------------------------------
#
# load('/home/zhouzw/Data_processing/20211028_zhulab_rp_library/lib_rt_211028.RData')
# load('/home/zhouzw/Data_processing/20211028_zhulab_rp_library/zhuRPlib_211028.RData')
# load('/home/zhouzw/Data_processing/20211028_zhulab_rp_library/rt_ref_211028.RData')
# load('/home/zhouzw/Data_processing/20211028_zhulab_rp_library/md_zhuRPLib_211028.RData')

# usethis::use_data(lib_rt, overwrite = TRUE)
# usethis::use_data(zhuRPlib, overwrite = TRUE)
# usethis::use_data(rt_ref, overwrite = TRUE)
# usethis::use_data(md_zhuRPLib, overwrite = TRUE)

################################################################################
# # 20211101 update zhuRPLib ---------------------------------------------------
#
# load('/home/zhouzw/Data_processing/20211028_zhulab_rp_library/zhuRPlib_211028.RData')
# zhuRPlib$meta$pos$`SCE20-30-40` %>% as_tibble()
# zhuRPlib$meta$pos$`SNCE20-30-40`$mz <- as.numeric(zhuRPlib$meta$pos$`SNCE20-30-40`$mz)
# zhuRPlib$meta$neg$`SNCE20-30-40`$mz <- as.numeric(zhuRPlib$meta$neg$`SNCE20-30-40`$mz)
# zhuRPlib$meta$compound$monoisotopic_mass <- as.numeric(zhuRPlib$meta$compound$monoisotopic_mass)
#
# usethis::use_data(zhuRPlib, overwrite = TRUE)
# zhuRPlib$meta$compound$name <- stringr::str_replace_all(zhuRPlib$meta$compound$name, pattern = ';', replacement = '^')
# usethis::use_data(zhuRPlib, overwrite = TRUE)

################################################################################
# # 20211129 update zhuRPLib ---------------------------------------------------
#
# load('/home/zhouzw/Data_processing/20211129_zhulab_rp_library_update/lib_rt_211129.RData')
# load('/home/zhouzw/Data_processing/20211129_zhulab_rp_library_update/zhuRPlib_211129.RData')
#
# usethis::use_data(lib_rt, overwrite = TRUE)
# usethis::use_data(zhuRPlib, overwrite = TRUE)
