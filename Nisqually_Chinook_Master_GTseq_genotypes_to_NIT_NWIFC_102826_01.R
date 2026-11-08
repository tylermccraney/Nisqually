library(sequoia)
library(ASRgenomics)
library(haplo.stats)
library(readxl)
library(tidyverse)

setwd("~/Google Drive/My Drive/Nisqually/Nisqually/")

name_repair <-
  function(names){
    names <- 
      tolower(names) |>
      gsub(pattern="[^A-Za-z0-9_%]+", replacement="_") |> 
      gsub(pattern="%", replacement="percent", fixed=TRUE)
    is_dup <- 
      names %in% names[duplicated(names)]
    if(any(is_dup)){
      suffixes <- 
        ave(names, names, FUN=function(x) paste0("_", seq_along(x)))
      names[is_dup] <- 
        paste0(names[is_dup], suffixes[is_dup])}
    return(names)}

excel_sheets(path="../Nisqually_Chinook_Master_GTseq_genotypes_to_NIT_NWIFC_102826.xlsx")

dat <- 
  read_excel(
    path="../Nisqually_Chinook_Master_GTseq_genotypes_to_NIT_NWIFC_102826.xlsx", 
    sheet="genotypes", 
    .name_repair=name_repair,
    col_types="text")

glimpse(dat)

table(dat$code, useNA="ifany")

table(dat$green_db_collection_name, useNA="ifany")

table(dat$genotyping_year, useNA="ifany")

table(dat$life_history, useNA="ifany")

table(dat$species_id, useNA="ifany")

dat$species_id <- 
  replace_values(
    dat$species_id, 
    "no data" ~ NA, 
    "No data" ~ NA, 
    "non-Chinook" ~ "not Chinook",
    "Not Chinook" ~ "not Chinook")

table(dat$speciesid, useNA="ifany")

dat$speciesid <- 
  replace_values(
    dat$speciesid, 
    "no data" ~ NA)

table(dat$ploidy, useNA="ifany")

dat$ploidy <- 
  replace_values(
    dat$ploidy, 
    "no data" ~ NA, 
    "U" ~ NA) # assuming "U" stands for "Unknown"

table(dat$sex_id, useNA="ifany")

dat$sex_id <- 
  replace_values(
    dat$sex_id, 
    "no data" ~ NA, 
    "No data" ~ NA)

table(dat$qc, useNA="ifany")

table(dat$qc_status, useNA="ifany")



dat <-
  filter(
    dat,
    !qc_status %in% names(table(dat$qc_status))[-8],
    !qc %in% names(table(dat$qc))[-9],
    !species_id %in% "not Chinook",
    !speciesid %in% c("Oke", "Oki")) |>
  select(
    id=use_sample_name, 
    alt_id,
    green_db_collection_name,
    ifi,
    life_history, 
    ploidy,
    sex=sex_id,
    27:358) |>
  mutate(ifi=as.numeric(ifi))

glimpse(dat)



excel_sheets(
  path="../data-from-Sayre/11-27 2022-2023 Genetics Individual_wd.xlsx",
)

dat_11_27_2022_2023_Genetics_Individual_wd <- 
  read_excel(
    path="../data-from-Sayre/11-27 2022-2023 Genetics Individual_wd.xlsx",
    sheet="Sheet1", 
    .name_repair=name_repair,
    col_types="text")

glimpse(dat_11_27_2022_2023_Genetics_Individual_wd)

read_excel(
  path="../data-from-Sayre/11-27 2022-2023 Genetics Individual_wd.xlsx",
  sheet="Sheet2", 
  col_names="V1",
  col_types="text")

# unlabeled vector of genetic_label's 
# not imported




excel_sheets(
  path="../data-from-Sayre/2018-2019 Nisqually_Chinook_genotyping_report-With Bio Data (2020-06-18).xlsx",
)

read_excel(
  path="../data-from-Sayre/2018-2019 Nisqually_Chinook_genotyping_report-With Bio Data (2020-06-18).xlsx",
  sheet="report", 
  .name_repair=name_repair,
  col_types="text")

dat_2018_2019_Nisqually_Chinook_genotyping_report_With_Bio_Data <- 
  read_excel(
    path="../data-from-Sayre/2018-2019 Nisqually_Chinook_genotyping_report-With Bio Data (2020-06-18).xlsx",
    sheet="report", 
    .name_repair=name_repair,
    col_types="text")

glimpse(dat_2018_2019_Nisqually_Chinook_genotyping_report_With_Bio_Data)



excel_sheets(
  path="../data-from-Sayre/2020-2023 Genetics Individual&Batch 12-5.xlsx",
)

read_excel(
  path="../data-from-Sayre/2020-2023 Genetics Individual&Batch 12-5.xlsx",
  sheet="2021-23 Individual Gen",
  .name_repair=name_repair,
  col_types="text")

dat_1_2020_2023_Genetics_Individual_Batch_12_5 <- 
  read_excel(
    path="../data-from-Sayre/2020-2023 Genetics Individual&Batch 12-5.xlsx",
    sheet="2021-23 Individual Gen", 
    .name_repair=name_repair,
    col_types="text")

glimpse(dat_1_2020_2023_Genetics_Individual_Batch_12_5)


dat_2_2020_2023_Genetics_Individual_Batch_12_5 <- 
  read_excel(
    path="../data-from-Sayre/2020-2023 Genetics Individual&Batch 12-5.xlsx",
    sheet="2020-2021 Batch Gen",
    .name_repair=name_repair,
    col_types="text")

glimpse(dat_2_2020_2023_Genetics_Individual_Batch_12_5)









excel_sheets(
  path="../data-from-Sayre/2025 Nisqually CH Genetics.xlsx",
)

read_excel(
  path="../data-from-Sayre/2025 Nisqually CH Genetics.xlsx",
  # sheet="2021-23 Individual Gen",
  .name_repair=name_repair,
  col_types="text")

dat_2025_Nisqually_CH_Genetics <- 
  read_excel(
    path="../data-from-Sayre/2025 Nisqually CH Genetics.xlsx",
    # sheet="2021-23 Individual Gen", 
    .name_repair=name_repair,
    col_types="text")

glimpse(dat_2025_Nisqually_CH_Genetics)









excel_sheets(
  path="../data-from-Sayre/Adult to Smolt Gen. Study- B.Y. 2017 & 2019 (List of Samples To Be Analyzed).xlsx",
)

read_excel(
  path="../data-from-Sayre/Adult to Smolt Gen. Study- B.Y. 2017 & 2019 (List of Samples To Be Analyzed).xlsx",
  sheet="2017-18-19-20 Adult&Juv. Gen.",
  .name_repair=name_repair,
  col_types="text")

dat_2017_18_19_20_Adult_Juv__Gen_ <- 
  read_excel(
    path="../data-from-Sayre/Adult to Smolt Gen. Study- B.Y. 2017 & 2019 (List of Samples To Be Analyzed).xlsx",
    sheet="2017-18-19-20 Adult&Juv. Gen.",
    .name_repair=name_repair,
    col_types="text")

glimpse(dat_2017_18_19_20_Adult_Juv__Gen_)





excel_sheets(
  path="../data-from-Sayre/genetics sample numbers checking v2.xlsx",
)

read_excel(
  path="../data-from-Sayre/genetics sample numbers checking v2.xlsx",
  sheet="spawning ground",
  .name_repair=name_repair,
  col_types="text")

dat_spawning_ground <- 
  read_excel(
    path="../data-from-Sayre/genetics sample numbers checking v2.xlsx",
    sheet="spawning ground",
    .name_repair=name_repair,
    col_types="text")

glimpse(dat_spawning_ground)


dat_weir <- 
  read_excel(
    path="../data-from-Sayre/genetics sample numbers checking v2.xlsx",
    sheet="weir",
    .name_repair=name_repair,
    col_types="text")

glimpse(dat_weir)


dat_1_2020_2023_Genetics_Individual_Batch_12_5 |> names() |> sort() 
dat_11_27_2022_2023_Genetics_Individual_wd |> names() |> sort() 
dat_2017_18_19_20_Adult_Juv__Gen_ |> names() |> sort() 
dat_2018_2019_Nisqually_Chinook_genotyping_report_With_Bio_Data |> names() |> sort() 
dat_2025_Nisqually_CH_Genetics |> names() |> sort() 
dat_spawning_ground |> names() |> sort() 
dat_weir |> names() |> sort()

dat_bio <-
  bind_rows(
    `2020-2023 Genetics Individual&Batch 12-5.xlsx`=dat_1_2020_2023_Genetics_Individual_Batch_12_5, 
    `11-27 2022-2023 Genetics Individual_wd.xlsx`=dat_11_27_2022_2023_Genetics_Individual_wd, 
    `Adult to Smolt Gen. Study- B.Y. 2017 & 2019 (List of Samples To Be Analyzed).xlsx`=rename(dat_2017_18_19_20_Adult_Juv__Gen_, date=date_of_sample), 
    `2018-2019 Nisqually_Chinook_genotyping_report-With Bio Data (2020-06-18).xlsx`=dat_2018_2019_Nisqually_Chinook_genotyping_report_With_Bio_Data, 
    `2025 Nisqually CH Genetics.xlsx`=dat_2025_Nisqually_CH_Genetics, 
    `genetics sample numbers checking v2$spawning ground.xlsx`=dat_spawning_ground, 
    `genetics sample numbers checking v2$weir.xlsx`=rename(dat_weir, comment=comments, genetic_sample=genetic_sample_), 
    .id="spreadsheet")

dat_bio <- 
  mutate(
    dat_bio, 
    sample_datetime_date=as.numeric(sample_datetime) |> as.Date(origin="1899-12-30"),
    date_date=as.numeric(date) |> as.Date(origin="1899-12-30"),
    delivery_date_date=as.numeric(delivery_date) |> as.Date(origin="1899-12-30"),
    trucked_date_date=as.numeric(trucked_date) |> as.Date(origin="1899-12-30"),
    j_t_rec_date_date=as.numeric(j_t_rec_date) |> as.Date(origin="1899-12-30"))


names(dat_bio) |> sort()

distinct(dat_bio, genetic_label)
distinct(dat_bio, genetic_sample)
distinct(dat_bio, genetic_label, nit_sample_name)

distinct(dat_bio, genetic_label, nit_gen_sample_id_code, nit_sample_name, wdfw_id)

semi_join(
  x=distinct(dat, id, alt_id), 
  y=distinct(dat_bio, wdfw_id), 
  by=join_by(id==wdfw_id)) |> 
  distinct(alt_id)




# lookup1 <- 
#   read_tsv(file="Green_DB_collection_name__origin.tsv")
# 
# dat <- 
#   mutate(
#     dat, 
#     origin=recode_values(
#       Green_DB_collection_name, 
#       from=lookup1$Green_DB_collection_name, 
#       to=lookup1$origin),
#     lifestage=Life_history) |>
#   mutate(
#     origin=if_else(
#       condition=grepl(pattern="Unmarked", Alt_ID), 
#       true="natural", 
#       false=origin), 
#     lifestage=if_else(
#       condition=grepl(pattern="Jack", Alt_ID), 
#       true="jack", 
#       false=lifestage)) |>
#   mutate(
#     origin=if_else(
#       condition=grepl(pattern="AdiposeClip", Alt_ID), 
#       true="hatchery", 
#       false=origin))
# 
# 
# distinct(dat, Green_DB_collection_name, origin, lifestage) %>% 
#   print(n=nrow(.))
# 
# dat <- 
#   left_join(x=dat, y=lookup1[-2])
# 
# table(dat$sample_type, dat$lifestage, useNA="ifany")
# table(dat$Green_DB_collection_name, useNA="ifany")
# 
# gt331 <- 
#   filter(
#     dat, 
#     !QC_status %in% names(table(dat$QC_status))[-8], 
#     !species_ID %in% "not Chinook", 
#     !speciesID %in% c("Oke", "Oki"),
#     !sample_type %in% c("baseline", "return")) |> 
#   select(
#     ID=USE_sample_name, 
#     Green_DB_collection_name, 
#     origin, 
#     sample_type, 
#     lifestage, 
#     IFI, 
#     Sex=Sex_ID, 
#     28:358) |>
#   mutate(IFI=as.numeric(IFI))
# 
# glimpse(gt331)
# 
# table(gt331$sample_type, gt331$lifestage)
# 
# 
# gt331 %>% 
#   select(where(~any(grepl("-", .x, fixed=TRUE)))) %>% 
#   glimpse()
# 
# gt331$Ots_128757_61R <- 
#   gsub(pattern="-", replacement="G", x=gt331$Ots_128757_61R)
# 
# gt331$Ots_GDH_81x <- 
#   gsub(pattern="-", replacement="T", x=gt331$Ots_GDH_81x)
# 
# gt331$Ots_Ots311_101x <- 
#   gsub(pattern="-", replacement="G", x=gt331$Ots_Ots311_101x)
# 
# gt331$Ots_RFC2_558 <- 
#   gsub(pattern="-", replacement="G", x=gt331$Ots_RFC2_558)
# 
# gt331$Ots_U2305_63 <- 
#   gsub(pattern="-", replacement="C", x=gt331$Ots_U2305_63)
# 
# gt331$Ots_pigh_105 <- 
#   gsub(pattern="-", replacement="G", x=gt331$Ots_pigh_105)
# 
# gt331$Ots_u07_17_373 <- 
#   gsub(pattern="-", replacement="G", x=gt331$Ots_u07_17_373)
# 
# gt331$Ots_vatf_251 <- 
#   gsub(pattern="-", replacement="A", x=gt331$Ots_vatf_251)
# 
# gt331 <- 
#   relocate(
#     gt331, 
#     Ots28_11023212, Ots28_11025336, Ots28_11033282, Ots28_11062192, Ots37124_12267397, 
#     Ots37124_12272852, Ots28_11070757, Ots28_11071377, Ots28_11072994, Ots28_11073102, 
#     Ots37124_12277401, Ots28_11073668, Ots28_11075348, Ots28_11075712, Ots28_11077016, 
#     Ots28_11077172, Ots37124_12281207, Ots28_11077576, Ots28_11095755, Ots37124_12310649, 
#     Ots28_11143508, Ots28_11160599, Ots28_11164637, Ots28_11186543, Ots28_11201129, 
#     Ots28_11202190, Ots28_11202400, Ots28_11202863, Ots28_11205423, Ots28_11205993, 
#     Ots28_11206740, Ots28_11207428, .before=Ots28_11210919)
# 
# lh331 <- 
#   select(gt331, lifestage, sample_type, origin, ID1=ID, Sex, BirthYear=ID) |> 
#   rename(ID=ID1) |>
#   mutate(BirthYear=str_extract(BirthYear, "^\\d{2}") |> as.integer()) |>
#   mutate(BirthYear=if_else(lifestage=="juvenile", true=BirthYear+1999, false=NA))
# 
# table(lh331$sample_type, lh331$BirthYear)
# 
# 
# gt33 <- 
#   select(gt331, ID, Ots28_11023212, Ots28_11025336, Ots28_11033282, Ots28_11062192, Ots37124_12267397, Ots37124_12272852, Ots28_11070757, Ots28_11071377, Ots28_11072994, Ots28_11073102, Ots37124_12277401, Ots28_11073668, Ots28_11075348, Ots28_11075712, Ots28_11077016, Ots28_11077172, Ots37124_12281207, Ots28_11077576, Ots28_11095755, Ots37124_12310649, Ots28_11143508, Ots28_11160599, Ots28_11164637, Ots28_11186543, Ots28_11201129, Ots28_11202190, Ots28_11202400, Ots28_11202863, Ots28_11205423, Ots28_11205993, Ots28_11206740, Ots28_11207428, Ots28_11210919) |> 
#   column_to_rownames(var="ID") |>
#   filter(!if_all(everything(), is.na))
# 
# glimpse(gt33)
# 
# gm33 <- 
#   snp.recode(
#     M=as.matrix(gt33), 
#     recoding="ATGCto012")
# 
# snpstats33 <- 
#   SnpStats(gm33$Mrecode, Plot=FALSE)
# 
# gm33map <- 
#   mutate(gm33$map, snp=names(gt33))
# 
# 
# gt33 <- 
#   separate_wider_position(gt33, cols=everything(), widths=c(a1=1, a2=1), names_sep="_") |>
#   mutate(ID=rownames(gt33)) |>
#   column_to_rownames(var="ID") |>
#   filter(!if_any(everything(), is.na)) |>
#   as.matrix()
# 
# 
# ht33 <- 
#   haplo.em(geno=gt33, locus.label=tail(names(gt331), n=33))
# 
# ht33haps <- 
#   ht33$haplotype |> 
#   as_tibble()
# 
# EL01 <- 
#   mutate(
#     ht33haps, 
#     Ots28_11023212=recode_values(Ots28_11023212, "G" ~ 0, "A" ~ 1),
#     Ots28_11025336=recode_values(Ots28_11025336, "C" ~ 0, "A" ~ 1),
#     Ots28_11033282=recode_values(Ots28_11033282, "A" ~ 0, "G" ~ 1),
#     Ots28_11062192=recode_values(Ots28_11062192, "G" ~ 0, "C" ~ 1),
#     Ots37124_12267397=recode_values(Ots37124_12267397, "C" ~ 0, "T" ~ 1),
#     Ots37124_12272852=recode_values(Ots37124_12272852, "C" ~ 0, "T" ~ 1),
#     Ots28_11070757=recode_values(Ots28_11070757, "G" ~ 0, "A" ~ 1),
#     Ots28_11071377=recode_values(Ots28_11071377, "C" ~ 0, "T" ~ 1),
#     Ots28_11072994=recode_values(Ots28_11072994, "T" ~ 0, "C" ~ 1),
#     Ots28_11073102=recode_values(Ots28_11073102, "A" ~ 0, "T" ~ 1),
#     Ots37124_12277401=recode_values(Ots37124_12277401, "A" ~ 0, "T" ~ 1),
#     Ots28_11073668=recode_values(Ots28_11073668, "A" ~ 0, "T" ~ 1),
#     Ots28_11075348=recode_values(Ots28_11075348, "A" ~ 0, "G" ~ 1),
#     Ots28_11075712=recode_values(Ots28_11075712, "T" ~ 0, "C" ~ 1),
#     Ots28_11077016=recode_values(Ots28_11077016, "T" ~ 0, "C" ~ 1),
#     Ots28_11077172=recode_values(Ots28_11077172, "A" ~ 0, "G" ~ 1),
#     Ots37124_12281207=recode_values(Ots37124_12281207, "A" ~ 0, "T" ~ 1),
#     Ots28_11077576=recode_values(Ots28_11077576, "G" ~ 0, "A" ~ 1),
#     Ots28_11095755=recode_values(Ots28_11095755, "A" ~ 0, "T" ~ 1),
#     Ots37124_12310649=recode_values(Ots37124_12310649, "T" ~ 0, "A" ~ 1),
#     Ots28_11143508=recode_values(Ots28_11143508, "G" ~ 0, "A" ~ 1),
#     Ots28_11160599=recode_values(Ots28_11160599, "G" ~ 0, "T" ~ 1),
#     Ots28_11164637=recode_values(Ots28_11164637, "C" ~ 0, "A" ~ 1),
#     Ots28_11186543=recode_values(Ots28_11186543, "A" ~ 0, "T" ~ 1),
#     Ots28_11201129=recode_values(Ots28_11201129, "T" ~ 0, "G" ~ 1),
#     Ots28_11202190=recode_values(Ots28_11202190, "C" ~ 0, "T" ~ 1),
#     Ots28_11202400=recode_values(Ots28_11202400, "C" ~ 0, "T" ~ 1),
#     Ots28_11202863=recode_values(Ots28_11202863, "C" ~ 0, "A" ~ 1),
#     Ots28_11205423=recode_values(Ots28_11205423, "G" ~ 0, "A" ~ 1),
#     Ots28_11205993=recode_values(Ots28_11205993, "T" ~ 0, "C" ~ 1),
#     Ots28_11206740=recode_values(Ots28_11206740, "C" ~ 0, "T" ~ 1),
#     Ots28_11207428=recode_values(Ots28_11207428, "T" ~ 0, "G" ~ 1),
#     Ots28_11210919=recode_values(Ots28_11210919, "C" ~ 0, "T" ~ 1))
# 
# 
# EL01 <- 
#   as.matrix(EL01)
# 
# EL01s <- 
#   scale(EL01)
# 
# EL01sd <- 
#   dist(EL01s)
# 
# EL01sdh <- 
#   hclust(EL01sd, method="ward.D2")
# 
# plot(EL01sdh)
# 
# EL01sdh$order
# 
# rownames(EL01) <- 
#   1:58
# 
# df <- 
#   as.data.frame(as.table(EL01)) |> 
#   mutate(Freq=factor(Freq, levels=c(1, 0), labels=c("Early", "Late")))
# 
# ggplot(df, aes(x=Var2, y=Var1, fill=Freq)) + 
#   geom_tile() +
#   scale_fill_discrete(name=NULL) +
#   labs(x="SNP", y="Haplotype") +
#   theme(axis.text.x=element_text(angle=270, hjust=0))
# 
