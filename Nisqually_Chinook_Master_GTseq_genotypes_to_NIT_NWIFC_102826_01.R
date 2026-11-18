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

odat <-
  dat

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
    locus_count=count,
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
  .name_repair=name_repair,
  col_types="text")

dat_2025_Nisqually_CH_Genetics <- 
  read_excel(
    path="../data-from-Sayre/2025 Nisqually CH Genetics.xlsx",
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
    sample_datetime=as.numeric(sample_datetime) |> as.Date(origin="1899-12-30"),
    date=as.numeric(date) |> as.Date(origin="1899-12-30"),
    delivery_date=as.numeric(delivery_date) |> as.Date(origin="1899-12-30"),
    trucked_date=as.numeric(trucked_date) |> as.Date(origin="1899-12-30"),
    j_t_rec_date=as.numeric(j_t_rec_date) |> as.Date(origin="1899-12-30"))

dat_bio <- 
  relocate(
    dat_bio, sample_datetime, date, delivery_date, trucked_date, j_t_rec_date, 
    yearjtcolor_code, sampling_year, week_number, time, day_month, 
    .before=year)

dat_bio <- 
  relocate(
    dat_bio, 
    wdfw_id, match, code, nit_gen_sample_id_code, other_useful_information, 
    nit_sample_name, dna_number, genetic_label, genetic_sample, comment, 
    wdfw_genetic_sample_code_id, gs_number, gs_status, gs_comment, 
    scale_sample_, gilbert_rich, age, life_stage, length, length_cm_,
    .after=management_week)

dat_bio <- 
  relocate(
    dat_bio,
    sex, sex_1, sex_2, sex_id, 
    .after=fork_length_cm)

dat_bio <- 
  relocate(
    dat_bio,
    ad_clip, cwt, source_hatchery, hatchery_or_natural, origin, 
    .after=mark_types)

dat_bio <- 
  mutate(
    dat_bio, 
    sample_date=if_else(!is.na(sample_datetime), true=sample_datetime, false=date),
    sample_date=if_else(is.na(sample_date), true=trucked_date, false=sample_date),
    sample_date=if_else(is.na(sample_date), true=j_t_rec_date, false=sample_date), 
    sample_year=if_else(!is.na(sample_date), true=year(sample_date), false=as.double(year))) |>
  relocate(sample_year, sample_date, .before=sample_datetime) |>
  select(-sample_datetime, -date, -year, -sampling_year, -week_number, -time, -code, )


dat_bio <- 
  mutate(
    dat_bio, 
    sample_year=if_else(
      condition=is.na(sample_year), 
      true=as.double(str_extract(wdfw_id, "^\\d{2}"))+2000, 
      false=sample_year))

dat_bio <-
  filter(dat_bio, !is.na(sample_year)) |>
  arrange(sample_year, sample_date)

dat_bio <-
  filter(dat_bio, !gs_status %in% "ambiguous source")

dat_bio <-
  mutate(
    dat_bio, 
    dna_number=as.numeric(dna_number), 
    genetic_sample=as.numeric(genetic_sample), 
    gs_number=as.numeric(gs_number), 
    jaw_tag_number=as.numeric(jaw_tag_number), 
    length_cm=as.numeric(length), 
    fork_length_cm=as.numeric(fork_length_cm), 
    length_cm_=as.numeric(length_cm_))

dat_bio <-
  select(dat_bio, -day_month, -management_week, -spcode, -species, -speciesquantity)

dat_bio <-
  mutate(
    dat_bio, 
    length_cm=if_else(condition=is.na(length_cm), true=fork_length_cm, false=length_cm), 
    length_cm=if_else(condition=is.na(length_cm), true=length_cm_, false=length_cm)) |>
  relocate(length_cm, .before=length) |>
  select(-fork_length_cm, -length_cm_)


table(dat_bio$sex, useNA="ifany")

dat_bio <-
  mutate(dat_bio, lh=ifelse(test=sex=="Jack", yes="jack", no=NA))

dat_bio$sex <- 
  replace_values(
    dat_bio$sex, 
    c("1 Male", "Jack", "Male") ~ "M", 
    c("2 Female", "Female") ~ "F", 
    "Unknown" ~ NA)

table(dat_bio$sex_1, useNA="ifany")

dat_bio$sex_1 <- 
  replace_values(
    dat_bio$sex_1, 
    "1" ~ "M", 
    "2" ~ "F")

table(dat_bio$sex_id, useNA="ifany")

dat_bio$sex_id <- 
  replace_values(
    dat_bio$sex_id, 
    "no data" ~ NA)

table(dat_bio$sex_type_description, useNA="ifany")

dat_bio$sex_type_description <- 
  replace_values(
    dat_bio$sex_type_description, 
    "Unknown" ~ NA)

dat_bio <-
  mutate(dat_bio, lh=ifelse(test=sex_type_description=="Precocious Male", yes="jack", no=lh))


dat_bio <-
  mutate(
    dat_bio, 
    sex=if_else(condition=is.na(sex), true=sex_1, false=sex), 
    sex=if_else(condition=is.na(sex), true=sex_2, false=sex), 
    sex=if_else(condition=is.na(sex), true=sex_id, false=sex), 
    sex=if_else(condition=is.na(sex), true=sex_type_description, false=sex)) |>
  select(-sex_1, -sex_2, -sex_id)


table(dat_bio$sex, useNA="ifany")

dat_bio$sex <- 
  replace_values(
    dat_bio$sex, 
    c("Male", "Precocious Male") ~ "M", 
    "Female" ~ "F")


table(dat_bio$mark_types, useNA="ifany")

dat_bio$mark_types <- 
  replace_values(
    dat_bio$mark_types, 
    c("Ad Mark", "Adipose Clip, Unknown Coded Wire Tag") ~ "Adipose Clip",
    "Adipose Clip, Coded Wire Tag" ~ "Adipose Clip, CWT",
    "No Adipose Clip, Coded Wire Tag" ~ "CWT",
    c("Unknown", "Unknown Adipose Clip", "Unknown Adipose Clip, Unknown Coded Wire Tag") ~ NA)


table(dat_bio$ad_clip, useNA="ifany")

dat_bio$ad_clip <- 
  replace_values(
    dat_bio$ad_clip, 
    c("Clip Or CWT", "Marked (fin clip)", "Yes") ~ "TRUE",
    c("No", "Unmarked (no fin clip)") ~ "FALSE",
    c("unknown", "Unknown") ~ NA) |> 
  as.logical()


table(dat_bio$cwt, useNA="ifany")

dat_bio$cwt <- 
  replace_values(
    dat_bio$cwt, 
    c("Clip Or CWT", "CWT", "Yes") ~ "TRUE",
    c("No", "No CWT") ~ "FALSE",
    c("Unknown", "Unknown/No Head") ~ NA) |> 
  as.logical()


table(dat_bio$source_hatchery, useNA="ifany")

dat_bio$source_hatchery <- 
  replace_values(
    dat_bio$source_hatchery, 
    c("Unknown") ~ NA)


table(dat_bio$hatchery_or_natural, useNA="ifany")

table(dat_bio$origin, useNA="ifany")

dat_bio$origin <- 
  replace_values(
    dat_bio$origin, 
    "H" ~ "Hatchery", 
    c("?", "U", "Unknown") ~ NA)


dat_bio <- 
  mutate(
    dat_bio, 
    ad_clip=if_else(mark_types %in% c("Adipose Clip", "Adipose Clip, CWT"), true=TRUE, false=ad_clip),
    cwt=if_else(mark_types %in% c("Adipose Clip, CWT", "CWT"), true=TRUE, false=cwt))

dat_bio <- 
  mutate(
    dat_bio, 
    origin=if_else(is.na(origin), true=hatchery_or_natural, false=origin),
    origin=if_else(source_hatchery %in% "Clear Creek Hatchery", true="Hatchery", false=origin)) |>
  select(-hatchery_or_natural)

# 2012
weir2012 <- 
  filter(dat_bio, spreadsheet=="genetics sample numbers checking v2$weir.xlsx", sample_year==2012) |>
  arrange(genetic_sample)

table(dat$green_db_collection_name)

weir2012$wdfw_id <- 
  filter(dat, green_db_collection_name=="2012 Nisqually Weir Chinook") |> 
  pull(id)

sgs2012 <- 
  filter(dat_bio, spreadsheet=="genetics sample numbers checking v2$spawning ground.xlsx", sample_year==2012) |>
  arrange(gs_number)

sgs2012$wdfw_id <- 
  filter(dat, green_db_collection_name=="2012 Nisqually Spawning Grounds Chinook") |> 
  pull(id)


# 2013
#######################################################################################################
# weir incomplete
#######################################################################################################
weir2013 <- 
  filter(dat_bio, spreadsheet=="genetics sample numbers checking v2$weir.xlsx", sample_year==2013) |>
  filter(genetic_sample!=2) |>
  arrange(genetic_sample)

wd2013 <- 
  filter(dat, green_db_collection_name=="2013  Nisqually Weir Chinook") |> 
  select(1:2, 5, 7)

weir2013$sex[c(1, 3:31)]
wd2013$sex[1:30]

wd2013 %>% 
  print(n=nrow(.))

#######################################################################################################
# weir incomplete
#######################################################################################################

sgs2013 <- 
  filter(dat_bio, spreadsheet=="genetics sample numbers checking v2$spawning ground.xlsx", sample_year==2013) |>
  arrange(gs_number, comment)


sgs2013$wdfw_id[1:44] <- 
  filter(dat, green_db_collection_name=="2013 Nisqually Spawning Ground Survey Chinook") |> 
  pull(id)

sgs2013$comment <-
  gsub(pattern="13IL4", replacement="13IL004", sgs2013$comment) |> 
  gsub(pattern="13l00038", replacement="13IL0038")

sgs2013 <- 
  mutate(sgs2013, wdfw_id=if_else(is.na(wdfw_id), true=comment, false=wdfw_id))

# 2014
#######################################################################################################
# incomplete
#######################################################################################################
weir2014 <- 
  filter(dat_bio, spreadsheet=="genetics sample numbers checking v2$weir.xlsx", sample_year==2014) |>
  # filter(!release_condition %in% c("5 No Signs Of Life (Dead)", "Not Released")) |>
  # filter(!is.na(genetic_sample)) |>
  # filter(!is.na(jaw_tag_number)) |>
  # filter(grepl(pattern="^y", gs_status)) |>
  arrange(genetic_sample, sample_date)

weir2014$comment[138:144]

weir2014$genetic_sample[138:144] <-
  c(388:389, 391, 394, 390, 392:393)

weir2014 <- 
  arrange(weir2014, genetic_sample)


weir2014 <- 
  slice(weir2014, -15, -31, -40, -46, -138, -140, -144)

table(weir2014$release_condition, weir2014$release_action, useNA="ifany")
table(weir2014$release_condition, useNA="ifany")

# weir2014$scale_sample_ <- 
#   factor(weir2014$scale_sample_, levels=c("1-1", "1-2", "1-3", "1-4", "1-5", "1-6", "1-7", "1-8", "1-9", "1-10", "1-11", "1-12", "2-1", "2-2", "2-3", "2-4", "2-5", "2-6", "3-1", "3-2", "3-3", "3-4", "3-5", "3-6", "3-7", "3-8", "4-1", "4-2", "4-3", "4-4", "4-5", "4-6", "5-1", "5-2", "6-1", "6-2", "6-3", "6-4", "6-5", "6-6", "6-7", "6-8", "6-9", "6-10", "6-11", "9-1", "9-2", "9-3"), ordered=TRUE)

# weir2014 <- 
#   arrange(weir2014, sample_date, jaw_tag_number, genetic_sample)

weir2014$sex[1:20]
filter(dat, green_db_collection_name=="2014 Nisqually Weir Chinook")$sex[1:20]


# weir2014$wdfw_id <-
wd2014 <- 
  filter(dat, green_db_collection_name=="2014 Nisqually Weir Chinook") |> 
  select(1:2, 5, 7)

wd2014 %>% 
  print(n=nrow(.))


#######################################################################################################
# incomplete
#######################################################################################################

# no 2014 SGS bio data or genotypes
filter(dat_bio, spreadsheet=="genetics sample numbers checking v2$spawning ground.xlsx", sample_year==2014) |>
  arrange(gs_number, comment)


# 2015
# all genotypes from juveniles
# no bio data
filter(dat, grepl(pattern="^15", id)) |> 
  select(1:2, 5) %>% 
  print(n=nrow(.))


# 2016
# no genotypes or bio data
filter(dat, grepl(pattern="^16", id)) |> 
  select(1:2) %>% 
  print(n=nrow(.))


# 2017
a2017 <- 
  filter(dat_bio, sample_year==2017) |>
  arrange(nit_gen_sample_id_code)

a2017$nit_gen_sample_id_code[1714] <- "3800"

a2017$nit_gen_sample_id_code[c(1:17, 48:1942)] <-
  paste0("GEN-HRY->SG '17 # ", a2017$nit_gen_sample_id_code[c(1:17, 48:1942)])

a2017 <- 
  mutate(
    a2017, 
    nit_gen_sample_id_code=if_else(
      condition=is.na(nit_gen_sample_id_code), 
      true=wdfw_genetic_sample_code_id, 
      false=nit_gen_sample_id_code))


d2017 <- 
  filter(dat, grepl(pattern="^17", id)) |> 
  select(1:2)


a2017$nit_gen_sample_id_code <- 
  gsub(pattern="11KS75", replacement="11KS 75", a2017$nit_gen_sample_id_code) |>
  gsub(pattern="14mm 76", replacement="14MM 76") |>
  gsub(pattern="24MM 78", replacement="14MM 78")



a2017 <- 
  left_join(x=a2017[-8], y=rename(d2017, wdfw_id=id), by=join_by(nit_gen_sample_id_code==alt_id)) |> 
  relocate(wdfw_id, .before=match)


a2017$nit_gen_sample_id_code[18:45]

a2017$wdfw_id[18:45] <-
  d2017$id[87:114]

anti_join(x=d2017, y=a2017, by=join_by(id==wdfw_id)) %>% 
  print(n=nrow(.))



# 2018
a2018 <- 
  filter(dat_bio, sample_year==2018) |> 
  # mutate(nit_gen_sample_id_code=as.numeric(nit_gen_sample_id_code)) |>
  arrange(wdfw_id, nit_gen_sample_id_code)


d2018 <- 
  filter(dat, grepl(pattern="^18", id)) |> 
  select(1:2, 5, 7)

d2018 %>% 
  print(n=nrow(.))

# 2019
a2019 <- 
  filter(dat_bio, sample_year==2019)

d2019 <- 
  filter(dat, grepl(pattern="^19", id)) |> 
  select(1:2)

a2019 <- 
  mutate(
    a2019, 
    nit_gen_sample_id_code=if_else(
      condition=is.na(nit_gen_sample_id_code), 
      true=wdfw_genetic_sample_code_id, 
      false=nit_gen_sample_id_code))


table(a2019$wdfw_genetic_sample_code_id)

anti_join(x=d2019, y=a2019, by=join_by(alt_id==nit_gen_sample_id_code)) %>% 
  print(n=nrow(.))

a2019$nit_gen_sample_id_code <-
  gsub(pattern="GEN_SpnGnd_19- ", replacement="GEN_SpnGnd_19-", a2019$nit_gen_sample_id_code) |>
  gsub(pattern="MASH100219", replacement="MASH 100219")


a2019$nit_gen_sample_id_code[3263:3267]
a2019$wdfw_id[3263:3267] <-
  c("19QU0003", "19QU0004", "19QU0005", "19QU0006", "19QU0007")


a2019$nit_gen_sample_id_code[3388:3399]
a2019$wdfw_id[3388:3399] <-
  c("19QU0008", "19QU0009", "19QU0010", "19QU0011", "19QU0012", "19QU0013", 
    "19QU0014", "19QU0015", "19QU0016", "19QU0017", "19QU0018", "19QU0019")


a2019$nit_gen_sample_id_code[3400]
a2019$wdfw_id[3400] <- 
  "19QU0001"


a2019$nit_gen_sample_id_code[3401:3402]
a2019$wdfw_id[3401:3402] <- 
  c("19QU0025", "19QU0026")



a2019$nit_gen_sample_id_code[3146]
a2019$wdfw_id[3146] <- 
  "19QU0002"

a2019$nit_gen_sample_id_code[3403:3407]
a2019$wdfw_id[3403:3407] <- 
  c("19QU0020", "19QU0021", "19QU0022", "19QU0023", "19QU0024")

a2019$nit_gen_sample_id_code[3409]
a2019$wdfw_id[3409] <- 
  "19QU0027"


a2019_na <- 
  filter(a2019, is.na(wdfw_id))

a2019_no <- 
  filter(a2019, !is.na(wdfw_id))

a2019_na <- 
  left_join(x=a2019_na[-8], y=rename(d2019, wdfw_id=id), by=join_by(nit_gen_sample_id_code==alt_id)) |> 
  relocate(wdfw_id, .before=match)

a2019 <- 
  bind_rows(a2019_na, a2019_no) |> 
  arrange(wdfw_id)




# 2021
a2021 <- 
  filter(dat_bio, sample_year==2021)

d2021 <- 
  filter(dat, grepl(pattern="^21", id)) |> 
  select(1:2)


a2021$genetic_label <-
  gsub(pattern="GEN_SG_21_", replacement="Gen_SG_21-", a2021$genetic_label)

a2021 <- 
  left_join(x=a2021[-8], y=rename(d2021, wdfw_id=id), by=join_by(genetic_label==alt_id)) |> 
  relocate(wdfw_id, .before=match)



# 2022
a2022 <- 
  filter(dat_bio, sample_year==2022)

d2022 <- 
  filter(dat, grepl(pattern="^22", id)) |> 
  select(1:2)

a2022$genetic_label <- 
  gsub(pattern=" ", replacement="", a2022$genetic_label)

a2022 <- 
  distinct(a2022, genetic_label, .keep_all=TRUE)

a2022 <- 
  left_join(x=a2022[-8], y=rename(d2022, wdfw_id=id), by=join_by(genetic_label==alt_id)) |> 
  relocate(wdfw_id, .before=match)


# 2023
a2023 <- 
  filter(dat_bio, sample_year==2023)

d2023 <- 
  filter(dat, grepl(pattern="^23", id)) |> 
  select(1:2)

a2023$genetic_label <- 
  gsub(pattern=" ", replacement="", a2023$genetic_label) |>
  gsub(pattern="CH.Gen.Trad2023", replacement="CH.Gten.Trad2023", fixed=TRUE)

a2023 <- 
  distinct(a2023, genetic_label, .keep_all=TRUE)

"CH.Gten.Trad2023-8"
'CH.Gen.Trad2023-8'

a2023 <- 
  left_join(x=a2023[-8], y=rename(d2023, wdfw_id=id), by=join_by(genetic_label==alt_id)) |> 
  relocate(wdfw_id, .before=match)

a2023$genetic_label[355]
  
a2023$wdfw_id[355] <- 
  "23NF0052"

a2023$genetic_label[382]

a2023$wdfw_id[382] <- 
  "23NF0056"

# 2024
# no bio data
filter(dat, grepl(pattern="^24", id)) |> 
  select(1:2) %>% 
  print(n=nrow(.))

# 2025
a2025 <- 
  filter(dat_bio, sample_year==2025) |>
  arrange(nit_sample_name)

a2025$nit_sample_name <- 
  gsub(pattern=" ", replacement="", a2025$nit_sample_name) |>
  gsub(pattern="SpnGnd25", replacement="SpnGnd_25")

d2025 <- 
  filter(dat, grepl(pattern="^25", id)) |> 
  select(1:2, 5, 7)

a2025 <- 
  left_join(
    x=a2025[-8], 
    y=select(d2025, wdfw_id=id, alt_id), 
    by=join_by(nit_sample_name==alt_id)) |> 
  relocate(wdfw_id, .before=match)

anti_join(x=d2025, y=a2025, by=join_by(alt_id==nit_sample_name)) %>% 
  print(n=nrow(.))

a2025_na <- 
  filter(a2025, is.na(wdfw_id))

################################################################################
# looks like the n=49 samples in a2025_na were not genotyped
################################################################################

a2025_na_SpnGnd <- 
  filter(a2025_na, grepl(pattern="SpnGnd", nit_sample_name)) |> 
  arrange(dna_number)

a2025_na_CIR <- 
  filter(a2025_na, grepl(pattern="CIR", nit_sample_name)) |> 
  arrange(dna_number)

a2025_na_Trad <- 
  filter(a2025_na, grepl(pattern="Trad", nit_sample_name)) |> 
  arrange(dna_number)

################################################################################
# 
################################################################################


table(dat$life_history, useNA="ifany")

table(dat$green_db_collection_name, useNA="ifany") |> 
  enframe() |> 
  mutate(value=as.numeric(value)) %>% 
  print(n=nrow(.))


# 2013 & 2014 weir bio data couldn't be matched to genetic data.
# I can match #'s of samples but can't get the row orders set, 
# so phenotypic sex wouldn't match genotypic sex.

# 2015 no bio data (all genetic data is from juveniles)

# 2016 no genotypes or bio data

# 2024 no bio data

# 2019 only year with bio data for smolts


bio <-
  bind_rows(weir2012, sgs2012, weir2013, sgs2013, weir2014, a2017, a2018, a2019, a2021, a2022, a2023, a2025) |>
  filter(!is.na(wdfw_id)) |>
  arrange(wdfw_id) |> 
  distinct(wdfw_id, .keep_all=TRUE)

8444/8772

bio[2787, 'length_cm'] <- 
  NA

bio <- 
  mutate(bio, wdfw_code=str_extract(wdfw_id, "^\\w{4}"))

gen <- 
  filter(dat, grepl(pattern="^12|^13|^14|^15|^16|^17|^18|^19|^20|^21|^22|^23|^24|^25", id)) |> 
  select(1:3, 5, 7:8) |>
  arrange(id)

gen$alt_id[gen$alt_id=="NA"] <- 
  NA

gen <- 
  mutate(gen, code=str_extract(id, "^\\w{4}"))


left_join(
  x=group_by(gen, code, green_db_collection_name) |> count() |> ungroup(), 
  y=group_by(bio, wdfw_code) |> count() |> ungroup(), 
  by=join_by(code==wdfw_code), 
  suffix=c("_genetic", "_biodata"))  %>% 
  write_csv("../data/genetic-biodata-summary.csv", na="0")

8444/19489

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
