library(sequoia)
library(ASRgenomics)
library(haplo.stats)
library(readxl)
library(tidyverse)

setwd("~/Google Drive/My Drive/Nisqually/Nisqually/")

excel_sheets(path="../Nisqually_Chinook_Master_GTseq_genotypes_to_NIT_NWIFC_102826.xlsx")

dat <- 
  read_excel(
    path="../Nisqually_Chinook_Master_GTseq_genotypes_to_NIT_NWIFC_102826.xlsx", 
    sheet="genotypes", 
    col_types="text")

glimpse(dat)

names(dat) <- 
  names(dat) |> 
  gsub(pattern="-", replacement="_") |> 
  gsub(pattern=" ", replacement="_") |> 
  gsub(pattern="'", replacement="_") |> 
  gsub(pattern=".", replacement="_", fixed=TRUE) |> 
  gsub(pattern="%", replacement="percent", fixed=TRUE)


table(dat$code, useNA="ifany")

table(dat$Green_DB_collection_name, useNA="ifany")

table(dat$genotyping_year, useNA="ifany")

table(dat$Life_history, useNA="ifany")

table(dat$species_ID, useNA="ifany")

dat$species_ID <- 
  replace_values(
    dat$species_ID, 
    "no data" ~ NA, 
    "No data" ~ NA, 
    "non-Chinook" ~ "not Chinook",
    "Not Chinook" ~ "not Chinook")

table(dat$speciesID, useNA="ifany")

dat$speciesID <- 
  replace_values(
    dat$speciesID, 
    "no data" ~ NA)

table(dat$ploidy, useNA="ifany")

dat$ploidy <- 
  replace_values(
    dat$ploidy, 
    "no data" ~ NA, 
    "U" ~ NA) # assuming "U" stands for "Unknown"

table(dat$Sex_ID, useNA="ifany")

dat$Sex_ID <- 
  replace_values(
    dat$Sex_ID, 
    "no data" ~ NA, 
    "No data" ~ NA)

table(dat$QC, useNA="ifany")

table(dat$QC_status, useNA="ifany")

lookup1 <- 
  read_tsv(file="Green_DB_collection_name__origin.tsv")

dat <- 
  mutate(
    dat, 
    origin=recode_values(
      Green_DB_collection_name, 
      from=lookup1$Green_DB_collection_name, 
      to=lookup1$origin),
    lifestage=Life_history) |>
  mutate(
    origin=if_else(
      condition=grepl(pattern="Unmarked", Alt_ID), 
      true="natural", 
      false=origin), 
    lifestage=if_else(
      condition=grepl(pattern="Jack", Alt_ID), 
      true="jack", 
      false=lifestage)) |>
  mutate(
    origin=if_else(
      condition=grepl(pattern="AdiposeClip", Alt_ID), 
      true="hatchery", 
      false=origin))


distinct(dat, Green_DB_collection_name, origin, lifestage) %>% 
  print(n=nrow(.))

dat <- 
  left_join(x=dat, y=lookup1[-2])

table(dat$sample_type, dat$lifestage, useNA="ifany")
table(dat$Green_DB_collection_name, useNA="ifany")

gt331 <- 
  filter(
    dat, 
    !QC_status %in% names(table(dat$QC_status))[-8], 
    !species_ID %in% "not Chinook", 
    !speciesID %in% c("Oke", "Oki"),
    !sample_type %in% c("baseline", "return")) |> 
  select(
    ID=USE_sample_name, 
    Green_DB_collection_name, 
    origin, 
    sample_type, 
    lifestage, 
    IFI, 
    Sex=Sex_ID, 
    28:358) |>
  mutate(IFI=as.numeric(IFI))

glimpse(gt331)

table(gt331$sample_type, gt331$lifestage)


gt331 %>% 
  select(where(~any(grepl("-", .x, fixed=TRUE)))) %>% 
  glimpse()

gt331$Ots_128757_61R <- 
  gsub(pattern="-", replacement="G", x=gt331$Ots_128757_61R)

gt331$Ots_GDH_81x <- 
  gsub(pattern="-", replacement="T", x=gt331$Ots_GDH_81x)

gt331$Ots_Ots311_101x <- 
  gsub(pattern="-", replacement="G", x=gt331$Ots_Ots311_101x)

gt331$Ots_RFC2_558 <- 
  gsub(pattern="-", replacement="G", x=gt331$Ots_RFC2_558)

gt331$Ots_U2305_63 <- 
  gsub(pattern="-", replacement="C", x=gt331$Ots_U2305_63)

gt331$Ots_pigh_105 <- 
  gsub(pattern="-", replacement="G", x=gt331$Ots_pigh_105)

gt331$Ots_u07_17_373 <- 
  gsub(pattern="-", replacement="G", x=gt331$Ots_u07_17_373)

gt331$Ots_vatf_251 <- 
  gsub(pattern="-", replacement="A", x=gt331$Ots_vatf_251)

gt331 <- 
  relocate(
    gt331, 
    Ots28_11023212, Ots28_11025336, Ots28_11033282, Ots28_11062192, Ots37124_12267397, 
    Ots37124_12272852, Ots28_11070757, Ots28_11071377, Ots28_11072994, Ots28_11073102, 
    Ots37124_12277401, Ots28_11073668, Ots28_11075348, Ots28_11075712, Ots28_11077016, 
    Ots28_11077172, Ots37124_12281207, Ots28_11077576, Ots28_11095755, Ots37124_12310649, 
    Ots28_11143508, Ots28_11160599, Ots28_11164637, Ots28_11186543, Ots28_11201129, 
    Ots28_11202190, Ots28_11202400, Ots28_11202863, Ots28_11205423, Ots28_11205993, 
    Ots28_11206740, Ots28_11207428, .before=Ots28_11210919)

lh331 <- 
  select(gt331, lifestage, sample_type, origin, ID1=ID, Sex, BirthYear=ID) |> 
  rename(ID=ID1) |>
  mutate(BirthYear=str_extract(BirthYear, "^\\d{2}") |> as.integer()) |>
  mutate(BirthYear=if_else(lifestage=="juvenile", true=BirthYear+1999, false=NA))

table(lh331$sample_type, lh331$BirthYear)


gt33 <- 
  select(gt331, ID, Ots28_11023212, Ots28_11025336, Ots28_11033282, Ots28_11062192, Ots37124_12267397, Ots37124_12272852, Ots28_11070757, Ots28_11071377, Ots28_11072994, Ots28_11073102, Ots37124_12277401, Ots28_11073668, Ots28_11075348, Ots28_11075712, Ots28_11077016, Ots28_11077172, Ots37124_12281207, Ots28_11077576, Ots28_11095755, Ots37124_12310649, Ots28_11143508, Ots28_11160599, Ots28_11164637, Ots28_11186543, Ots28_11201129, Ots28_11202190, Ots28_11202400, Ots28_11202863, Ots28_11205423, Ots28_11205993, Ots28_11206740, Ots28_11207428, Ots28_11210919) |> 
  column_to_rownames(var="ID") |>
  filter(!if_all(everything(), is.na))

glimpse(gt33)

gm33 <- 
  snp.recode(
    M=as.matrix(gt33), 
    recoding="ATGCto012")

snpstats33 <- 
  SnpStats(gm33$Mrecode, Plot=FALSE)

gm33map <- 
  mutate(gm33$map, snp=names(gt33))


gt33 <- 
  separate_wider_position(gt33, cols=everything(), widths=c(a1=1, a2=1), names_sep="_") |>
  mutate(ID=rownames(gt33)) |>
  column_to_rownames(var="ID") |>
  filter(!if_any(everything(), is.na)) |>
  as.matrix()


ht33 <- 
  haplo.em(geno=gt33, locus.label=tail(names(gt331), n=33))

ht33haps <- 
  ht33$haplotype |> 
  as_tibble()

EL01 <- 
  mutate(
    ht33haps, 
    Ots28_11023212=recode_values(Ots28_11023212, "G" ~ 0, "A" ~ 1),
    Ots28_11025336=recode_values(Ots28_11025336, "C" ~ 0, "A" ~ 1),
    Ots28_11033282=recode_values(Ots28_11033282, "A" ~ 0, "G" ~ 1),
    Ots28_11062192=recode_values(Ots28_11062192, "G" ~ 0, "C" ~ 1),
    Ots37124_12267397=recode_values(Ots37124_12267397, "C" ~ 0, "T" ~ 1),
    Ots37124_12272852=recode_values(Ots37124_12272852, "C" ~ 0, "T" ~ 1),
    Ots28_11070757=recode_values(Ots28_11070757, "G" ~ 0, "A" ~ 1),
    Ots28_11071377=recode_values(Ots28_11071377, "C" ~ 0, "T" ~ 1),
    Ots28_11072994=recode_values(Ots28_11072994, "T" ~ 0, "C" ~ 1),
    Ots28_11073102=recode_values(Ots28_11073102, "A" ~ 0, "T" ~ 1),
    Ots37124_12277401=recode_values(Ots37124_12277401, "A" ~ 0, "T" ~ 1),
    Ots28_11073668=recode_values(Ots28_11073668, "A" ~ 0, "T" ~ 1),
    Ots28_11075348=recode_values(Ots28_11075348, "A" ~ 0, "G" ~ 1),
    Ots28_11075712=recode_values(Ots28_11075712, "T" ~ 0, "C" ~ 1),
    Ots28_11077016=recode_values(Ots28_11077016, "T" ~ 0, "C" ~ 1),
    Ots28_11077172=recode_values(Ots28_11077172, "A" ~ 0, "G" ~ 1),
    Ots37124_12281207=recode_values(Ots37124_12281207, "A" ~ 0, "T" ~ 1),
    Ots28_11077576=recode_values(Ots28_11077576, "G" ~ 0, "A" ~ 1),
    Ots28_11095755=recode_values(Ots28_11095755, "A" ~ 0, "T" ~ 1),
    Ots37124_12310649=recode_values(Ots37124_12310649, "T" ~ 0, "A" ~ 1),
    Ots28_11143508=recode_values(Ots28_11143508, "G" ~ 0, "A" ~ 1),
    Ots28_11160599=recode_values(Ots28_11160599, "G" ~ 0, "T" ~ 1),
    Ots28_11164637=recode_values(Ots28_11164637, "C" ~ 0, "A" ~ 1),
    Ots28_11186543=recode_values(Ots28_11186543, "A" ~ 0, "T" ~ 1),
    Ots28_11201129=recode_values(Ots28_11201129, "T" ~ 0, "G" ~ 1),
    Ots28_11202190=recode_values(Ots28_11202190, "C" ~ 0, "T" ~ 1),
    Ots28_11202400=recode_values(Ots28_11202400, "C" ~ 0, "T" ~ 1),
    Ots28_11202863=recode_values(Ots28_11202863, "C" ~ 0, "A" ~ 1),
    Ots28_11205423=recode_values(Ots28_11205423, "G" ~ 0, "A" ~ 1),
    Ots28_11205993=recode_values(Ots28_11205993, "T" ~ 0, "C" ~ 1),
    Ots28_11206740=recode_values(Ots28_11206740, "C" ~ 0, "T" ~ 1),
    Ots28_11207428=recode_values(Ots28_11207428, "T" ~ 0, "G" ~ 1),
    Ots28_11210919=recode_values(Ots28_11210919, "C" ~ 0, "T" ~ 1))


EL01 <- 
  as.matrix(EL01)

EL01s <- 
  scale(EL01)

EL01sd <- 
  dist(EL01s)

EL01sdh <- 
  hclust(EL01sd, method="ward.D2")

plot(EL01sdh)

EL01sdh$order

rownames(EL01) <- 
  1:58

df <- 
  as.data.frame(as.table(EL01)) |> 
  mutate(Freq=factor(Freq, levels=c(1, 0), labels=c("Early", "Late")))

ggplot(df, aes(x=Var2, y=Var1, fill=Freq)) + 
  geom_tile() +
  scale_fill_discrete(name=NULL) +
  labs(x="SNP", y="Haplotype") +
  theme(axis.text.x=element_text(angle=270, hjust=0))

