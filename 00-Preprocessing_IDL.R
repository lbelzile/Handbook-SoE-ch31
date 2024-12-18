###########################################################
#####   Preprocessing of the englandwales database   ######
###########################################################
# This script shows how to create the 'englandwales' data
# from the file contianing all of the records
# https://www.supercentenarians.org//en/data-and-metadata/?idl_userDownloadFile=57
# and turn it into a data frame that is suitable
# for other analyses using the R package longevity.
setwd(this.path::here())
library(lubridate)
library(tidyverse)

# 1) Download the IDL data
# 2) Select relevant columns
# 3) Rename columns
# 4) Extract country ID from the IDNumber (first three characters)
# 5) Cast columns to factor or Date, whenever relevant

# Note: user must login with credentials to "supercentenarians.org"
#  in order to download the file. If it does not work, sign in
#  and download a copy, then modify the code to fetch a local version
url <- "https://www.supercentenarians.org/en/data-and-metadata/?idl_userDownloadFile=57"
rawdata <- read.table(
  file = url,
  header = TRUE,
  sep = ";"
)
# Replace url with absolute path of downloaded
# file if necessary
ew <- rawdata |>
  as_tibble() |>
  select(
    IDL_ID, AGEDAYS, AGEYEARS, SEX, BIRTH_DAY, BIRTH_MONTH, BIRTH_YEAR,
    DEATH_YEAR, DEATH_MONTH, DEATH_DAY, VALIDATION
  ) |>
  mutate(
    DDATE = paste0(DEATH_YEAR, "-", DEATH_MONTH, "-", DEATH_DAY),
    BDATE = paste0(BIRTH_YEAR, "-", BIRTH_MONTH, "-", BIRTH_DAY)
  ) |>
  rename(
    ndays = AGEDAYS,
    gender = SEX,
    bdate = BDATE,
    ddate = DDATE,
    ageyear = AGEYEARS,
    valid = VALIDATION
  ) |>
  mutate(country = factor(substr(IDL_ID, 1, 3))) |>
  filter(country == "EAW") |>
  select(!c(
    IDL_ID, DEATH_YEAR, DEATH_MONTH, DEATH_DAY,
    BIRTH_YEAR, BIRTH_MONTH, BIRTH_DAY
  )) |>
  # USA doesn't record dates, only years
  mutate(
    bdate = ymd(bdate),
    ddate = ymd(ddate),
    valid = factor(valid),
    gender = factor(gender)
  )

sampling_frame <- matrix(
  ncol = 4L,
  byrow = TRUE,
  data = c(
    "EAW", "110+", "01-01-1968", "31-12-2020",
    "EAW", "105-109", "28-12-1999", "31-12-2014"
  )
) |>
  as_tibble(.name_repair = "universal") |>
  rename(
    country = ...1,
    group = ...2,
    ldate = ...3,
    rdate = ...4
  ) |>
  mutate(
    country = factor(country),
    group = factor(group),
    ldate = dmy(ldate),
    rdate = dmy(rdate)
  )
sf_wide <- tidyr::pivot_wider(sampling_frame,
  names_from = "group",
  values_from = c("ldate", "rdate")
) |>
  rename(
    c1 = "ldate_110+",
    d1 = "ldate_105-109",
    c2 = "rdate_110+",
    d2 = "rdate_105-109"
  )
ew <- ew |>
  mutate(
    country = factor(country, exclude = NULL),
    group = factor(ifelse(ageyear < 110, "105-109", "110+"))
  ) |>
  left_join(y = sf_wide, by = c("country")) |>
  mutate(
    x110 = if_else(mday(bdate) == 29 & month(bdate) == 2,
      bdate + days(1) + years(110), # leap year and people born on February 29th
      bdate + years(110)
    ),
    x105 = if_else(mday(bdate) == 29 & month(bdate) == 2,
      bdate + days(1) + years(105),
      bdate + years(105)
    )
  ) |>
  select(!group) |>
  filter(if_else(ageyear >= 110,
    ddate >= c1 & ddate <= c2,
    ddate >= d1 & ddate <= d2
  )) |>
  arrange(country, ndays) |>
  mutate(gender = factor(ifelse(gender == "F", "female", "male")))

# Pivot the table and create columns with left and right truncation
# potentially double because of interval truncation
ew$ltrunc1 <- NA
ew$rtrunc1 <- NA
ew$ltrunc2 <- NA
ew$rtrunc2 <- NA
# This for loop is inefficient, but the code is more digestable
for (i in seq_along(ew$ndays)) {
  if (is.na(ew$d1[i])) {
    # For countries with no semisupercentenarian
    ew$ltrunc1[i] <- as.integer(max(ew$x110[i], ew$c1[i]) - ew$bdate[i])
    # [110+) collection first, person reached 110 after [105,110) period
  } else if (ew$x105[i] > ew$d2[i]) {
    ew$ltrunc1[i] <- as.integer(max(ew$x110[i], ew$c1[i]) - ew$bdate[i])
    # Person not observed before 110
  } else if (ew$c1[i] <= ew$d1[i] & ew$d1[i] < ew$x110[i]) {
    ew$ltrunc1[i] <- as.integer(max(ew$x105[i], ew$d1[i]) - ew$bdate[i])
    # Person reached 110 before [105,110) period
  } else if (ew$c1[i] <= ew$d1[i] & ew$x110[i] <= ew$d1[i]) {
    ew$ltrunc1[i] <- as.integer(max(ew$x110[i], ew$c1[i]) - ew$bdate[i])
    # [105, 110) period first, person reached 110 after that period
  } else if (ew$d1[i] < ew$c1[i] & ew$x110[i] >= ew$c1[i]) {
    ew$ltrunc1[i] <- as.integer(max(ew$x105[i], ew$d1[i]) - ew$bdate[i])
    # [105, 110) period first, but person was already above 110 before then
  } else if (ew$d1[i] < ew$c1[i] & ew$x105[i] < ew$c1[i]) {
    ew$ltrunc1[i] <- as.integer(max(ew$d1[i], ew$x105[i]) - ew$bdate[i])
    if (ew$x110[i] < ew$c1[i]) {
      # Discontinuous trajectory
      ew$ltrunc2[i] <- as.integer(ew$c1[i] - ew$bdate[i])
    }
  }
  # # Lifeline snaps at upper end, [105, 110) collection period ends before 110+
  # # person could have appeared in [105, 110) window
  # if(ew$d2[i] < ew$c2[i] & ew$x110[i] > ew$d2[i] & ew$x105[i] < ew$d2[i]){
  #   ew$ltrunc1[i] <- as.integer(max(ew$d1[i], ew$x105[i]) - ew$bdate[i])
  #   ew$ltrunc2[i] <- as.integer(ew$x110[i] - ew$bdate[i])
  # }
  # Right truncation
  if (is.na(ew$d2[i])) { # no semisupercentenarians
    ew$rtrunc1[i] <- as.integer(ew$c2[i] - ew$bdate[i])
  } else if (ew$d2[i] < ew$c2[i] & ew$x110[i] > ew$d2[i] & ew$x105[i] < ew$d2[i]) {
    # Lifeline snaps at upper end, [105, 110) collection period ends before 110+
    # person could have appeared in [105, 110) window
    ew$rtrunc1[i] <- as.integer(ew$d2[i] - ew$bdate[i])
    # ew$rtrunc2[i] <- as.integer(ew$c2[i] - ew$bdate[i])
  } else if (ew$c2[i] < ew$x110[i] & ew$c2[i] <= ew$d2[i]) {
    ew$rtrunc1[i] <- as.integer(min(ew$x110[i], ew$d2[i]) - ew$bdate[i])
  } else {
    ew$rtrunc1[i] <- as.integer(ew$c2[i] - ew$bdate[i])
  }
  # Double truncation
  if (!is.na(ew$d1[i]) && !is.na(ew$d2[i])) {
    if (ew$d1[i] < ew$c1[i] & ew$x110[i] > ew$d1[i] & ew$x110[i] < ew$c1[i]) {
      ew$ltrunc1[i] <- as.integer(max(ew$x105[i], ew$d1[i]) - ew$bdate[i])
      ew$ltrunc2[i] <- as.integer(ew$c1[i] - ew$bdate[i])
      ew$rtrunc1[i] <- as.integer(ew$x110[i] - ew$bdate[i])
      ew$rtrunc2[i] <- as.integer(ew$c2[i] - ew$bdate[i])
    } else if (ew$d2[i] < ew$c2[i] & ew$d2[i] < ew$x110[i] & ew$x105[i] < ew$d2[i]) {
      ew$ltrunc1[i] <- as.integer(max(ew$x105[i], ew$d1[i]) - ew$bdate[i])
      ew$rtrunc1[i] <- as.integer(ew$d2[i] - ew$bdate[i])
      ew$ltrunc2[i] <- as.integer(ew$x110[i] - ew$bdate[i])
      ew$rtrunc2[i] <- as.integer(ew$c2[i] - ew$bdate[i])
    }
  }
}
englandwales <- ew |>
  mutate(
    ltrunc1 = as.integer(ltrunc1),
    rtrunc1 = as.integer(rtrunc1),
    ltrunc2 = as.integer(ltrunc2),
    rtrunc2 = as.integer(rtrunc2)
  ) |>
  select(!c(c1, c2, d1, d2, x105, x110, country))
save(englandwales, file = "../Data/englandwales.rda")
