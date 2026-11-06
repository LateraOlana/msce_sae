# Raw inputs

The raw inputs are not part of the repository. Place them in this folder as
listed below (or anywhere else and point `MCV1_DATA_DIR` or `R/config_local.R`
at it; see `R/00_config.R`). Every path is read from `cfg$paths`.

| folder | file(s) | source | used by |
|---|---|---|---|
| `survey/` | `mini_edhs_2019_children.csv` | 2019 Ethiopia Mini Demographic and Health Survey, children's recode (KR), DHS Program, after registration at https://dhsprogram.com | stage 02 |
| `dhs_gps/` | `ETGE81FL.shp` and its side files | 2019 mini-EDHS cluster locations (geographic data set ETGE81FL), DHS Program | stage 02 |
| `boundaries/` | `Eth_Wereda_Bofed_CSA_2021.shp` and its side files | 2021 administrative boundaries of Ethiopia at woreda (third) level, Central Statistical Agency, with the fields `ADM1_EN`, `ADM1_PCODE`, `ADM2_EN`, `ADM2_PCODE`, `ADM3_EN`, `ADM3_PCODE` | stage 02 |
| `ghsl/` | `GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0_<tile>.tif`, `GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0_<tile>.tif` for the tiles `R8_C22`, `R8_C23`, `R9_C22`, `R9_C23` | Global Human Settlement Layer, release R2023A, epoch 2020, 1 km, Mollweide: https://human-settlement.emergency.copernicus.eu/download.php | stage 04 |
| `ghsl/` | the same four tiles of `GHS_SMOD_E2030_...` | only if `cfg$urban$epoch` is set to 2030 | stage 04 |
| `mbg/` | `agg_summaries_all_admin2.csv`, `agg_summaries_all.csv` | IHME model-based geostatistical estimates of vaccine coverage, aggregated to the 2021 second- and third-level units; columns `ADM2_PCODE` (or `ADM3_PCODE`), `year`, `vaccine`, `mean`, `lower`, `upper` | stage 06 (the comparison with the geostatistical estimates) |

## The survey extract

`mini_edhs_2019_children.csv` has one row per child under three years of age
and the columns below. They are standard DHS children's-recode variables; the
names are the ones the scripts use.

| column | content |
|---|---|
| `cluster` | DHS cluster number (`v001`), matching `DHSCLUST` in the GPS file |
| `householdID` | household identifier |
| `weight` | sampling weight of the child (`v005`, any constant scaling) |
| `urban` | 1 urban, 0 rural (`v025`) |
| `strata` | sampling stratum (`v023`; 21 strata, region by urban/rural) |
| `age_year` | age in completed years (the analysis keeps `age_year == 1`, that is 12 to 23 months) |
| `mcv1` | first dose of measles-containing vaccine received: 1 yes, 0 no, missing if unknown |
| `has_vacc_card` | `"seen"` when the interviewer saw a vaccination card, otherwise another value or missing |
| `mcv_dose_from_card` | number of measles doses recorded on the card (missing if no card) |
| `mcv_dose_from_recall` | number of measles doses reported by the caregiver (missing if not asked) |

DHS Program data may not be redistributed, and cluster coordinates may not be
published. Everything in `output/` is therefore aggregated to administrative
units.
