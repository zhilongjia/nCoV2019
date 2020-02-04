

clinial_raw <- readxl::read_xlsx("../documents/李明锟/clinical information.xlsx")
sample_rawmeta <- readxl::read_xlsx("../data/sample_meta-v2.xlsx")


sample_meta_v4 <- dplyr::full_join(sample_rawmeta, clinial_raw, by=)


readr::write_csv(sample_meta_v4, path="../data/sample_meta_v4.csv", na="")

