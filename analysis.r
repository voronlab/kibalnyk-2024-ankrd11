suppressPackageStartupMessages({
  library(Seurat)
  library(magrittr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(forcats)
})


mk_seurat <- function(dir_path) {
  meta <- paste(dir_path, "meta.csv", sep = "/") %>%
    read.csv(
      row.names = 1,
      numerals = "no.loss"
    ) %>%
    mutate(
      condition = ifelse(isKO == "True", "ncko", "nchet"),
      within.bb.heart = within.bb.heart == "True"
    ) %>%
    select(!isKO)

  so <- paste(dir_path, "counts.csv", sep = "/") %>%
    read.csv(row.names = 1, check.names = FALSE) %>%
    as.matrix() %>%
    CreateSeuratObject(assay = "vizgen") %>%
    AddMetaData(meta) %>%
    AddMetaData(
      paste(dir_path, "centroids.csv", sep = "/") %>%
        read.csv(
          row.names = 1,
          numerals = "no.loss"
        ) %>%
        mutate(y_num = as.numeric(x), x_num = as.numeric(y)) %>%
        mutate(y = y_num, x = x_num) %>%
        CreateFOV(
          key = paste(
            "spatial",
            meta[[1, "source"]] %>%
              str_remove_all("-") %>%
              str_to_lower(),
            "_",
            sep = ""
          ),
          type = "centroids",
          assay = "vizgen"
        ),
      meta[[1, "source"]]
    )

  so
}

sp <- function(so, algo = "leiden") {
  so %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(algorithm = algo) %>%
    RunUMAP(dims = 1:30)
}

strip_fov <- function(so) {
  Images(so) %>%
    reduce(
      function(so, x) {
        AddMetaData(so, NULL, x)
      },
      .init = so
    )
}

sample_dirs <- list(
  "AH59.5" = "202309211049_MsEmbryo-VS147-AH59-2-5-A_VMSC18002/region_1",
  "AH59.6" = "202309151103_MsEmbryo-VS147-AH59-4-6-C_VMSC17602_reanalysis/region_1",
  "AI04.2" = "20231204_Ankrd11_OFT/Region1",
  "AH50.1" = "202309181306_MsEmbryo-VS147-AH50-59-1-A_VMSC07101/region_0",
  "AH59.2" = "202309151104_MsEmbryo-VS147-AH59-2-5-B_VMSC16102/region_1",
  "AI13.3" = "20231204_Ankrd11_OFT/Region0"
)

so <- map(sample_dirs, compose(
  mk_seurat,
  partial(paste, "out", ... = , sep = "/")
)) %>%
  reduce(merge) %>%
  JoinLayers()

# ensure nchet grouping comes first in plots
so <- AddMetaData(
  so,
  mutate(
    so[["condition"]],
    condition = fct_relevel(condition, "nchet", "ncko")
  )
)

so_heart <- subset(so, within.bb.heart) %>% sp()
so <- AddMetaData(so, so_heart[["seurat_clusters"]])

saveRDS(so, "out/all_cells.rds")
saveRDS(so_heart, "out/heart_only.rds")
write.csv(select(so[[]], !orig.ident), "out/metadata.csv")
