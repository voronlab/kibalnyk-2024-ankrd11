suppressPackageStartupMessages({
  library(Seurat)
  library(magrittr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(forcats)
})

SAMPLE_DIRS <- list(
  "AH59.5" = "202309211049_MsEmbryo-VS147-AH59-2-5-A_VMSC18002/region_1",
  "AH59.6" = "202309151103_MsEmbryo-VS147-AH59-4-6-C_VMSC17602_reanalysis/region_1",
  "AI04.2" = "20231204_Ankrd11_OFT/Region1",
  "AH50.1" = "202309181306_MsEmbryo-VS147-AH50-59-1-A_VMSC07101/region_0",
  "AH59.2" = "202309151104_MsEmbryo-VS147-AH59-2-5-B_VMSC16102/region_1",
  "AI13.3" = "20231204_Ankrd11_OFT/Region0"
)

SO <- readRDS("out/all_cells.rds")
SO_HEART <- readRDS("out/heart_only.rds")
SO_CNCC <- subset(SO_HEART, seurat_clusters %in% c(2, 3, 7))


DISCRETE_PALETTE <- "glasbey"
CONTINUOUS_PALETTE <- "plasma"
FONT_SIZE <- 20

get_heart_bb <- function(dir) {
  points <- paste(
    "out",
    dir,
    "bb.csv",
    sep = "/"
  ) %>%
    read.csv(
      row.names = 1
    ) %>%
    .[["bb.heart", "coords"]] %>%
    str_extract_all(
      .,
      "\\-?[0-9\\.]+, \\-?[0-9\\.]+"
    ) %>%
    .[[1]] %>%
    str_split(",") %>%
    map(as.numeric)

  append(points, list(points[[1]])) %>%
    data.frame() %>%
    t() %>%
    data.frame() %>%
    `colnames<-`(c("x", "y")) %>%
    `rownames<-`(seq(1, length(points) + 1))
}

mk_highlight <- function(so,
                         fov,
                         size = 0.5,
                         linewidth = 2,
                         pal = NULL) {
  ImageDimPlot(
    so,
    fov = fov,
    dark.background = FALSE,
    group.by = "seurat_clusters",
    size = size,
    cols = pal,
  ) + coord_fixed() + geom_path(
    data = get_heart_bb(SAMPLE_DIRS[[fov]]),
    aes(x = x, y = y),
    inherit.aes = FALSE,
    color = "red",
    linewidth = linewidth
  )
}

feature_plot <- function(size, fov, gene) {
  ImageFeaturePlot(
    SO_HEART,
    features = gene,
    fov = fov,
    size = size,
    dark.background = FALSE,
  ) + coord_fixed()
}


# figure 1
wrap_plots(
  ... =
    list(
      wrap_plots(
        mk_highlight(
          SO,
          "AH59.6",
          size = 1.5,
          pal = DISCRETE_PALETTE
        ) + theme(legend.position = "none"),
        ImageDimPlot(
          SO_HEART,
          fov = "AH59.6",
          dark.background = FALSE,
          group.by = "seurat_clusters",
          size = 3,
          cols = DISCRETE_PALETTE,
        ) + coord_fixed() + theme(legend.position = "none"),
        mk_highlight(
          SO,
          "AH59.2",
          size = 1.5,
          pal = DISCRETE_PALETTE
        ) + theme(legend.position = "none"),
        ImageDimPlot(
          SO_HEART,
          fov = "AH59.2",
          dark.background = FALSE,
          group.by = "seurat_clusters",
          size = 3,
          cols = DISCRETE_PALETTE,
        ) +
          coord_fixed() +
          theme(legend.position = "bottom") +
          guides(fill = guide_legend(title = "cluster")),
        design = "
            AAABBCCC
            AAABBCCC
            AAADDCCC
            AAADDCCC
          "
      ),
      DimPlot(
        SO_HEART,
        group.by = "seurat_clusters",
        cols = DISCRETE_PALETTE
      ) + labs(title = "leiden clusters"),
      DimPlot(SO_HEART, group.by = "condition"),
      DoHeatmap(
        SO_HEART,
        features = rev(c("Postn", "Sox9", "Pdgfra", "Penk")),
        label = FALSE,
        group.colors = DiscretePalette(
          DISCRETE_PALETTE,
          n = length(levels(SO_HEART))
        )
      ) +
        scale_colour_manual(
          values = DiscretePalette(
            DISCRETE_PALETTE,
            n = length(levels(SO_HEART))
          ),
          name = "cluster"
        ) +
        scale_fill_viridis_c(
          option = CONTINUOUS_PALETTE,
          name = "expression"
        ),
      map(c("nchet", "ncko"), function(co) {
        data.frame(
          count = SO_HEART[[c("seurat_clusters", "condition")]] %>%
            table() %>%
            .[c(2, 3, 7), co],
          cluster = c("2", "3", "7"),
          condition = co
        ) %>%
          mutate(pct = count / sum(count) * 100)
      }) %>%
        reduce(rbind) %>%
        ggplot() +
        aes(y = pct, x = condition, fill = cluster) +
        geom_bar(position = "stack", stat = "identity") +
        theme_classic() +
        labs(y = "percentage of cells") +
        theme(axis.title.y = element_blank()) +
        scale_y_reverse() +
        scale_fill_manual(
          values = DiscretePalette(
            DISCRETE_PALETTE,
            n = length(levels(SO_HEART))
          )[c(2, 3, 7)]
        ) +
        coord_flip(),
      wrap_plots(
        append(
          map(
            c("Postn", "Sox9", "Pdgfra", "Penk"),
            partial(feature_plot, 2.5, "AH59.6")
          ),
          map(
            c("Postn", "Sox9", "Pdgfra", "Penk"),
            partial(feature_plot, 2.5, "AH59.2")
          )
        ),
        ncol = 4,
        nrow = 2,
        guides = "collect"
      ) &
        theme(legend.position = "bottom") &
        scale_fill_viridis_c(
          option = CONTINUOUS_PALETTE,
          name = "expression",
          limits = c(0, 7.5)
        )
    ),
  design = "
      AAAAAAAA
      AAAAAAAA
      AAAAAAAA
      AAAAAAAA
      BBBDDDDD
      BBBDDDDD
      CCCDDDDD
      CCCEEEEE
      FFFFFFFF
      FFFFFFFF
      FFFFFFFF
    "
) & theme(
  text = element_text(size = FONT_SIZE)
)
ggsave("out/fig_1.pdf", w = 210 * 30, h = 297 * 30, units = "px")

wrap_plots(
  ... =
    list(
      wrap_plots(
        ... = map(c(2, 3, 7), function(cluster) {
          m <- FindMarkers(
            SO_CNCC %>% subset(seurat_clusters == cluster),
            ident.1 = "ncko",
            ident.2 = "nchet",
            group.by = "condition",
            method = "DESeq2"
          ) %>%
            subset(p_val_adj <= 0.05) %>%
            arrange(abs(avg_log2FC))

          rbind(
            data.frame(
              pct = m[["pct.1"]],
              gene = rownames(m),
              condition = "ncko",
              avg_log2FC = m[["avg_log2FC"]]
            ),
            data.frame(
              pct = m[["pct.2"]],
              gene = rownames(m),
              condition = "nchet",
              avg_log2FC = 0
            )
          ) %>%
            mutate(
              gene = do.call(
                partial(fct_relevel, gene),
                as.list(rownames(m))
              )
            ) %>%
            ggplot() +
            aes(x = condition, y = gene, colour = avg_log2FC, size = pct) +
            geom_point() +
            theme_classic() +
            theme(plot.title = element_text(face = "bold")) +
            scale_colour_viridis_c(
              limits = c(-1.5, 1.5),
              breaks = seq(-1.5, 1.5, 0.75),
              labels = seq(-1.5, 1.5, 0.75),
              option = CONTINUOUS_PALETTE,
              name = "average\nlogFC"
            ) +
            scale_radius(
              limits = c(0.1, 1),
              name = "percent\nexpression"
            ) +
            labs(
              title = sprintf("cluster %d", cluster)
            )
        }),
        guides = "collect",
        axes = "collect_x",
        design = "
            A
            A
            A
            B
            B
            C
          "
      ) + theme(legend.position = "right"),
      wrap_plots(
        ... = map(c(2, 3, 7), function(cluster) {
          so <- SO_CNCC %>% subset(seurat_clusters == cluster)
          m <- FindMarkers(
            so,
            ident.1 = "ncko",
            ident.2 = "nchet",
            group.by = "condition",
            method = "DESeq2"
          ) %>%
            subset(p_val_adj <= 0.05) %>%
            arrange(desc(abs(avg_log2FC)))

          DoHeatmap(
            so,
            rownames(m),
            group.by = "condition",
            label = FALSE
          ) + labs(
            title = sprintf("cluster %d", cluster)
          ) + scale_fill_viridis_c(
            limits = c(-2.5, 2.5),
            option = CONTINUOUS_PALETTE,
            name = "expression"
          ) + theme(
            plot.title = element_text(face = "bold")
          ) + scale_colour_discrete(name = "condition")
        }),
        guides = "collect",
        design = "
            A
            A
            A
            B
            B
            C
          "
      ) & theme(legend.position = "bottom"),
      wrap_plots(
        ... = imap(c(2, 3, 7), function(cluster, i) {
          so <- subset(SO_HEART, seurat_clusters == cluster)
          m <- FindMarkers(
            so,
            ident.1 = "ncko",
            ident.2 = "nchet",
            group.by = "condition",
            method = "DESeq2"
          ) %>%
            subset(p_val_adj <= 0.05) %>%
            arrange(desc(abs(avg_log2FC))) %>%
            rownames()

          wrap_plots(
            VlnPlot(
              so,
              features = head(m, 8),
              group.by = "condition",
              combine = FALSE,
              pt.size = FALSE
            ) %>% imap(function(p, j) {
              free(p) + labs(
                y = if (j == 1) {
                  sprintf("cluster %d", cluster)
                }
              ) + theme(
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = if (j != 1) {
                  element_blank()
                } else {
                  element_text(angle = 45, vjust = 0.5)
                }
              )
            }),
            nrow = 1
          )
        }),
        guides = "collect",
        nrow = 3
      )
    ),
  design = "
      AABBBBBB
      AABBBBBB
      AABBBBBB
      AABBBBBB
      AACCCCCC
      AACCCCCC
    "
) & theme(
  text = element_text(size = FONT_SIZE)
)
ggsave("out/fig_2.pdf", w = 210 * 30, h = 297 * 30, units = "px")

DoHeatmap(
  SO_HEART,
  features = FindAllMarkers(
    SO_HEART,
    verbose = FALSE,
    only.pos = TRUE,
    method = "DESeq2"
  ) %>%
    subset(p_val_adj <= 0.05) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 5) %>%
    .[["gene"]],
  label = TRUE,
  group.colors = DiscretePalette(
    DISCRETE_PALETTE,
    n = length(levels(SO_HEART))
  )
) +
  scale_colour_manual(
    values = DiscretePalette(
      DISCRETE_PALETTE,
      n = length(levels(SO_HEART))
    ),
    name = "cluster"
  ) +
  scale_fill_viridis_c(option = CONTINUOUS_PALETTE) +
  theme(legend.position = "bottom") +
  guides(colour = FALSE)

ggsave("out/fig_s5.pdf", w = 210 * 30 * 2, h = 297 * 30, units = "px")
