suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

get_heart_bb <- function(fov) {
  points <- paste(
    "data/mtx",
    str_replace(fov, fixed("."), "-"),
    "bb.csv",
    sep = "/"
  ) |>
    read.csv(
      row.names = 1
    ) |>
    _[["bb.heart", "coords"]] |>
    str_extract_all(
      string = _,
      "\\-?[0-9\\.]+, \\-?[0-9\\.]+"
    ) |>
    _[[1]] |>
    str_split(",") |>
    map(as.numeric)

  append(points, list(points[[1]])) |>
    data.frame() |>
    t() |>
    data.frame() |>
    `colnames<-`(c("x", "y")) |>
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
    data = get_heart_bb(fov),
    aes(x = x, y = y),
    inherit.aes = FALSE,
    color = "red",
    linewidth = linewidth
  )
}

fig.1.1 <- function(so_full, so_heart, d.pal) {
  wrap_plots(
    mk_highlight(
      so_full,
      "AH59.6",
      size = 1.5,
      pal = d.pal
    ) + theme(legend.position = "none"),
    ImageDimPlot(
      so_heart,
      fov = "AH59.6",
      dark.background = FALSE,
      group.by = "seurat_clusters",
      size = 3,
      cols = d.pal,
    ) + coord_fixed() + theme(legend.position = "none"),
    mk_highlight(
      so_full,
      "AH59.2",
      size = 1.5,
      pal = d.pal
    ) + theme(legend.position = "none"),
    ImageDimPlot(
      so_heart,
      fov = "AH59.2",
      dark.background = FALSE,
      group.by = "seurat_clusters",
      size = 3,
      cols = d.pal,
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
  )
}

fig.1.2 <- function(so_heart, d.pal) {
  DimPlot(
    so_heart,
    group.by = "seurat_clusters",
    cols = d.pal
  ) + labs(title = "leiden clusters")
}

fig.1.3 <- function(so_heart) {
  DimPlot(SO_HEART, group.by = "condition")
}

fig.1.4 <- function(so_heart, d.pal, c.pal) {
  d.pal <- DiscretePalette(
    d.pal,
    n = length(levels(so_heart))
  )
  DoHeatmap(
    so_heart,
    features = rev(c("Postn", "Sox9", "Pdgfra", "Penk")),
    label = FALSE,
    raster = FALSE,
    group.colors = d.pal
  ) +
    scale_colour_manual(
      values = d.pal,
      name = "cluster"
    ) +
    scale_fill_viridis_c(
      option = c.pal,
      name = "expression"
    )
}

fig.1.5 <- function(so_heart, d.pal) {
  d.pal <- DiscretePalette(
    d.pal,
    n = length(levels(so_heart))
  )[c(2, 3, 7)]
  map(c("nchet", "ncko"), function(co) {
    data.frame(
      count = so_heart[[c("seurat_clusters", "condition")]] |>
        table() |>
        _[c(2, 3, 7), co],
      cluster = c("2", "3", "7"),
      condition = co
    ) |>
      mutate(pct = count / sum(count) * 100)
  }) |>
    reduce(rbind) |>
    ggplot() +
    aes(y = pct, x = condition, fill = cluster) +
    geom_bar(position = "stack", stat = "identity") +
    theme_classic() +
    labs(y = "percentage of cells") +
    theme(axis.title.y = element_blank()) +
    scale_y_reverse() +
    scale_fill_manual(
      values = d.pal
    ) +
    coord_flip()
}

fig.1.6 <- function(so_heart, c.pal) {
  feature_plot <- function(fov, gene) {
    ImageFeaturePlot(
      so_heart,
      features = gene,
      fov = fov,
      size = 2.5,
      dark.background = FALSE,
    ) + coord_fixed()
  }

  wrap_plots(
    append(
      map(
        c("Postn", "Sox9", "Pdgfra", "Penk"),
        partial(feature_plot, "AH59.6")
      ),
      map(
        c("Postn", "Sox9", "Pdgfra", "Penk"),
        partial(feature_plot, "AH59.2")
      )
    ),
    ncol = 4,
    nrow = 2,
    guides = "collect"
  ) &
    theme(legend.position = "bottom") &
    scale_fill_viridis_c(
      option = c.pal,
      name = "expression",
      limits = c(0, 7.5)
    )
}

fig.2.1 <- function(marker.df.l, c.pal) {
  wrap_plots(
    ... = imap(marker.df.l, function(m, cluster) {
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
      ) |>
        mutate(
          gene = do.call(
            partial(fct_relevel, gene),
            as.list(rownames(m))
          )
        ) |>
        ggplot() +
        aes(x = condition, y = fct_rev(gene), colour = avg_log2FC, size = pct) +
        geom_point() +
        theme_classic() +
        theme(plot.title = element_text(face = "bold")) +
        scale_colour_viridis_c(
          limits = c(-1.5, 1.5),
          breaks = seq(-1.5, 1.5, 0.75),
          labels = seq(-1.5, 1.5, 0.75),
          option = c.pal,
          name = "average\nlogFC"
        ) +
        scale_radius(
          limits = c(0.1, 1),
          name = "percent\nexpression"
        ) +
        labs(
          title = sprintf("cluster %s", cluster),
          y = "gene"
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
  ) + theme(legend.position = "right")
}

fig.2.2 <- function(so_cncc, marker.df.l, c.pal) {
  wrap_plots(
    ... = imap(marker.df.l, function(m, cluster) {
      so <- subset(so_cncc, seurat_clusters == cluster)

      DoHeatmap(
        so,
        rownames(m),
        group.by = "condition",
        raster = FALSE,
        label = FALSE
      ) + labs(
        title = sprintf("cluster %s", cluster)
      ) + scale_fill_viridis_c(
        limits = c(-2.5, 2.5),
        option = c.pal,
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
  ) & theme(legend.position = "bottom")
}

fig.2.3 <- function(so_cncc, marker.df.l, c.pal) {
  wrap_plots(
    ... = imap(marker.df.l, function(m, cluster) {
      so <- subset(so_cncc, seurat_clusters == cluster)

      wrap_plots(
        VlnPlot(
          so,
          features = rownames(m) |>
            head(8),
          group.by = "condition",
          combine = FALSE,
          pt.size = FALSE
        ) |> imap(function(p, j) {
          free(p) + labs(
            y = if (j == 1) {
              sprintf("cluster %s", cluster)
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
}

fig.s <- function(so_heart, d.pal, c.pal) {
  d.pal <- DiscretePalette(
    DISCRETE_PALETTE,
    n = length(levels(so_heart))
  )
  DoHeatmap(
    so_heart,
    features = FindAllMarkers(
      so_heart,
      verbose = FALSE,
      only.pos = TRUE
    ) |>
      subset(p_val_adj < 0.05) |>
      group_by(cluster) |>
      slice_max(avg_log2FC, n = 5) |>
      _[["gene"]],
    label = TRUE,
    raster = FALSE,
    group.colors = d.pal
  ) +
    scale_colour_manual(
      values = d.pal,
      name = "cluster"
    ) +
    scale_fill_viridis_c(option = c.pal) +
    theme(legend.position = "bottom") +
    guides(colour = FALSE)
}

if (TRUE) {
  DISCRETE_PALETTE <- "glasbey"
  CONTINUOUS_PALETTE <- "plasma"
  FONT_SIZE <- 20
  OUT_FORMAT <- "pdf"

  SO <- readRDS("out/all_cells.rds")
  SO_HEART <- readRDS("out/heart_only.rds")
  SO_CNCC <- subset(SO_HEART, seurat_clusters %in% c(2, 3, 7))

  CNCC_CLUSTER_MARKERS <- map(
    setNames(nm = c(2, 3, 7)),
    function(cluster) {
      FindMarkers(
        subset(SO_CNCC, seurat_clusters == cluster),
        ident.1 = "ncko",
        ident.2 = "nchet",
        group.by = "condition"
      ) |>
        subset(p_val_adj < 0.05) |>
        arrange(desc(abs(avg_log2FC)))
    }
  )

  wrap_plots(
    fig.1.1(SO, SO_HEART, DISCRETE_PALETTE),
    fig.1.2(SO_HEART, DISCRETE_PALETTE),
    fig.1.3(SO_HEART),
    fig.1.4(SO_HEART, DISCRETE_PALETTE, CONTINUOUS_PALETTE),
    fig.1.5(SO_HEART, DISCRETE_PALETTE),
    fig.1.6(SO_HEART, CONTINUOUS_PALETTE),
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
  ggsave(sprintf("out/fig_1.%s", OUT_FORMAT), w = 210 * 30, h = 297 * 30, units = "px")

  wrap_plots(
    fig.2.1(CNCC_CLUSTER_MARKERS, CONTINUOUS_PALETTE),
    fig.2.2(SO_CNCC, CNCC_CLUSTER_MARKERS, CONTINUOUS_PALETTE),
    fig.2.3(SO_CNCC, CNCC_CLUSTER_MARKERS, CONTINUOUS_PALETTE),
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
  ggsave(sprintf("out/fig_2.%s", OUT_FORMAT), w = 210 * 30, h = 297 * 30, units = "px")

  fig.s(SO_HEART, DISCRETE_PALETTE, CONTINUOUS_PALETTE)
  ggsave(sprintf("out/fig_s.eps", OUT_FORMAT), w = 210 * 30 * 2, h = 297 * 30, units = "px")
}
