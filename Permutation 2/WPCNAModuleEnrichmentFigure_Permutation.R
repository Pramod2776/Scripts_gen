library(tidyverse)
library(ggdendro)
library(cowplot)

# Data load
gene_lists <- lapply(
    setNames(nm = readxl::excel_sheets('src/ASDRelevantGeneListsFromLiterature.xlsx')[-13]),
    function(sheet) {
        readxl::read_xlsx(
            'src/ASDRelevantGeneListsFromLiterature.xlsx',
            sheet = sheet
        )[, 1] %>%
            distinct()
    }
)
gene_lists <- gene_lists[!names(gene_lists) %in% c(
    "EichlerDNM_LGD_MIS_ASD_BD_SCZ", 
    "EichlerDNM_LGD_ASD_BD_SCZ",
    "VulnerableASD",
    "ASDSanders65"
)]

eigengene_files <- list.files(
    'src', pattern = "*eigengene*", full.names = TRUE
)
names(eigengene_files) <- lapply(eigengene_files, function(fn) {
    paste(
        str_split(basename(fn), pattern = "_", simplify = TRUE)[, 1],
        str_split(basename(fn), pattern = "_", simplify = TRUE)[, 2],
        sep = "_"
    )
})
eigengene_data <- lapply(
    eigengene_files,
    function(fn) {
        readxl::read_xlsx(fn, sheet = 1) %>%
            mutate(filename = basename(fn)) %>%
            mutate(
                Reg = str_split(filename, pattern = "_", simplify = TRUE)[, 1],
                Age = str_split(filename, pattern = "_", simplify = TRUE)[, 2]
            ) %>%
            mutate(Sample = str_split(Samplename, "_", simplify = TRUE)[, 1])
    }
)
module_membership_files <- list.files(
    'src', pattern = "*module_genes*", full.names = TRUE
)
names(module_membership_files) <- lapply(module_membership_files, function(fn) {
    paste(
        str_split(basename(fn), pattern = "_", simplify = TRUE)[, 4],
        str_split(basename(fn), pattern = "_", simplify = TRUE)[, 3],
        sep = "_"
    )
})
module_membership_data <- lapply(
    module_membership_files,
    function(fn) {
        lapply(
            setNames(nm = readxl::excel_sheets(fn)),
            function(sheet) {
                readxl::read_xlsx(fn, sheet = sheet) %>%
                    mutate(filename = basename(fn)) %>%
                    mutate(
                        Reg = str_split(filename, pattern = "_", simplify = TRUE)[, 4],
                        Age = str_split(filename, pattern = "_", simplify = TRUE)[, 3]
                    ) %>%
                    mutate(label = sheet)
            }
        ) %>%
            bind_rows()
    }
)


metadata <- lapply(
    Sys.glob("src/meta_group*\\.xlsx"),
    function(fn) {
        readxl::read_xlsx(fn, sheet = 1) %>%
            mutate(filename = fn) %>%
            mutate(
                Reg = str_split(filename, pattern = "_", simplify = TRUE)[, 4] %>%
                    str_replace(".xlsx", ""),
                Age = str_split(filename, pattern = "_", simplify = TRUE)[, 3]
            )
    }
) %>%
    bind_rows() %>%
    mutate(Genotype = sapply(SampleGroup, function(x) {
        if (x == "Cul3 HET") { return("HET") } else { return(x) }
    })) %>%
    filter(Genotype %in% c("WT", "HET")) %>%
    mutate(Genotype = relevel(factor(Genotype), ref = "WT"))

# Define Region-Age combinations
reg_age <- metadata[, c("Reg", "Age")] %>%
    arrange(Reg, Age) %>%
    distinct() %>%
    mutate(RegAge = paste(Reg, Age, sep = "_"))

all_plots <- list()
for (i in seq(1, nrow(reg_age))) {
    this_reg <- as.character(reg_age[i, "Reg"])
    this_age <- as.character(reg_age[i, "Age"])
    this_reg_age <- as.character(reg_age[i, "RegAge"])
    this_eigengene <- eigengene_data[[this_reg_age]] %>%
        inner_join(
            metadata, 
            by = c("Sample" = "SampleName", "Reg", "Age")
        )
    this_module_membership <- module_membership_data[[this_reg_age]]
    all_genes <- unique(c(toupper(this_module_membership), toupper(unlist(gene_lists))))
    # Calculate all data
    # Module-trait association
    module_trait <- lapply(
        setNames(nm = colnames(dplyr::select(this_eigengene, starts_with("ME")))),
        function(module) {
            this_formula <- paste(module, "~ Genotype")
            lm(as.formula(this_formula), data = this_eigengene) %>%
                broom::tidy() %>%
                mutate(label = module) %>%
                mutate(fdr = p.adjust(p.value, method = "BH")) %>%
                mutate(Genotype = c("WT", "HET")) %>%
                filter(Genotype == "HET")
        }
    ) %>%
        bind_rows() %>%
        mutate(Star = ifelse(fdr <= 0.05, "*", ""))
    # List enrichment
    list_enrichment <- lapply(
        setNames(nm = unique(this_module_membership$label)),
        function(ml) {
            lapply(
                names(gene_lists),
                function(gl) {
                    this_list <- toupper(pull(gene_lists[[gl]], gene_symbol))
                    module_list <- toupper(pull(filter(this_module_membership, label == ml), Gene_name))
                    actual_overlap <- length(intersect(this_list, module_list))
                    
                    null_dist_file <- list.files("./Permutation_after_WPCNA/Module_gene_overlap/", 
                                                this_reg_age, full.names = T)
                    null_dist <- read.csv(null_dist_file) %>%
                      filter(Module_number == ml) %>%
                      filter(GENELIST == gl) %>%
                      pull(Overlap)
                    
                    p_val <- (sum(null_dist > actual_overlap)+1)/(length(null_dist)+1)
                    
                    data.frame(Sample = c(this_reg_age),
                               Region = c(this_reg),
                               Period = c(this_age),
                               label = c(str_replace(ml, "M", "ME")),
                               GENELIST = c(gl),
                               P_Val = p_val
                    )
                    
                }
            ) %>%
                bind_rows() %>%
                mutate(fdr_bh = p.adjust(P_Val, method = "BH"),
                       fdr_bonferroni = p.adjust(P_Val, method = "bonferroni"))
        }
    ) %>%
        bind_rows() %>%
        mutate(Star_P_Val = ifelse(P_Val <= 0.05, "*", ""),
               Star_BH = ifelse(fdr_bh <= 0.1, "*", ""),
               Star_Bonferroni = ifelse(fdr_bonferroni <= 0.1, "*", ""))
    
    selected_modules <- unique(c(
        pull(filter(module_trait, Star == "*"), label),
        pull(filter(list_enrichment, Star_Bonferroni == "*"), label)
    ))
    selected_modules <- selected_modules[selected_modules != "ME0"]
    # Dendrogram
    cluster_MEs <- this_eigengene[, colnames(this_eigengene) %in% selected_modules] %>%
        dplyr::select(starts_with("ME")) %>%
        t() %>%
        as.matrix() %>%
        dist() %>%
        hclust()
    dendrogram_MEs <- as.dendrogram(cluster_MEs)
    dendrogram_data_MEs <- dendro_data(dendrogram_MEs)
    dendrogram_segments_MEs <- segment(dendrogram_data_MEs)
    # This data will be used to align ALL plots
    dendrogram_labels_MEs <- dendrogram_data_MEs$labels %>%
        mutate(module_label = as.numeric(str_replace(label, "ME", ""))) %>%
        mutate(module_colour = WGCNA::labels2colors(module_label))
    dendrogram_plt <- ggplot() +
        coord_equal() +
        geom_segment(
            data = dendrogram_segments_MEs,
            mapping = aes(x = x, y = y, xend = xend, yend = yend)
        ) +
        geom_point(
            data = dendrogram_labels_MEs,
            mapping = aes(x = x, y = y, fill = module_colour),
            shape = 21, show.legend = FALSE, size = 5
        ) +
        geom_text(
            data = dendrogram_labels_MEs,
            mapping = aes(x = x, y = y, label = label), 
            size = 5, nudge_y = -0.5, angle = -90, hjust = 0
        ) +
        geom_tile(
            data = dendrogram_labels_MEs,
            mapping = aes(x = x, y = y),
            alpha = 0
        ) +
        scale_fill_manual(values = setNames(nm = dendrogram_labels_MEs$module_colour)) +
        scale_y_continuous(expand = expand_scale(mult = c(0.75, 0.001))) +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank()
        )
    dendrogram_plt
    module_trait_plot_data <- left_join(
        module_trait, dendrogram_labels_MEs, by = "label"
    ) %>%
        filter(label != "ME0")
    module_trait_plt <- ggplot() +
        coord_equal() +
        geom_tile(
            data = module_trait_plot_data,
            mapping = aes(x = x, y = Genotype, fill = estimate),
            colour = "grey"
        ) +
        geom_text(
            data = module_trait_plot_data,
            mapping = aes(x = x, y = Genotype, label = Star),
            size = 10, vjust = 0.75
        ) +
        scale_fill_gradient2(
            low = "lightblue", mid = "white", high = "firebrick1", 
            midpoint = 0, limits = c(-1, 1)
        ) +
        guides(
            fill = guide_colourbar(title = "Linear Regression Beta")
        ) +
        theme(
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank()
        )
    module_trait_plt
    list_enrichment_plot_data <- left_join(
        list_enrichment, dendrogram_labels_MEs, by = "label"
    ) %>%
        filter(label != "ME0")
    list_enrichment_plt <- ggplot() +
        coord_equal() +
        geom_tile(
            data = list_enrichment_plot_data %>%
                filter(Star_Bonferroni == "*"),
            mapping = aes(x = x, y = GENELIST, fill = -log10(fdr_bonferroni)),
            colour = "grey"
        ) +
        geom_tile(
            data = list_enrichment_plot_data %>%
                filter(Star_Bonferroni == "") %>%
                filter(label %in% selected_modules),
            mapping = aes(x = x, y = GENELIST),
            colour = "grey", fill = "white"
        ) +
        # geom_tile(
        #     data = list_enrichment_plot_data,
        #     mapping = aes(x = x, y = gene_list, fill = -log10(fdr)),
        #     colour = "grey"
        # ) +
        geom_text(
            data = list_enrichment_plot_data,
            mapping = aes(x = x, y = GENELIST, label = Star_Bonferroni),
            size = 10, vjust = 0.75
        ) +
        scale_fill_gradient(high = "darkorchid", low = "white") +
        theme(
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank()
        )
    list_enrichment_plt
    # writexl::write_xlsx(
    #     x = list("ModuleTraitAssociation" = module_trait_plot_data, "ListEnrichment" = list_enrichment_plot_data),
    #     path = paste0("./Permutation_after_WPCNA/Output_files/WPCNA_TraitAssoc_ListEnr_", this_reg_age, ".xlsx")
    # )
    this_plot <- egg::ggarrange(
        dendrogram_plt,
        # module_trait_plt,
        list_enrichment_plt,
        ncol = 1
    )
    ggsave(
        filename = paste0("./Permutation_after_WPCNA/Figures/", this_reg_age, "_Purple_NegLog10FDR_WPCNA_FDR_Bonferroni_redo.pdf"),
        plot = this_plot,
        width = 12, height = 12, device = "pdf"
    )
    all_plots[[this_reg_age]] <- this_plot
}

### Region-based combine plots
for (reg in unique(reg_age$Reg)) {
    module_trait_data <- list()
    list_enrichment_data <- list()
    cluster_me_segments_data <- list()
    cluster_me_labels_data <- list()
    for (age in unique(reg_age$Age)) {
        this_reg_age <- paste(reg, age, sep = "_")
        this_eigengene <- eigengene_data[[this_reg_age]] %>%
            inner_join(
                metadata, 
                by = c("Sample" = "SampleName", "Reg", "Age")
            )
        this_eigengene <- eigengene_data[[this_reg_age]] %>%
            inner_join(
                metadata, 
                by = c("Sample" = "SampleName", "Reg", "Age")
            )
        this_module_membership <- module_membership_data[[this_reg_age]]
        all_genes <- unique(c(toupper(this_module_membership), toupper(unlist(gene_lists))))
        module_trait <- lapply(
            setNames(nm = colnames(dplyr::select(this_eigengene, starts_with("ME")))),
            function(module) {
                this_formula <- paste(module, "~ Genotype")
                lm(as.formula(this_formula), data = this_eigengene) %>%
                    broom::tidy() %>%
                    mutate(label = module) %>%
                    mutate(fdr = p.adjust(p.value, method = "BH")) %>%
                    mutate(Genotype = c("WT", "HET")) %>%
                    filter(Genotype == "HET")
            }
        ) %>%
            bind_rows() %>%
            mutate(Star = ifelse(fdr <= 0.05, "*", ""))
        list_enrichment <- lapply(
          setNames(nm = unique(this_module_membership$label)),
          function(ml) {
            lapply(
              names(gene_lists),
              function(gl) {
                this_list <- toupper(pull(gene_lists[[gl]], gene_symbol))
                module_list <- toupper(pull(filter(this_module_membership, label == ml), Gene_name))
                actual_overlap <- length(intersect(this_list, module_list))
                
                null_dist_file <- list.files("./Permutation_after_WPCNA/Module_gene_overlap/", 
                                             this_reg_age, full.names = T)
                null_dist <- read.csv(null_dist_file) %>%
                  filter(Module_number == ml) %>%
                  filter(GENELIST == gl) %>%
                  pull(Overlap)
                
                p_val <- (sum(null_dist > actual_overlap)+1)/(length(null_dist)+1)
                
                data.frame(Sample = c(this_reg_age),
                           Region = c(this_reg),
                           Period = c(this_age),
                           label = c(str_replace(ml, "M", "ME")),
                           GENELIST = c(gl),
                           P_Val = p_val
                )
                
              }
            ) %>%
              bind_rows() %>%
              mutate(fdr_bh = p.adjust(P_Val, method = "BH"),
                     fdr_bonferroni = p.adjust(P_Val, method = "bonferroni"))
          }
        ) %>%
          bind_rows() %>%
          mutate(Star_P_Val = ifelse(P_Val <= 0.05, "*", ""),
                 Star_BH = ifelse(fdr_bh <= 0.1, "*", ""),
                 Star_Bonferroni = ifelse(fdr_bonferroni <= 0.1, "*", ""))
        selected_modules <- unique(c(
            pull(filter(module_trait, Star == "*"), label),
            pull(filter(list_enrichment, Star_Bonferroni == "*"), label)
        ))
        selected_modules <- selected_modules[selected_modules != "ME0"]
        # Dendrogram
        cluster_MEs <- this_eigengene[, colnames(this_eigengene) %in% selected_modules] %>%
            dplyr::select(starts_with("ME")) %>%
            t() %>%
            as.matrix() %>%
            dist() %>%
            hclust()
        dendrogram_MEs <- as.dendrogram(cluster_MEs)
        dendrogram_data_MEs <- dendro_data(dendrogram_MEs)
        dendrogram_segments_MEs <- segment(dendrogram_data_MEs)
        # This data will be used to align ALL plots
        dendrogram_labels_MEs <- dendrogram_data_MEs$labels %>%
            mutate(module_label = as.numeric(str_replace(label, "ME", ""))) %>%
            mutate(module_colour = WGCNA::labels2colors(module_label))
        module_trait_plot_data <- left_join(
            module_trait, dendrogram_labels_MEs, by = "label"
        ) %>%
            filter(label != "ME0")
        list_enrichment_plot_data <- left_join(
            list_enrichment, dendrogram_labels_MEs, by = "label"
        ) %>%
            filter(label != "ME0")
        module_trait_data[[age]] <- module_trait_plot_data %>% mutate(Age = age)
        list_enrichment_data[[age]] <- list_enrichment_plot_data %>% mutate(Age = age)
        cluster_me_segments_data[[age]] <- dendrogram_segments_MEs %>% mutate(Age = age)
        cluster_me_labels_data[[age]] <- dendrogram_labels_MEs %>% mutate(Age = age)
    }
    module_trait_plot_data <- bind_rows(module_trait_data)
    list_enrichment_plot_data <- bind_rows(list_enrichment_data)
    cluster_me_segments_plot_data <- bind_rows(cluster_me_segments_data)
    cluster_me_labels_plot_data <- bind_rows(cluster_me_labels_data)
    
    module_trait_plt <- ggplot() +
        facet_grid(. ~ Age, scales = "free_x", space = "free_x") +
        geom_tile(
            data = module_trait_plot_data,
            mapping = aes(x = x, y = Genotype, fill = estimate),
            colour = "grey"
        ) +
        geom_text(
            data = module_trait_plot_data,
            mapping = aes(x = x, y = Genotype, label = Star),
            size = 10, vjust = 0.75
        ) +
        scale_fill_gradient2(
            low = "lightblue", mid = "white", high = "firebrick1", 
            midpoint = 0, limits = c(-1, 1)
        ) +
        guides(
            fill = guide_colourbar(title = "Linear Regression Beta")
        ) +
        theme(
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank()
        )
    module_trait_plt
    list_enrichment_plt <- ggplot() +
        facet_grid(. ~ Age, scales = "free_x", space = "free_x") +
        geom_tile(
            data = list_enrichment_plot_data,
            mapping = aes(x = x, y = GENELIST, fill = -log10(fdr_bonferroni)),
            colour = "grey"
        ) +
        geom_tile(
            data = list_enrichment_plot_data %>%
                filter(Star_Bonferroni == "") %>%
                filter(label %in% selected_modules),
            mapping = aes(x = x, y = GENELIST),
            colour = "grey", fill = "white"
        ) +
        # geom_tile(
        #     data = list_enrichment_plot_data,
        #     mapping = aes(x = x, y = gene_list, fill = -log10(fdr)),
        #     colour = "grey"
        # ) +
        geom_text(
            data = list_enrichment_plot_data,
            mapping = aes(x = x, y = GENELIST, label = Star_Bonferroni),
            size = 10, vjust = 0.75
        ) +
        scale_fill_gradient(high = "darkorchid", low = "white") +
        theme(
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank()
        )
    list_enrichment_plt
    dendrogram_plt <- ggplot() +
        facet_grid(. ~ Age, scales = "free_x", space = "free_x") +
        geom_segment(
            data = cluster_me_segments_plot_data,
            mapping = aes(x = x, y = y, xend = xend, yend = yend)
        ) +
        geom_point(
            data = cluster_me_labels_plot_data,
            mapping = aes(x = x, y = y, fill = module_colour),
            shape = 21, show.legend = FALSE, size = 5
        ) +
        geom_text(
            data = cluster_me_labels_plot_data,
            mapping = aes(x = x, y = y, label = label), 
            size = 5, nudge_y = -0.5, angle = -90, hjust = 0
        ) +
        geom_tile(
            data = cluster_me_labels_plot_data,
            mapping = aes(x = x, y = y),
            alpha = 0
        ) +
        scale_fill_manual(values = setNames(nm = cluster_me_labels_plot_data$module_colour)) +
        scale_y_continuous(expand = expand_scale(mult = c(0.75, 0.001))) +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank()
        )
    dendrogram_plt
    
    this_plot <- plot_grid(
        plot_grid(
            dendrogram_plt,
            module_trait_plt +
                theme(
                    strip.background = element_blank(),
                    strip.text = element_blank(),
                    legend.position = "none"
                ),
            list_enrichment_plt +
                theme(
                    strip.background = element_blank(),
                    strip.text = element_blank(),
                    legend.position = "none"
                ),
            ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 0.2, 0.8)
        ),
        plot_grid(
            get_legend(module_trait_plt), get_legend(list_enrichment_plt),
            ncol = 1
        ),
        rel_widths = c(1, 0.1)
    )
    this_plot
    ggsave(
        filename = paste0("./Permutation_after_WPCNA/Figures/", reg, "_Purple_NegLog10FDR_WPCNA_Combined_FDR_Bonferroni.pdf"),
        plot = this_plot,
        width = 18, height = 6, device = "pdf", useDingbats = FALSE
    )
}
