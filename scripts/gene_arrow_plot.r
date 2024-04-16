library(gggenes)
library(tidyverse)
library(patchwork)

shrink_segments <- function(starts, ends, shrink_factor) {
  # Calculate midpoints
  midpoints <- (starts + ends) / 2
  
  # Adjust start and end points
  new_starts <- midpoints + (starts - midpoints) * shrink_factor
  new_ends <- midpoints + (ends - midpoints) * shrink_factor
  
  list(starts = new_starts, ends = new_ends)
}

# Load Carlos' file
dir_bias <- read_csv('/g/scb/zeller/karcher/carlos_gene_plot/data/directionality_per_locus_same_over_total_2024-04-10.csv') %>%
    mutate(forward = strand == "+")  %>%
    mutate(feature_type =
        case_when(
            str_detect(locus_tag, "IG-") ~ "IG",
            # Carlos, here I assume IG is intergenic
            # and PM is the gene itself, but I wasn't quite sure what PM stands for
            str_detect(locus_tag, "PM") ~ "gene",
        )
    ) %>%
    # intergenic regions don't have a strandedness, this fixes that
    ## Update: I'm leaving this as you have the +/- strand info also for intergenic regions
    ## How do you define that? Based on the previous gene? I think we should discuss
    # mutate(forward = ifelse(feature_type == "IG", NA, forward)) %>%
    # Here I make it so that each gene and each IG region is exactly the same length
    arrange(species, media, feature_number) %>%
    group_by(species, media) %>%
    mutate(dummy_x_start = seq(1, 100000, length.out = n())) %>%
    mutate(dummy_x_end = seq(dummy_x_start[2], 100000+dummy_x_start[1], length.out = n())) %>%
    # Carlos, here I swich the dummy_x_end and dummy_x_start values for genes that have FALSE for the forward column
    mutate(dummy_x_start_SAVE = dummy_x_start) %>%
    mutate(dummy_x_end_SAVE = dummy_x_end) %>%
    mutate(dummy_x_start = ifelse(forward, dummy_x_start_SAVE, dummy_x_end_SAVE)) %>%
    mutate(dummy_x_end = ifelse(forward, dummy_x_end_SAVE, dummy_x_start_SAVE)) %>%
    select(-dummy_x_end_SAVE, dummy_x_start_SAVE) %>%
    identity()

# I want the arrows for the insertions to not touch each other...
shrunken_segments <- shrink_segments(dir_bias$dummy_x_start, dir_bias$dummy_x_end, 0.6)
dir_bias$dummy_x_start_aux <- shrunken_segments$start
dir_bias$dummy_x_end_aux <- shrunken_segments$end

shrunken_segments <- shrink_segments(dir_bias$dummy_x_start, dir_bias$dummy_x_end, 0.9)
dir_bias$dummy_x_start <- shrunken_segments$start
dir_bias$dummy_x_end <- shrunken_segments$end

# Just for Nic to see, can be deleted
# for (s in unique(dir_bias$species)) {
#     for (m in unique(dir_bias$media)) {
#         print(head(dir_bias %>% filter(species == s, media == m) %>% arrange(feature_number)))
#     }
# }

head(dir_bias)

get_gene_plot <- function(s, m, xstart, xend, arrow_length = 0.05, absolute_cutoff = 0.9) {

    get_aux_arrow_plot <- function(how = NULL) {
    

        arrow_data <- example_genes %>%
            filter(!is.na(ratio_dir))
        aux <- ggplot()

        if (how == "in_dir") {
            arrow_data <- arrow_data
            aux <- aux +
                geom_segment(data = arrow_data, 
                aes(
                    x = dummy_x_start_aux, 
                    xend = dummy_x_end_aux, 
                    y = 1, 
                    yend = 1, 
                    linewidth = ratio_dir,
                    alpha = ratio_dir), arrow = arrow(length = unit(arrow_length, "inches")), show.legend = TRUE)
        } else if (how == "reverse_dir") {
            arrow_data <- arrow_data %>% 
                mutate(ratio_dir = 1-ratio_dir)
            aux <- aux +
                geom_segment(data = arrow_data, 
                aes(
                    x = dummy_x_end_aux, 
                    xend = dummy_x_start_aux, 
                    y = 1, 
                    yend = 1, 
                    linewidth = ratio_dir,
                    alpha = ratio_dir), arrow = arrow(length = unit(arrow_length, "inches")), show.legend = FALSE)                
        } else {
            stop("how must be either 'in_dir' or 'reverse_dir'")
        }
        aux +
        theme_classic() +
        theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.line.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.line.y = element_blank(),
                axis.ticks.length.x = unit(0, "pt"),
                plot.margin = unit(c(0,0,0,0), "cm")) +
        scale_linewidth_continuous(range = c(0, 1.25)) +
        scale_alpha_continuous(range = c(0, 1)) +
        scale_x_continuous(limits = c(xmin, xmax)) +
        NULL
    }

    example_genes <- dir_bias %>% filter(species == s, media == m)
    example_genes <- example_genes %>% filter(dummy_x_start >= xstart, dummy_x_start <= xend)

    xmin <- min(example_genes$dummy_x_start)
    xmax <- max(example_genes$dummy_x_start)

    aux_arrow_plot_top <- get_aux_arrow_plot(how = 'in_dir')
    aux_arrow_plot_bottom <- get_aux_arrow_plot(how = 'reverse_dir')
                

    gene_plot <- ggplot() +
        # Plot the genes first
        geom_gene_arrow(data = example_genes %>% filter(feature_type == "gene"), 
        aes(xmin = dummy_x_start, xmax = dummy_x_end, y = 1, fill = ratio_dir), height = 0.5) + 
        # And now the IG regions
        geom_gene_arrow(data = example_genes %>% filter(feature_type == "IG"), 
        aes(xmin = dummy_x_start, xmax = dummy_x_end, y = 1, fill = ratio_dir), height = 0.5, arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm")) +         
        theme_minimal() +
        scale_y_continuous(expand = c(0, 0)) +
        #labs(title = str_c(s, ", ", m)) +
        theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.line.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                plot.background = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                plot.margin = unit(c(0,0,0,0), "cm")) +        
        
        scale_x_continuous(limits = c(xmin, xmax)) +
        NULL
    return(
        aux_arrow_plot_top + 
        plot_spacer() +
        aux_arrow_plot_bottom + 
        plot_spacer() +
        gene_plot + 
        
        plot_layout(ncol = 1, heights = c(0.5, 0, 0.5, 0, 2), guides = 'collect') +  plot_annotation(theme = theme(plot.margin = margin())))
}

gene_plot <- get_gene_plot(
    s = "merdae", 
    m = "liquid", 
    xstart = 1, 
    xend = 500)

ggsave(plot = gene_plot, filename = "../results/gene_plot.pdf", width = 10, height = 0.8)
