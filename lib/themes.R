## barplots ####
ggtheme_bar <- theme_bw() +
    theme(
        legend.text = element_text(family = "Times", size = 10),
        legend.title = element_text(family = "Times", size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(hjust = 1, angle = 45),
        strip.text = element_text(size = 10),
        plot.margin = unit(c(0.025,0.025,.025,0.025), "npc")
    )

ggtheme_reads <- theme_bw() +
    theme(
        legend.position = "bottom",
        legend.text = element_text(family = "Times", size = 10),
        legend.title = element_text(family = "Times", size = 12),
        axis.text = element_text(family = "Times", size = 10),
        axis.title = element_text(family = "Times", size = 12),
        strip.text = element_text(family = "Times", size = 10),
        legend.key.height = unit(1.4, "lines"),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
    )

ggtheme_alpha <- theme_bw() + theme(
    
    legend.position="top",
    legend.text = element_text(family = "Times",size = 12),
    legend.title = element_text(family = "Times",size = 12),
    axis.text = element_text(family = "Times", size = 12),
    axis.title = element_text(),
    strip.text = element_text(family = "Times", size = 12),
    plot.margin = unit(c(0.025,0.025,.025,0.025), "npc"),
    axis.text.x = element_text(family = "Times", size = rel(1), 
                               angle = 30, hjust = 1, vjust = 1),
    axis.title = element_text(family = "Times",  size = rel(1), 
                              lineheight = 1.5),
    legend.key = element_rect(colour = "white")
)

ggtheme_beta <- theme_bw() + theme(
    legend.text = element_text(family = "Times", size = 12),
    legend.title = element_text(family = "Times", size = 12),
    axis.text = element_text(family = "Times", size = 10),
    axis.title = element_text(family = "Times", size = 12),
    strip.text = element_text(family = "Times", size = 12),
    plot.margin = unit(c(0.025,0.025,.025,0.025), "npc")
    )

ggtheme_core <- theme_bw() + theme(
    legend.text = element_text(family="Times",size = 10),
    legend.title = element_text(family="Times",size = 12),
    axis.title = element_text(family="Times", size = 10),
    strip.text = element_text(family="Times", size = 10),
    plot.margin = unit(c(0.000005,0.000005,.000005,0.000005), "npc"),
    axis.text.x = element_text(family="Times", size=10, hjust = 1, vjust = 1),
    axis.ticks = element_blank(),
    axis.text.y = element_text(family="Times", size = 10),
    axis.title = element_text(family="Times",  size=10, lineheight=1.5),
    legend.key = element_rect(colour = "white"),
    #panel.grid.major = element_blank()
    panel.grid.major = element_line(size=0.025, colour="grey50")
)
