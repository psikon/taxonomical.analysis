# plot_overview_bar <- function(phyloseq, level, seperator = NULL, group_by = NULL) {
#     # remove "k__" schema from legend
#     phyloseq <- remove_Underscore(phyloseq)
#     # create the string for plotting function
#     ifelse(is.null(seperator), 
#         ifelse(is.null(group_by),
#             plot_bar <- paste0("plot_bar(phyloseq, fill = level)"),
#             plot_bar <- paste0("plot_bar(phyloseq, x = group_by,fill = level)")),
#         ifelse(is.null(group_by),
#             plot_bar <- paste0("plot_bar(phyloseq, fill = level, facet_grid=~", 
#                                 as.name(seperator), ")"),
#             plot_bar <- paste0("plot_bar(phyloseq, x = group_by, fill = level, facet_grid=~", 
#                                 as.name(seperator), ")")))
#     # evaluate the barplot
#     p = eval(parse(text=plot_bar))
#     # create string for modifing the barplot
#     geom = paste0("geom_bar(aes(color = ", as.name(level),
#                 ", fill = ", as.name(level), "), stat = 'identity', position = 'stack')")
#     # draw and modify the barplot
#     bar= p + eval(parse(text=geom))
#     bar
# }
# sd <- sample_data(phylum)
# class(sd@row.names) <- as.factor(sd@row.names)
# sd
# plot_bar(phylum,facet_grid=~Environment)
# p$labels
