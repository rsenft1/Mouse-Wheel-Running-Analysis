pretty_boxplot <- function(df, x_var, y_var, fill_color, title, color_palette=c("#00BBF9","#FDCA40"), x_lab=x_var, y_lab=y_var, color_lab=fill_color, axis_scale="None",save=TRUE, savename=".png"){
  g <- ggplot(phase_meta, aes(y=df[[y_var]], x=df[[x_var]], fill=df[[fill_color]]))+
    stat_boxplot(geom = "errorbar", width = 0.5, position=position_dodge(0.8))+
    geom_boxplot(position=position_dodge(0.8), outlier.shape = NA)+
    geom_beeswarm(size=2.5,cex=2.5, alpha=0.8, dodge.width=.8, aes(color=.data[[fill_color]]))+
    scale_fill_manual(values = color_palette) + 
    scale_color_manual(values = darken(color_palette, amount = 0.6))
    if(axis_scale=="percent"){
      g <- g + scale_y_continuous(labels = scales::percent)
    }
    g <- g + labs(
      x = x_lab,
      y = y_lab,
      color = color_lab,
      fill = color_lab,
      title = title
    )+
    theme_classic(base_size = 18)+
    theme(axis.title.x = element_text(margin=margin(t=15)), axis.title.y = element_text(margin=margin(r=15)))
  print(g)
  if(save){
    ggsave(here(dir_figs, savename))
  }
}

pretty_facet_graph <- function(data, x_var, y_var, grid_row_var=".", grid_col_var=".", color_palette=c("#00BBF9","#FDCA40"), save=TRUE){
  g3 <- ggplot(data, aes(x=.data[[x_var]], y=.data[[y_var]], fill=.data[[x_var]])) + 
    geom_boxplot(outlier.colour = NA) +
    geom_point(position=position_dodge(width=0.75), aes(color=.data[[x_var]]), alpha=0.5)+
    facet_grid(reformulate(grid_row_var,grid_col_var))+
    scale_fill_manual(values = color_palette) + 
    scale_color_manual(values = darken(color_palette, amount = 0.5))+
    scale_y_continuous(name= "Percent time running",labels = scales::percent)+
    ggtitle(paste(phase))+
    theme_bw(base_size = 14)+
    xlab("Light cycle phase")
  print(g3)
  if(save){
    ggsave(here(dir_figs, paste0(phaseforfilename, "_percent_time_running_bySexGeno_active.png")))
  }
}
