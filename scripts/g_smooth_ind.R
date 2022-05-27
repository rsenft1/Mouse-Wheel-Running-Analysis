smoothIndActivity <- function(mouseNum, micelist, genotypes, cycleMapped, downSample=5, subset=FALSE, w=11, h=4){
  #name for mouse and genotype
  myName <- micelist[mouseNum]
  myData <- get(paste0(myName,"dt"))
  if (subset!=FALSE){
    myData <- myData[myData$timefactor%in%subset, ]
    myData[with(myData, timefactor %in% subset), ] # can have multiple conditions after %in%
  }
  myGeno <- genotypes[mouseNum]
  myGenofileName <- str_replace_all(myGeno,"/","-")
  grob <- grobTree(textGrob(myName, x=0.05,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13, fontface="bold")))
  grobGeno <- grobTree(textGrob(myGeno, x=0.05,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=13, fontface="italic")))
  # the plot itself:
  p <- ggplot(data=cycleMapped, aes(xmin=start, xmax=end, ymin=0, ymax=Inf))+
    geom_rect(fill='gray', alpha=0.6)+
    geom_area(data=myData[seq(1, nrow(myData), downSample), ], inherit.aes = FALSE, aes(x=time, y=smoothSpeed), fill="cornflowerblue", alpha=0.5, color="black")+
    #geom_line(data=myData, inherit.aes = FALSE, aes(x=time, y=smoothSpeed, group=1),size=0.5)+
    scale_y_continuous(limits = c(0, 25), expand = c(0, 0))+
    scale_x_continuous(limits=c(min(myData$time), max(myData$time)), expand = c(0, 0), breaks=seq(1,max(as.numeric(myData$timefactor)),1))+
    xlab("Days")+
    ylab("Mean wheel speed (meter/min)")+
    annotation_custom(grob)+
    annotation_custom(grobGeno)+
    theme_classic()
  ggsave(
    filename = paste0(myName,"_",myGenofileName,".pdf"),
    plot = p,
    device = "pdf",
    path = here("figs"),
    scale = 1,
    width = w,
    height = h,
    units = "in",
    dpi = 300,
    limitsize = TRUE
  )
  #ggsave(path = here("figs"), filename = "fig1.png")
  #pdf(paste0(name,".pdf"))
  #print(p)
  #dev.off()
}