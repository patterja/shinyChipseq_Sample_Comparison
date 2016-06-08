shinyServer(function(input, output) {

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Depth of Coverage 

output$plt_gibb1_macs <- renderPlot({
  ggplot(k_macsjoin, aes(x=coverage)) +
   geom_histogram(binwidth=1) +
   geom_rect(data=data.frame(xmin=min(input$depth1_macs), xmax=max(input$depth1_macs)), 
             aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf),
             color="red", alpha=0.5, inherit.aes = FALSE) +
   ggtitle("Gibbon1 Coverage") +
   theme_bw() +
   theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),
         axis.title.x = element_text(size=16), axis.title.y=element_text(size=16)) +
   labs(x="depth", y="frequency") +
   ylim(0, 2000)
         
  })

output$plt_gibb2_macs <- renderPlot({
  ggplot(m_macsjoin, aes(x=coverage)) +
     geom_histogram(binwidth=1) +
     geom_rect(data=data.frame(xmin=min(input$depth2_macs), xmax=max(input$depth2_macs)), 
               aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf),
               color="red", alpha=0.5, inherit.aes = FALSE) +
     ggtitle("Gibbon2 Coverage") +
     theme_bw() +
     theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),
           axis.title.x = element_text(size=16), axis.title.y=element_text(size=16)) +
     labs(x="depth", y="frequency") +
     ylim(0, 2000)
})


output$plt_gibb1_bayp <- renderPlot({
  ggplot(k_baypjoin, aes(x=coverage)) +
     geom_histogram(binwidth=1) +
     geom_rect(data=data.frame(xmin=min(input$depth1_bayp), xmax=max(input$depth1_bayp)), 
               aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf),
               color="red", alpha=0.5, inherit.aes = FALSE) +
     theme_bw() +
     theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),
           axis.title.x = element_text(size=16), axis.title.y=element_text(size=16)) +
     labs(x="depth", y="frequency") +
     ylim(0, 2000)
})

output$plt_gibb2_bayp <- renderPlot({
  ggplot(m_baypjoin, aes(x=coverage)) +
     geom_histogram(binwidth=1) +
     geom_rect(data=data.frame(xmin=min(input$depth2_bayp), xmax=max(input$depth2_bayp)), 
               aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf),
               color="red", alpha=0.5, inherit.aes = FALSE) +
     theme_bw() +
     theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),
           axis.title.x = element_text(size=16), axis.title.y=element_text(size=16)) +
     labs(x="depth", y="frequency") +
     ylim(0, 2000)
  })
#-------------------------------------------------------------------------------------
gibbon1dat_macs<- reactive({
  k_macsjoin %>%
    filter(coverage >=min(input$depth1_macs) & coverage<=max(input$depth1_macs)) %>%
    select(chr, start, end, foldchange,neglog10pval,neglog10qval,coverage,prop_cov)
  
})

gibbon2dat_macs<- reactive({
  m_macsjoin %>% 
           filter(coverage >=min(input$depth2_macs) & coverage<=max(input$depth2_macs)) %>%
           select(chr, start, end, foldchange,neglog10pval,neglog10qval,coverage,prop_cov)
})



gibbon1dat_bayp<- reactive({
  k_baypjoin %>% 
           filter(coverage >=min(input$depth1_bayp) & coverage<=max(input$depth1_bayp)) %>%
           select("chr"=space, start, end, PP, coverage)
})

gibbon2dat_bayp<- reactive({
  m_baypjoin %>%
    filter(coverage >=min(input$depth2_bayp) & coverage<=max(input$depth2_bayp)) %>%
    select("chr"=space, start, end, PP, coverage)
})

#-------------------------------------------------------------------------------------
ov_filtmacs <-reactive({
  olRanges(query=makeGRangesFromDataFrame(gibbon1dat_macs()), subject=makeGRangesFromDataFrame(gibbon2dat_macs()), output="df") 
})

ov_filtbayp <-reactive({
  olRanges(query=makeGRangesFromDataFrame(gibbon1dat_bayp()), 
           subject=makeGRangesFromDataFrame(gibbon2dat_bayp()), output="df")
})
#-------------------------------------------------------------------------------------Venn Diagrams
output$macs_venn<-renderD3vennR({
  km<-gibbon1dat_macs()
  mm<-gibbon2dat_macs()
  ov<- ov_macs
  venn_tooltip2(
    d3vennR(
      data = list(
        list( sets = list(0), 
              label=paste("Gibbon1","MACS2", nrow(km), sep="\n"), 
              size=nrow(km), gene=""),
        list( sets = list(1), 
              label=paste("Gibbon2","MACS2", nrow(mm), sep=",\n"), 
              size=nrow(mm), gene=""), 
        list( sets = list(0,1), 
              size=length(unique(ov$Qindex)), 
              gene="")
      ),layoutFunction = '
function(d) { return venn.venn(d, { initialLayout: venn.classicMDSLayout });}
      '
      
    ))
})

output$bayp_venn<-renderD3vennR({
  km<-gibbon1dat_bayp()
  mm<-gibbon2dat_bayp()
  ov<- ov_bayp
  venn_tooltip2(
    d3vennR(
      data = list(
        list( sets = list(0), 
              label=paste("Gibbon1","BayesPeak", nrow(km), sep="\n"), 
              size=nrow(km), gene=""),
        list( sets = list(1), 
              label=paste("Gibbon2","BayesPeak", nrow(mm), sep=",\n"), 
              size=nrow(mm), gene=""), 
        list( sets = list(0,1), 
              size=length(unique(ov$Qindex)), 
              gene="")
      ),layoutFunction = '
function(d) { return venn.venn(d, { initialLayout: venn.classicMDSLayout });}
      '
      
    ))
})
#-------------------------------------------------------------------------------------STAT TABLE
output$macsstats <-renderTable({
  km<-gibbon1dat_macs()
  mm<-gibbon2dat_macs()
  ov<- ov_macs
  data.frame(        
    Metric = c("Number of Peaks",
               "Average Width",
               "Min Width",
               "Max Width",
               "Number of overlapping peaks",
               "Average % overlap of Gibbons"),
    
    Gibbon1 = c(paste0(nrow(km)," / ",nrow(k_macsjoin)),
                paste0(mean(km$end-km$start+1)," / ",mean(k_macsjoin$end-k_macsjoin$start+1)),
                paste0(min(km$end-km$start+1)),
                paste0(max(km$end-km$start+1)),
                paste0(length(unique(ov_filtmacs()$Qindex))," / 79664"),
                paste0(mean(ov$OLpercQ))),
                
                
               
    Gibbon2 = c(paste0(nrow(mm)," / ",nrow(m_macsjoin)),
                paste0(mean(mm$end-mm$start+1)," / ",mean(m_macsjoin$end-m_macsjoin$start+1)),
                paste0(min(mm$end-mm$start+1)),
                paste0(max(mm$end-mm$start+1)),
                paste0(length(unique(ov_filtmacs()$Sindex))," / 79664"),
                paste0(mean(ov$OLpercS)))
    
  )}, digits=3)

output$baypstats <-renderTable({
  km<-gibbon1dat_bayp()
  mm<-gibbon2dat_bayp()
  ov<- ov_bayp
  data.frame(        
    Metric = c("Number of Peaks",
               "Average Width",
               "Min Width",
               "Max Width",
               "Number of overlapping peaks Gibbon1&2",
               "Average % overlap of Gibbons"),
    
    Gibbon1 = c(paste0(nrow(km)," / ",nrow(k_baypjoin)),
                paste0(mean(km$end-km$start+1)," / ",mean(k_baypjoin$end-k_baypjoin$start+1)),
                paste0(min(km$end-km$start+1)),
                paste0(max(km$end-km$start+1)),
                paste0(length(unique(ov_filtbayp()$Qindex))," / 50114"),
                paste0(mean(ov$OLpercQ))),
    
    
    
    Gibbon2 = c(paste0(nrow(mm)," / ",nrow(m_baypjoin)),
                paste0(mean(mm$end-mm$start+1)," / ",mean(m_baypjoin$end-m_baypjoin$start+1)),
                paste0(min(mm$end-mm$start+1)),
                paste0(max(mm$end-mm$start+1)),
                paste0(length(unique(ov_filtbayp()$Sindex))," / 50114"),
                paste0(mean(ov$OLpercS)))
    
  )}, digits=3)
#-------------------------------------------------------------------------------------
output$gibbon1_macs<-renderDataTable({
  as.data.frame(gibbon1dat_macs())
})

output$gibbon2_macs<-renderDataTable({
  as.data.frame(gibbon2dat_macs())
})
output$gibbon1_bayp<-renderDataTable({
  as.data.frame(gibbon1dat_bayp())
})

output$gibbon2_bayp<-renderDataTable({
  as.data.frame(gibbon2dat_bayp())
})
output$ov_macstable<- renderDataTable({
  data.frame(cbind(gibbon1=k_macsjoin[ov_filtmacs()$Qindex,], 
                   gibbon2=m_macsjoin[ov_filtmacs()$Sindex,], ov_filtmacs())) %>%
    select("Overlap_Type"=OLtype, "overlap%gibbon1"=OLpercQ, "overlap%gibbon2"=OLpercS, 
           "chr"=space,gibbon1.start, gibbon1.end,
           gibbon2.start, gibbon2.end,
           gibbon1.coverage, gibbon1.width,
           gibbon2.coverage, gibbon2.width,
           gibbon1.neglog10pval, gibbon1.neglog10qval,
           gibbon2.neglog10pval, gibbon2.neglog10qval
    )
    
  
  
})

output$ov_bayptable <- renderDataTable({
  data.frame(cbind(gibbon1=k_baypjoin[ov_filtbayp()$Qindex,], 
                   gibbon2=m_baypjoin[ov_filtbayp()$Sindex,], ov_filtbayp())) %>%
    select("Overlap_Type"=OLtype, "overlap%gibbon1"=OLpercQ, "overlap%gibbon2"=OLpercS, 
           "chr"=space,gibbon1.start, gibbon1.end,gibbon1.width.x,
           gibbon2.start, gibbon2.end, gibbon2.width.x,
           gibbon1.PP, gibbon1.coverage,
           gibbon2.PP, gibbon2.coverage)
  

 })
})
