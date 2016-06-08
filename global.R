

library(idr)
library(GenomicRanges)
source("rangeoverlapper.R")
library(dplyr)
library(ggplot2)
library(d3vennR)


#-----------------------------------------------------Load MACS2 and BayesPeak data



karenina_macs <- read.delim("data/karenina_SSY_Siamang_CTCF_VS_inputctrl_peaks.narrowPeak", header= FALSE, sep="\t", stringsAsFactors = FALSE)
monty_macs <- read.delim("data/monty_SSY_Siamang_CTCF_VS_inputctrl_peaks.narrowPeak", header= FALSE, sep="\t", stringsAsFactors = FALSE)
karenina_bayp <- read.delim("data/karenina_thresh0.txt", header= TRUE, sep=" ", stringsAsFactors = FALSE)
monty_bayp <- read.delim("data/monty_thresh0.txt", header= TRUE, sep=" ", stringsAsFactors = FALSE)
#-----------------------------------------------------
cols <-c("chr","start","end","peaknum","visnum", ".","foldchange","neglog10pval","neglog10qval","relsummitpos_tostart")
colnames(karenina_macs) <- cols
colnames(monty_macs)<-cols

#-----------------------------------------------------Load DOC files
# doc_kmacs <- read.delim("C:/Users/Owner/Box Sync/chip_seq/coverage/karenina_coverage_macs.txt", header= FALSE, sep="\t", stringsAsFactors = FALSE)
# doc_mmacs <- read.delim("C:/Users/Owner/Box Sync/chip_seq/coverage/monty_coverage_macs.txt", header= FALSE, sep="\t", stringsAsFactors = FALSE)
# doc_kbayp <- read.delim("C:/Users/Owner/Box Sync/chip_seq/coverage/karenina_coverage_bayp-copy-1.txt", header= FALSE, sep="\t", stringsAsFactors = FALSE)
# doc_mbayp <- read.delim("C:/Users/Owner/Box Sync/chip_seq/coverage/monty_coverage_bayp-copy-1.txt", header= FALSE, sep="\t", stringsAsFactors = FALSE)

doc_kmacs <- read.delim("data/karenina_coverage_macs.txt", header= FALSE, sep="\t", stringsAsFactors = FALSE)
doc_mmacs <- read.delim("data/monty_coverage_macs.txt", header= FALSE, sep="\t", stringsAsFactors = FALSE)
doc_kbayp <- read.delim("data/karenina_coverage_bayp.txt", header= FALSE, sep="\t", stringsAsFactors = FALSE)
doc_mbayp <- read.delim("data/monty_coverage_bayp.txt", header= FALSE, sep="\t", stringsAsFactors = FALSE)
colnames(doc_kmacs)<-c("chr","start","end","coverage","bases_cov","width","prop_cov")
colnames(doc_mmacs)<-c("chr","start","end","coverage","bases_cov","width","prop_cov")
colnames(doc_kbayp)<-c("chr","start","end","coverage","bases_cov","width","prop_cov")
colnames(doc_mbayp)<-c("chr","start","end","coverage","bases_cov","width","prop_cov")
#-----------------------------------------------------
k_macsjoin<-inner_join(karenina_macs, doc_kmacs, by=c("chr"="chr", "start"="start", "end"="end"))
m_macsjoin<-inner_join(monty_macs, doc_mmacs, by=c("chr"="chr", "start"="start", "end"="end"))
k_baypjoin<-inner_join(karenina_bayp, doc_kbayp, by=c("space"="chr", "start"="start", "end"="end"))
m_baypjoin<-inner_join(monty_bayp, doc_mbayp, by=c("space"="chr", "start"="start", "end"="end"))

#----------------------------------------------------
# 
ov_macs <-olRanges(query=makeGRangesFromDataFrame(karenina_macs), subject=makeGRangesFromDataFrame(monty_macs), output="df") 

ov_bayp <-olRanges(query=makeGRangesFromDataFrame(karenina_bayp,seqnames.field ="space"), 
                   subject=makeGRangesFromDataFrame(monty_bayp,seqnames.field ="space"), output="df")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# uv_macs<- readRDS("C:/Users/Owner/Box Sync/chip_seq/idr/uv_macs_reprod.rds")
# uv_bayp<- readRDS("C:/Users/Owner/Box Sync/chip_seq/idr/uv_bayp_reprod.rds")
uv_macs<- readRDS("data/uv_macs_reprod.rds")
uv_bayp<- readRDS("data/uv_bayp_reprod.rds")

#Tooltip function
venn_tooltip2<- function( venn ){
  venn$x$tasks[length(venn$x$tasks)+1] <- list(
    htmlwidgets::JS('
                    function(){
                    console.log("here");
                    var div = d3.select(this);
                    
                    // add a tooltip
                    var tooltip = d3.select("body").append("div")
                    .attr("class", "venntooltip")
                    .style("position", "absolute")
                    .style("text-align", "center")
                    .style("width", 250)
                    .style("height", 500)
                    .style("background", "#333")
                    .style("color","#ddd")
                    .style("padding","2px")
                    .style("border","0px")
                    .style("border-radius","8px")
                    .style("opacity",0);
                    
                    div.selectAll("path")
                    .style("stroke-opacity", 0)
                    .style("stroke", "#fff")
                    .style("stroke-width", 0)
                    
                    // add listeners to all the groups to display tooltip on mousover
                    div.selectAll("g")
                    .on("mouseover", function(d, i) {
                    
                    // sort all the areas relative to the current item
                    venn.sortAreas(div, d);
                    
                    // Display a tooltip with the current size
                    tooltip.transition().duration(400).style("opacity", .9);
                    tooltip.html(d.size+"<br>"+d.gene);
                    
                    // highlight the current path
                    var selection = d3.select(this).transition("tooltip").duration(400);
                    selection.select("path")
                    .style("stroke-width", 3)
                    .style("fill-opacity", d.sets.length == 1 ? .4 : .1)
                    .style("stroke-opacity", 1);
                    })
                    
                    .on("mousemove", function() {
                    tooltip.style("left", (d3.event.pageX) + "px")
                    .style("top", (d3.event.pageY - 28) + "px");
                    })
                    
                    .on("mouseout", function(d, i) {
                    tooltip.transition().duration(400).style("opacity", 0);
                    var selection = d3.select(this).transition("tooltip").duration(400);
                    selection.select("path")
                    .style("stroke-width", 0)
                    .style("fill-opacity", d.sets.length == 1 ? .25 : .0)
                    .style("stroke-opacity", 0);
                    });
                    }
                    ')
    )
  venn
  }
