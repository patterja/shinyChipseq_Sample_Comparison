
shinyUI(navbarPage("ChIPseq Sample Comparison",
#----------------------------------------------------------------Between Tumors
tabPanel("MACS2",
  fluidRow(
     column(5,plotOutput("plt_gibb1_macs", height=350)),
     column(5,offset=0.5,plotOutput("plt_gibb2_macs", height=350))),
  fluidRow(
    column(5,offset=1,sliderInput('depth1_macs', 'Minimum Depth of Coverage threshold',
                       min = min(k_macsjoin$coverage, na.rm=TRUE),max = max(k_macsjoin$coverage, na.rm=TRUE),
                       value = c(min(k_macsjoin$coverage, na.rm=TRUE),max(k_macsjoin$coverage, na.rm=TRUE))),
           submitButton("Submit")), 
    column(5,sliderInput('depth2_macs', 'Minimum Depth of Coverage threshold',
                         min = min(m_macsjoin$coverage, na.rm=TRUE),max = max(m_macsjoin$coverage, na.rm=TRUE),
                         value = c(min(m_macsjoin$coverage),max(m_macsjoin$coverage))), 
           submitButton("Submit"))),
  fluidRow(
    column(6, d3vennROutput("macs_venn")),
    column(6, tableOutput("macsstats"))),
  tabsetPanel(
    tabPanel("Gibbon1 MAC2",dataTableOutput("gibbon1_macs")),
    tabPanel("Gibbon2 MACS2",dataTableOutput("gibbon2_macs")), 
    tabPanel("Overlaps", dataTableOutput("ov_macstable")))
  ),
tabPanel("BayesPeak",
         fluidRow(
           column(5,plotOutput("plt_gibb1_bayp"),
                  sliderInput('depth1_bayp', 'Minimum Depth of Coverage threshold',
                              min = min(k_baypjoin$coverage, na.rm=TRUE),max = max(k_baypjoin$coverage, na.rm=TRUE),
                              value = c(min(k_baypjoin$coverage, na.rm=TRUE),max(k_baypjoin$coverage, na.rm=TRUE))),
                  submitButton("Submit")),
           column(5,plotOutput("plt_gibb2_bayp"), 
                  sliderInput('depth2_bayp', 'Minimum Depth of Coverage threshold',
                              min = min(m_baypjoin$coverage, na.rm=TRUE),max = max(m_baypjoin$coverage, na.rm=TRUE),
                              value = c(min(m_baypjoin$coverage),max(m_baypjoin$coverage))), 
                  submitButton("Submit"))
         ), 
         fluidRow(
           column(6, d3vennROutput("bayp_venn")),
           column(6, tableOutput("baypstats"))),
         tabsetPanel(
           tabPanel("Gibbon1",dataTableOutput("gibbon1_bayp")),
           tabPanel("Gibbon2",dataTableOutput("gibbon2_bayp")),
           tabPanel("Overlaps",dataTableOutput("ov_bayptable")))
)
))
