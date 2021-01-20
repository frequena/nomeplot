library(shinydashboard)
library(shinycssloaders)
library(dplyr)
library(shinydashboardPlus)
library(shinyjs)
library(shinyWidgets)

tags$head(tags$link(rel="shortcut icon", href="favicon.ico"))
header <- dashboardHeader(
  title = "NOMePlot - A Webtool for the analysis of NOMe and Bisulfite data",
  titleWidth = 750
)

mm9 <- list('chr1' = 1, 'chr2' = 2,'chr3' = 3,'chr4' = 4,'chr5' = 5,'chr6' = 6,'chr7' = 7,'chr8' = 8,'chr9' = 9,'chr10' = 10,'chr11' = 11,'chr12' = 12,'chr13' = 13,
            'chr14' = 14,'chr15' = 15,'chr16' = 16,'chr17' = 17,'chr18' = 18,'chr19' = 19,'chrX' = 'X','chrY' = 'Y')

options(shiny.maxRequestSize=200*1024^2) 

sidebar <- dashboardSidebar(
  
  useShinyjs(),
  width = 250,
  sidebarMenu(id = 'tab',
              menuItem("Welcome", tabName = "welcome", icon = icon("info-circle")),
              menuItem("Sanger sequencing", tabName = "", icon = icon("flask"),
                       menuSubItem("BS-PCR", tabName = "bisulfite", icon = icon("dna")),
                       menuSubItem("NOMe-PCR", tabName = "nome", icon = icon("dna"))),
              menuItem("High-throughput sequencing", tabName = "", icon = icon("flask"),
                       menuSubItem("BS-seq", tabName = "gw_bi", icon = icon("dna")),
                      menuSubItem("NOMe-seq", tabName = "gw", icon = icon("dna")))
              
           
  ),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  hr(),
  # email sharing link
  menuItem("Feedback & suggestion", icon = icon("envelope-o"),
           href = "mailto:francisco.requena@institutimagine.org, davidlandeira@ugr.es"),
  br(),
  menuItem("Questions & discussion", icon = icon("users"),
           href = "https://groups.google.com/forum/#!forum/nomeplot"),
  br(),
  menuItem("The code is available at Github", icon = icon("github"),
           href = "https://github.com/frequena/nomeplot"))
  

 body <-  dashboardBody(
   
   setShadow("box"),
   tabItems(
     tabItem(tabName = "welcome",
             fluidRow(
               column(width = 8,
                      tabBox(title = tagList(shiny::icon("gear"), "Welcome to NOMePlot"), width = 200, side = "right", 
                             selected = tags$b("Info"),
                             
                             tabPanel(tags$b("Contact us..."),
                                      tags$br(),
                                      tags$b('If you have questions or comments, please write to:'),
                                      tags$br(),
                                      tags$br(),
                                      p(tags$b('- Francisco Requena'), '(francisco.requena@institutimagine.org)'),
                                      p(tags$b('- David Landeira'), '(davidlandeira@ugr.es)'),
                                      tags$br(),
                                      p('If you want to join', tags$b('the NOMePlot Google group:')),
                                      p('1. Sign in to Google Groups.'),
                                      p('2. In the box at the top, enter the subject "nomeplot".'),
                                      p('3. To join the group, click Join group or Apply to join group.')
                             ),
                                      
                                      
                            
                             tabPanel(tags$b("FAQ"), 
                                      
                                      
                                      tags$b('-	What if my genome assembly is not available?'),
                                      tags$p('Send us an email with the genome assembly that you need and we will update the app with this new option.'),
                                    tags$b('-	What inputs does NOMePlot accept?'),
                                    tags$p('The program accepts nucleic acid sequences with FASTA format (Sanger sequencing section) and BAM format (High-throughput sequencing section). Please, make sure that your FASTA files have one of the next terminations: .seq, .fasta, .txt'),
                                    tags$b('-	Is there a limit on the number of sequences (Sanger sequencing section) that I upload to NOMePlot?'),
                                    tags$p('Yes, the limit is 60 sequences files.'),
                                    tags$b('-	Is there a limit on the size of the BAM file that I submit to NOMePlot?'),
                                    tags$p('Yes, the maximum size allowed is 50MB. If your file is bigger, you can ask your bioinformatician to use samtools to reduce the size: samtools view -b your_file.bam "Chr1:start-end" > output.bam'),
                                    tags$b('-	Can I download a .png image from the graphics?'),
                                    tags$p('Yes, you can download it just making right-click on the image and “Save image as…”. You can also download a report (.html format) with all the plots of your data, just make click on “Download NOMe/Bisulfite/Genome wide analysis” button that you will find in every section.'),
                                    tags$b('-	What is the .svg format?'),
                                    tags$p('It is a vector image format and it can be useful for your lollipop graphics for two reasons: 1) Independently of the scale that you apply, your image will conserve the same resolution. 2) Every element of the lollipop can be modified independently. These two elements confer publication-ready quality to your lollipops.'),
                                    tags$b('-	How can I open a lollipop plot with .svg format?'),
                                    tags$p('A simple web browser such as Chrome can be used to visualize your plot. But if you want to modify the plot, you need to use another program such as Inkscape. It is a free and open-source vector-graphic editor, that you can download directly from their website (inkscape.org)'),
                                    tags$b('-	How can I upload the sample data?'),
                                    tags$p('Please, follow the instruction that you will find on the compressed file, or go to the section “Instructions to upload the sample data” in this page.'),
                                    tags$b('- Why does the application turn gray when I am using it?'),
                                    tags$p('NOMePlot browser sessions running for more than 15 minutes without activity  are terminated.'),
                                    tags$b('- Why does I lose my uploaded files?'),
                                    tags$p('To keep a good performance of the app, each time that you select a different kind of analysis, the application
                                           remove the input data submitted by the user.'),
                                    tags$b('- Can I run NOMePlot locally?'),
                                    tags$p('Yes, the code of NOMePlot is stored in a a git repository on GitHub, therefore you can run the app directly. From a R session and with the shiny library installed, you can run the following command: shiny::runGitHub("#####", "####").'),
                                    tags$b('-	Where can I get further help for NOMePlot?'),
                                    tags$p('We recommend that you use the Google Group available to ask your questions. Besides, you can contact us by email. You can find both the link of the group as the email address in the sidebar of the app or in the “Contact us” section on the this page.')
                                    
  
                                      
                                      ),
                             
                             

                             
                             tabPanel(tags$b("Info"),
                                               fluidRow(
                                                 column(width = 12,
                                                        p("NOMePlot is a webtool that facilitates the analysis of Nucleosome Occupancy and Methylome (NOMe) and Bisulfite (BS) sequencing data generated using Sanger or High-throughput sequencing (NOMe-PCR, NOMe-seq, BS-PCR and BS-seq)."),
                                                        
                                                        p("This web application is easy to use and does not required bioinformatic skills. Input sequences for sanger sequencing analysis must be provided in fasta format while input data for high-throughput sequencing must be provided as a bam file. The user-interface is self-explanatory."),
                                                        
                                                        p("NOMePlot has been developed with R/Bioconductor packages and Shiny. Source code for local installation is available from GitHub."),
                                                          p("NOMePlot operates under the",tags$b("GNU General Public License version 3.")),
                                                        tags$br(),
                                                        tags$br(),
                                                        p(tags$b('Citing NOMePlot:')),
                                                        p('Requena F, Asenjo HG, Barturen G, Martorell-Marugán J, Carmona-Sáez P and David Landeira.', 
                                                          tags$b('NOMePlot:  analysis of DNA methylation and nucleosome occupancy at the single molecule.'), 'Scientific Reports. May 31;9(1):8140. doi: 10.1038/s41598-019-44597-2.  2019.', tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/31148571", "Abstract")),
                                                        
                                                        
                                                        column(12,  align="center",
                                                               tags$hr(),
                                                        tags$img(src='home.png')),
                                                        tags$br(),
                                                        tags$br(),
                                                        tags$b('Current version: 0.80'),
                                                        tags$br(),
                                                        tags$br()
                                                        

                                                        )

))

)),


column(width = 4,
       
       
       
       box(title = tagList(shiny::icon("download"), "Download BS-PCR data"), width = 200, 
           column(width = 12, 
                  fluidRow(
                    tags$b("Please, click the box below if you want to download BS-PCR data to check NOMePlot functionalities:"),
                    tags$br(),
                    tags$hr(),
                    column(width = 12, align="center",
                           
                           downloadButton('sample_seq3', 'Sample data obtained by bisulfite sanger sequencing',
                                          style="color: #fff; width: 360px; background-color: #337ab7; border-color: #2e6da4;
                                          align: 'center'"))))),
       
       box(title = tagList(shiny::icon("download"), "Download NOMe-PCR data"), 
           width = 200,
           column(width = 12, 
                  fluidRow(
                    tags$b("Please, click the box below if you want to download example NOMe-PCR sequences to check NOMePlot functionalities:"),
                    tags$br(),
                    tags$hr(),
                    column(width = 12, align="center",
                           
                           downloadButton('sample_seq', 'Download sample data obtained by NOMe-PCR',
                                          style="color: #fff; width: 330px; background-color: #337ab7; border-color: #2e6da4;
                                                        align: 'center'"))))),
       box(title = tagList(shiny::icon("download"), "Download NOMe-seq data (15MB)")
           , width = 200, 
           column(width = 12, 
                  fluidRow(
                    tags$b("Please, click the box below if you want to download example NOMe-seq data to check NOMePlot functionalities:"),
                    tags$br(),
                    tags$hr(),
                    column(width = 12, align="center",
                           
                           downloadButton('sample_seq2', 'Download sample data obtained by NOMe-seq',
                                          style="color: #fff; width: 330px; background-color: #337ab7; border-color: #2e6da4;
                                   align: 'center'")))))
       
       
       
                           )





)),

   
     tabItem(tabName = "nome",
    fluidRow(
      column(width = 3,
             box(title = "Input sequences", width = 200, solidHeader = TRUE, status = "primary",
                 fileInput("file_reference", label = h5(strong("Reference genomic sequence")), 
                           accept = c(".seq", "text", ".fasta", ".txt")),
                 fileInput("files_problem", label = h5(strong("NOMe-treated sequences")), multiple = TRUE,
                           accept = c(".seq", "text", ".fasta", ".txt")),
                 radioButtons("primers", label = h5(strong("Define region of interest using primer sequences")),
                              choices = list("Yes" = 1,  "No" = 2), 
                              selected = 2),
                 conditionalPanel(
                   condition = "input.primers == 1" ,
                   textInput("fw_primer", label = ("") , value = "Forward primer bisulfite converted"),
                   textInput("rv_primer", label = ("") , value = "Reverse primer bisulfite converted"))),
             box(title = "Arrange sequences", width = 200, solidHeader = TRUE, status = "primary",
                 radioButtons("type_arrang", label = h5(strong("Define criterion of arrangement")),
                              choices = list("Use GpC sites" = 1,  "Use CpG sites" = 2), 
                              selected = 1),
                 conditionalPanel(
                  
                 condition = "input.type_arrang == 1" ,
                 sliderInput("arrang_gpc", 
                             label = h5(strong("GpC clustering window (nucleotide length)")),
                             min = 20, max = 160, value = 60, step = 10)),
                 conditionalPanel(
                   
                   condition = "input.type_arrang == 2" ,
                   sliderInput("arrang_cpg", 
                               label = h5(strong("CpG clustering window (nucleotide length)")),
                               min = 20, max = 160, value = 60, step = 10)))),
      column(width = 9, solidHeader = TRUE,
             box(width = NULL, solidHeader = TRUE,
                 tabsetPanel(type = "tabs", 
                             tabPanel("Lollipop graphic",  plotOutput("instant_plot") %>% withSpinner()),
                             tabPanel("Heatmap",  plotOutput("heatmap_web") %>% withSpinner()),
                             tabPanel("Line Chart",  plotOutput("plot_window_iowa") %>% withSpinner()))
             ),
             fixedRow(
               column(width = 12, solidHeader = TRUE,
             box(title = "Personalize Lollipop graphic",
                 width = 12, solidHeader = TRUE, status = "primary",
                 column(width = 3, solidHeader = TRUE,
                        
             radioButtons("red_lines", label = h5(strong("Identify nucleosome occupied regions")),
                          choices = list("Yes" = TRUE,  "No" = 2), 
                          selected = 2),
             conditionalPanel(
               condition = "input.red_lines != 2",
               sliderInput("size_draw",
                           label = h5(strong("Region size length")),
                           min = 100, max = 160, value = 140, step = 10),
               sliderInput("mean_nucleosome",
                           label = h5(strong("Stringency")),
                           min = 0.5, max = 1, value = 0.7, step = 0.1))),
               
             column(width = 2, solidHeader = TRUE,
             radioButtons("tss", label = h5(strong("Label TSS")),
                          choices = list("Yes" = 1,  "No" = 2), 
                          selected = 2),
             conditionalPanel(
               condition = "input.tss == 1" ,
               numericInput("num_tss", label = h5(strong("Position of TSS")), value = 0))),
             column(width = 3, solidHeader = TRUE,
               selectInput("color_circles", label =  h5(strong("Nucleosome lines color")), 
                           choices = list("Red" = "red", "Blue" = "blue", "Black" = "black",
                                          "Cyan" = "cyan", "Magenta" = "magenta", "Green" = "green"), 
                           selected = 1)),
              
               
               
             column(width = 4, solidHeader = TRUE,
                    selectInput("color_lines", label = h5(strong("GpC circles color")), 
                                choices = list("Blue" = "blue", "Red" = "red", "Black" = "black",
                                               "Cyan" = "cyan", "Magenta" = "magenta", "Green" = "green"), 
                                selected = 1),
                    tags$hr(),
               downloadButton('downloadData_nome', 'Download NOMe analysis'),
               tags$hr(),
               downloadButton('downloadData_graphic', 'Download Lollipop Graphic (.svg)'),
               tags$hr(),
               downloadButton('downloadData_technical', 'Download technical report'))
               ))))              
      )
  ),
  
  tabItem(tabName = "gw",
          fluidRow(
            column(width = 3,
                   box(title = "Reference sequence", width = 200, solidHeader = TRUE, status = "primary",
                      
                         
                       selectInput("file_reference_gw", label =  h5(strong("Select genomic assembly")), 
                                  
                                   choices = list_genomes,
                                   selected = 'mm9'),
                       selectInput("chr_select", label =  h5(strong("Choose chromosome")), 
                                   choices = list_genomes),
                       textInput("start_interval", label = ("Genomic Interval - Start (max. 2.000 nt. distance)") , value = 0), 
                       textInput("end_interval", label = ("Genomic Interval - End") , value = 0),
                       actionButton("run", label = "Run analysis", icon("paper-plane"), style="width: 270px")),
                       
                   box(title = "Problem sequences", width = 200, solidHeader = TRUE, status = "primary",
                    
                         fileInput("bam_gw", label = h5(strong("Load BAM file (max. 50MB)")), 
                                   accept = c(".bam"))),
                   
                   box(title = "Arrange sequences", width = 200, solidHeader = TRUE, status = "primary",
                       radioButtons("type_arrang_gw", label = h5(strong("Define criterion of arrangement")),
                                    choices = list("Use GpC sites" = 1,  "Use CpG sites" = 2), 
                                    selected = 1),
                       conditionalPanel(
                         
                         condition = "input.type_arrang_gw == 1" ,
                         sliderInput("arrang_gpc_gw", 
                                     label = h5(strong("GpC clustering window (nucleotide length)")),
                                     min = 20, max = 160, value = 60, step = 10)),
                       conditionalPanel(
                         
                         condition = "input.type_arrang_gw == 2" ,
                         sliderInput("arrang_cpg_gw", 
                                     label = h5(strong("CpG clustering window (nucleotide length)")),
                                     min = 20, max = 160, value = 60, step = 10))),
               
                 

                   
                  box(title = "Discards SNPs in the region with a VCF file?", width = 200, solidHeader = TRUE, status = "primary",
                       radioButtons("snp_gw", label = h5(strong("")),
                                    choices = list("No" = 1,  "Yes" = 2) , 
                                                   
                                    selected = 1),
                       conditionalPanel(
                         condition = "input.snp_gw == 2" ,
                         fileInput("file_vcf", label = h5(strong("Load VCF file")), 
                                   accept = c(".seq", "text", ".fasta", ".txt")))),
                  
                  
                  box(title = "Download", width = 200, solidHeader = TRUE, status = "primary",
                      downloadButton('downloadData_gw', 'Download Genome Wide analysis'),
                      tags$hr(),
                      downloadButton('downloadData_graphic_gw', 'Download Lollipop Graphic (.svg format)')
                      
                      
                  )
                   
                  ),
            column(width = 9, solidHeader = TRUE,
                   box(width = 300, height = 500, solidHeader = TRUE,
                       tabsetPanel(type = "tabs", 
                                   tabPanel("Genomic Interval", plotOutput("genomic_plot") %>% withSpinner()),
                                   tabPanel("GpC values",  plotOutput("lollipop_gw_gpc") %>% withSpinner()),
                                   tabPanel("CpG values",  plotOutput("lollipop_gw_cpg") %>% withSpinner()))
                      
                   )),
                   column(width = 9, solidHeader = TRUE,
                          box(width = NULL, height = 500, solidHeader = TRUE,
                              tabsetPanel(type = "tabs", 
                                          
                                          tabPanel("Overlap",  plotOutput("plot_combined") %>% withSpinner()),
                                          tabPanel("1 - GpC methylation (%)",  plotOutput("plot_gpc_gw") %>% withSpinner()),
                                          tabPanel("CpG methylation (%)",  plotOutput("plot_cpg_gw") %>% withSpinner())))
                   
            
                   ))),
 
  tabItem(tabName = "bisulfite",
          fluidRow(
            column(width = 3,
                   box(title = "Input sequences", width = 200, solidHeader = TRUE, status = "primary",
                       fileInput("file_reference_bi", label = h5(strong("Reference genomic sequence")), 
                                 accept = c(".seq", "text", ".fasta", ".txt")),
                       fileInput("files_problem_bi", label = h5(strong("Bisulfite-treated sequences")), multiple = TRUE,
                                 accept = c(".seq", "text", ".fasta", ".txt")),
                       radioButtons("primers_bi", label = h5(strong("Define region of interest using primer sequences")),
                                    choices = list("Yes" = 1,  "No" = 2), 
                                    selected = 2),
                       conditionalPanel(
                         condition = "input.primers_bi == 1" ,
                         textInput("fw_primer_bi", label = ("") , value = "Forward primer bisulfite converted"),
                         textInput("rv_primer_bi", label = ("") , value = "Reverse primer bisulfite converted"))),
                   box(title = "Arrange sequences", width = 200, solidHeader = TRUE, status = "primary",
                       sliderInput("wd_size_bi", 
                                   label = h5(strong("CpG clustering window (nucleotide length)")),
                                   min = 20, max = 160, value = 60, step = 10))),
            column(width = 9, solidHeader = TRUE,
                   box(width = NULL, solidHeader = TRUE,
                       tabsetPanel(type = "tabs", 
                                   tabPanel("Lollipop graphic",  plotOutput("instant_plot_bi") %>% withSpinner()),
                                   tabPanel("Heatmap",  plotOutput("heatmap_web_bi") %>% withSpinner()),
                                   tabPanel("Line Chart",  plotOutput("plot_window_iowa_bi") %>% withSpinner()))
                   ),
                   fixedRow(
                     column(width = 12, solidHeader = TRUE,
                            box(title = "Personalize Lollipop graphic",
                                width = 12, solidHeader = TRUE, status = "primary",
                                
                                
                                column(width = 2, solidHeader = TRUE,
                                       radioButtons("tss_bi", label = h5(strong("Label TSS")),
                                                    choices = list("Yes" = 1,  "No" = 2), 
                                                    selected = 2),
                                       conditionalPanel(
                                         condition = "input.tss_bi == 1" ,
                                         numericInput("num_tss_bi", label = h5(strong("Position of TSS")), value = 0))),
                                
                                
                                
                                column(width = 2, solidHeader = TRUE,
                                       selectInput("color_circles_bi", label = h5(strong("CpG circles color")), 
                                                   choices = list("Black" = "black", "Red" = "red", "Blue" = "blue",
                                                                  "Cyan" = "cyan", "Magenta" = "magenta", "Green" = "green"), 
                                                   selected = 1)),
                                column(width = 4, solidHeader = TRUE,
                                       tags$hr(),
                                       downloadButton('downloadData_nome_bi', 'Download Bisulfite analysis'),
                                       tags$hr(),
                                       downloadButton('downloadData_graphic_bi', 'Download Lollipop Graphic (.svg)'),
                                       tags$hr(),
                                       downloadButton('downloadData_technical_bi', 'Download technical report'))
                            )))))),
  
  tabItem(tabName = "gw_bi",
          fluidRow(
            column(width = 3,
                   box(title = "Reference sequence", width = 200, solidHeader = TRUE, status = "primary",
                       
                       
                       selectInput("file_reference_gw_bi", label = h5(strong("Select genomic assembly")),
                                    choices = list_genomes,
                                    selected = 'mm9'),
                       selectInput("chr_select_bi", label =  h5(strong("Choose chromosome")), 
                                   choices = list_genomes),
                       textInput("start_interval_bi", label = ("Genomic Interval - Start (max. 2.000 nt. distance)") , value = 0),
                       textInput("end_interval_bi", label = ("Genomic Interval - End") , value = 0),
                       actionButton("run_bi", label = "Run analysis", icon("paper-plane"), style="width: 270px")),
                   
                   box(title = "Problem sequences", width = 200, solidHeader = TRUE, status = "primary",
                       
                       fileInput("bam_gw_bi", label = h5(strong("Load BAM file (max. 50MB)")), 
                                 accept = c(".bam"))),
                   
                   box(title = "Arrange sequences", width = 200, solidHeader = TRUE, status = "primary",
                       radioButtons("type_arrang_gw_bi", label = h5(strong("Define criterion of arrangement")),
                                    choices = list("Use CpG sites" = 1), 
                                    selected = 1),
                       sliderInput("arrang_cpg_gw_bi", 
                                   label = h5(strong("CpG clustering window (nucleotide length)")),
                                   min = 20, max = 160, value = 60, step = 10)),
                   
                   
                   
                   
                   box(title = "Discards SNPs in the region with a VCF file?", width = 200, solidHeader = TRUE, status = "primary",
                       radioButtons("snp_gw_bi", label = h5(strong("")),
                                    choices = list("No" = 1,  "Yes" = 2) , 
                                    
                                    selected = 1),
                       conditionalPanel(
                         condition = "input.snp_gw_bi == 2" ,
                         fileInput("file_vcf_bi", label = h5(strong("Load VCF file")), 
                                   accept = c(".seq", "text", ".fasta", ".txt")))),
                   
                   
                   box(title = "Download", width = 200, solidHeader = TRUE, status = "primary",
                       downloadButton('downloadData_gw_bi', 'Download Genome Wide analysis'),
                       tags$hr(),
                       downloadButton('downloadData_graphic_gw_bi', 'Download Lollipop Graphic (.svg format)')
                       
                       
                   )
                   
            ),
            column(width = 9, solidHeader = TRUE,
                   box(width = NULL, height = 500, solidHeader = TRUE,
                       tabsetPanel(type = "tabs", 
                                   tabPanel("Genomic Interval",  plotOutput("genomic_plot_bi") %>% withSpinner()),
                                   tabPanel("CpG values",  plotOutput("lollipop_gw_cpg_bi") %>% withSpinner()))
                       
                   )),
            column(width = 9, solidHeader = TRUE,
                   box(width = NULL, height = 500, solidHeader = TRUE,
                       tabsetPanel(type = "tabs", 

                                   tabPanel("Methylation of CpG sites",  plotOutput("plot_cpg_gw_bi") %>% withSpinner()
                                            )))
                   
                   
            )))
  ))

dashboardPage(
  header,
  sidebar,
  body
)
