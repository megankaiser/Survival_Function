# SURVIVAL FUNCTION 
# OCTOBER 2018

#===========================================================================================#
# STEP 1: Mount Isilon if doing a CoMMpass Survival Analysis                                #
#===========================================================================================#

#===========================================================================================#
# STEP 2: Load required libraries for survival analysis and file merging                    #
#===========================================================================================#
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('muhaz')) install.packages('muhaz'); library('muhaz')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survMisc')) install.packages('survMisc'); library('survMisc')
if (!require('survminer')) install.packages('survminer'); library('survminer')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('plotrix')) install.packages('plotrix'); library('plotrix')
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')  
#===========================================================================================#  

#===========================================================================================#
# STEP 3: Compile the survival function                                                     #
#===========================================================================================#
Survival_Function <- function(COMMPASS_SURVIVAL_ANALYSIS,IA_RELEASE,  RENAME_PATIENT_ID_HEADER = c("public_id"), RENAME_OS_TIME_HEADER = c("ttcos"),
                              RENAME_OS_CENSOR_HEADER = c("censos"), RENAME_PFS_TIME_HEADER = c("ttcpfs"), RENAME_PFS_CENSOR_HEADER = c("censpfs"),
                              YOUR_CLINICAL_FILE_PATH, NEW_STRATA, STRATA_VARIABLE, NEW_STRATA_FILE, STRATA_FILE_PATH, STRATIFY_FROM_EXISTING_DATA_FRAME, YOUR_DATA_FRAME, 
                              STRATA_PT_ID_COL_HEADER= c("Patient_ID"),  SUBSET, SUBSET_HEADER = c(""), SUBSET_VALUE, PFS_OS= c(""), DEFAULT_LEGEND_LABELS, 
                              CUSTOM_LEGEND_LABELS = c(""), DEFAULT_TITLE, CUSTOM_TITLE = c(""), X_LABEL = c("Days"), Y_LABEL = c("Survival Probability"),
                              RISK_TABLE, RISK_TABLE_HEIGHT = 0.25, RISK_TABLE_FONT = 5, PVAL_MATRIX, PRINT_PVAL_MATRIX, SAVE_PVAL_MATRIX, P_VALUE, 
                              PVAL_HEADER_SIZE = 1.0, PVAL_FONT_SIZE = 0.75,CONFIDENCE_INTERVAL, CENSOR_TICK, MEDIAN_SURVIVAL_LINE, PLOT_ANNOTATIONS_CALL,
                              PLOT_ANNOTATIONS_TEXT = c(""),PLOT_ANNOTATIONS_X = 250, PLOT_ANNOTATIONS_Y = 0.15, PLOT_ANNOTATIONS_SIZE = 5, OS_TABLE, PFS_TABLE, 
                              OS_COXPH_TABLE, PFS_COXPH_TABLE, SAVE_PLOT, FILE_NAME = c("")) {
#===========================================================================================#
  
#===========================================================================================#  
# SCROLL TO STEP 4 (~ line 500)
#===========================================================================================#   
  
#===========================================================================================#   
# TEST VARIABLES
  # COMMPASS_SURVIVAL_ANALYSIS = T
  # IA_RELEASE = 14 # Enter CoMMpass release number
  # 
  # YOUR_CLINICAL_FILE_PATH = c("/Volumes/MMRF/commpass/clinical_flat_files/IA14/CoMMpass_IA14_FlatFiles/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL.csv") # not necessary for coMMpass
  # RENAME_PATIENT_ID_HEADER = c("public_id") # insert your patient id variable
  # RENAME_OS_TIME_HEADER = c("ttt2line") # insert your overall survival time variable
  # RENAME_OS_CENSOR_HEADER = c("censt2line") # insert your overall survival flag variable
  # RENAME_PFS_TIME_HEADER = c("") # insert your progression free survival time variable
  # RENAME_PFS_CENSOR_HEADER = c("") # insert your profression free survival flag variable
  # 
  # NEW_STRATA = T # T or F
  # STRATA_VARIABLE = c("D_PT_race") #insert the variable you want to stratify by
  # STRATA_PT_ID_COL_HEADER= c("PUBLIC_ID")
  # NEW_STRATA_FILE= T
  # STRATA_FILE_PATH = c("/Volumes/MMRF/commpass/clinical_flat_files/IA14/CoMMpass_IA14_FlatFiles/MMRF_CoMMpass_IA14_PER_PATIENT.csv") # not needed if you want to compare info from the survival table
  # STRATIFY_FROM_EXISTING_DATA_FRAME = F
  # YOUR_DATA_FRAME = df
  # 
  # SUBSET = F  # T or F
  # SUBSET_HEADER = c("") # if Subset = T must be filled in
  # SUBSET_VALUE = 1 # insert a number if subset = T
  # 
  # PFS_OS= c("PFS")   # Set the type of survival you want plotted: PFS = Progression-free survival, OS = Overall survival
  # 
  # DEFAULT_LEGEND_LABELS = T # T or F
  # CUSTOM_LEGEND_LABELS = c("")
  # DEFAULT_TITLE = F # T or F
  # CUSTOM_TITLE = c("TITLE")
  # X_LABEL = c("Days") #can be days, months, years, etc
  # Y_LABEL = c("Survival Probability") # specify Y axis label
  # 
  # RISK_TABLE = T # T or F
  # RISK_TABLE_HEIGHT = 0.25
  # RISK_TABLE_FONT = 5
  # 
  # PVAL_MATRIX = T # T or F
  # PRINT_PVAL_MATRIX = T # T or F
  # SAVE_PVAL_MATRIX = F # T or F
  # 
  # P_VALUE = T # T or F
  # PVAL_HEADER_SIZE = 1.0 # DEFAULT = 1
  # PVAL_FONT_SIZE = 0.75 # DEFAULT =.75
  # 
  # CONFIDENCE_INTERVAL = F # T or F
  # CENSOR_TICK = F #  T or F
  # MEDIAN_SURVIVAL_LINE = T  # T or F
  # 
  # PLOT_ANNOTATIONS_CALL = F # T or F
  # PLOT_ANNOTATIONS_TEXT = c("")
  # PLOT_ANNOTATIONS_X = 250  # DEFAULT = 250
  # PLOT_ANNOTATIONS_Y = 0.15 # DEFAULT = .15
  # PLOT_ANNOTATIONS_SIZE = 5 # DEFAULT= 5
  # OS_TABLE = T # T or F
  # PFS_TABLE = T # T or F
  # OS_COXPH_TABLE = T  # T or F
  # PFS_COXPH_TABLE = T  # T or F
  # SAVE_PLOT  = F # T or F
  # FILE_NAME = c("")  # name the file if SAVE_PLOT= T
#===========================================================================================#


  # User provided column headers list
  RENAME_HEADERS <- list( "Patient_ID" = RENAME_PATIENT_ID_HEADER, 'Time_To_Censored_OS' = RENAME_OS_TIME_HEADER, 'OS_Censor_Flag' = RENAME_OS_CENSOR_HEADER,'Time_To_Censored_PFS' =  RENAME_PFS_TIME_HEADER, "PFS_Censor_Flag" = RENAME_PFS_CENSOR_HEADER )
  RENAME_HEADERS <- unlist(RENAME_HEADERS)
  
  # Handles non-CoMMpass analysis and renames strata patient ID to patient id
  if(COMMPASS_SURVIVAL_ANALYSIS == F && NEW_STRATA_FILE == F && STRATA_PT_ID_COL_HEADER == c("")){
    STRATA_PT_ID_COL_HEADER <- RENAME_PATIENT_ID_HEADER}
  
  
# CLINICAL FILE SECTION
     # COMMPASS ANALYSIS 
     #reads in the stand alone survival table, then renames columns
    if(COMMPASS_SURVIVAL_ANALYSIS == T){
        IA <- as.character(IA_RELEASE)
        CLINICAL_DATA_MATRIX_ORIG<- read.csv(paste("/Volumes/MMRF/commpass/clinical_flat_files/IA",IA,"/CoMMpass_IA",IA,"_FlatFiles/MMRF_CoMMpass_IA",IA,"_STAND_ALONE_SURVIVAL.csv", sep = ""), header = T, na.strings=c("","NA"))
        CLINICAL_DATA_MATRIX <- subset(CLINICAL_DATA_MATRIX_ORIG, select = c('public_id', 'ttcos', 'censos', 'ttcpfs', 'censpfs'))
        colnames(CLINICAL_DATA_MATRIX) <- c('Patient_ID', 'Time_To_Censored_OS', 'OS_Censor_Flag', 'Time_To_Censored_PFS', 'PFS_Censor_Flag')
      # Handles stratifying from the stand alone survival table without reading in again
      if(NEW_STRATA==T && NEW_STRATA_FILE == F && STRATIFY_FROM_EXISTING_DATA_FRAME == F){
       CLINICAL_DATA_MATRIX_ORIG <- subset(CLINICAL_DATA_MATRIX_ORIG, select = c(STRATA_VARIABLE, 'public_id'))}
    # NON-COMMPASS ANALYSIS
    }else if (COMMPASS_SURVIVAL_ANALYSIS == F){
        end <-  str_sub(YOUR_CLINICAL_FILE_PATH, start=-3)
          # read in clinical file path
          if(end == "csv"){
            CLINICAL_DATA_MATRIX_ORIG <- read.csv(YOUR_CLINICAL_FILE_PATH, header = T, na.strings=c("", "NA"))
          }else if(end == "txt"){
            CLINICAL_DATA_MATRIX_ORIG <- read.csv(YOUR_CLINICAL_FILE_PATH, header = T, sep = "\t")
          }
        # rename all 5 headers (PT_ID, OS TIME, OS CENSOR, PFS TIME, PFS CENSOR)
        if(!is.null(RENAME_PATIENT_ID_HEADER) && !is.null(RENAME_OS_TIME_HEADER) && !is.null(RENAME_OS_CENSOR_HEADER) && (RENAME_PFS_TIME_HEADER != c("") | RENAME_PFS_CENSOR_HEADER != c(""))){
          CLINICAL_DATA_MATRIX <- subset(CLINICAL_DATA_MATRIX_ORIG, select = c(RENAME_PATIENT_ID_HEADER, RENAME_OS_TIME_HEADER, RENAME_OS_CENSOR_HEADER,  RENAME_PFS_TIME_HEADER, RENAME_PFS_CENSOR_HEADER) )
          colnames(CLINICAL_DATA_MATRIX) <- names(RENAME_HEADERS)
        # only renames PT_ID, OS TIME AND OS CENSOR
        }else if(!is.null(RENAME_PATIENT_ID_HEADER) && !is.null(RENAME_OS_TIME_HEADER) && !is.null(RENAME_OS_CENSOR_HEADER) && (RENAME_PFS_TIME_HEADER == c("")  | RENAME_PFS_CENSOR_HEADER == c(""))){
          CLINICAL_DATA_MATRIX <- subset(CLINICAL_DATA_MATRIX_ORIG, select = c(RENAME_PATIENT_ID_HEADER, RENAME_OS_TIME_HEADER, RENAME_OS_CENSOR_HEADER))
          colnames(CLINICAL_DATA_MATRIX) <- names(unlist(RENAME_HEADERS[1:3]))
      # Stratifies from the clinical file path without reading file in again.
      }else if (NEW_STRATA == T && NEW_STRATA_FILE == F && STRATIFY_FROM_EXISTING_DATA_FRAME == F){
        CLINICAL_DATA_MATRIX_ORIG <- subset(CLINICAL_DATA_MATRIX_ORIG, select = c(STRATA_VARIABLE, RENAME_PATIENT_ID_HEADER))
        colnames(CLINICAL_DATA_MATRIX_ORIG)[ grep( RENAME_PATIENT_ID_HEADER, colnames(CLINICAL_DATA_MATRIX_ORIG) ) ] <- STRATA_PT_ID_COL_HEADER 
      } 
    }
      
    
 
# STRATA SECTION
    if ( NEW_STRATA == T){
      if(NEW_STRATA_FILE == F && STRATIFY_FROM_EXISTING_DATA_FRAME == F){
        # IF STRATIFYING FROM STAND ALONE SURVIVAL, ADD STRATA FIELD
        if(COMMPASS_SURVIVAL_ANALYSIS==T){
          CLINICAL_DATA_MATRIX <- merge(CLINICAL_DATA_MATRIX, CLINICAL_DATA_MATRIX_ORIG, by.x="Patient_ID", by.y=STRATA_PT_ID_COL_HEADER)
        # IF STRAFIFYING FROM CLINICAL FILE PATH, ADD STRATA FIELD. 
         }else if (COMMPASS_SURVIVAL_ANALYSIS==F){
          CLINICAL_DATA_MATRIX <- merge(CLINICAL_DATA_MATRIX, CLINICAL_DATA_MATRIX_ORIG, by.x="Patient_ID", by.y=STRATA_PT_ID_COL_HEADER)
         }
      # IF A STRATA FILE PATH IS ENTERED, READ IN FILE AND STRATIFY FROM THAT FILE
      }else if (NEW_STRATA_FILE == T && STRATIFY_FROM_EXISTING_DATA_FRAME == F){
          var <-  str_sub(STRATA_FILE_PATH, start=-3)
        if (var == "csv"){
          STRATA_FILE <- read.csv(STRATA_FILE_PATH, header = T, na.strings= c("", "NA"))
        }else if ( var == "txt"){
          STRATA_FILE <- read.csv(STRATA_FILE_PATH, header = T, sep = "\t")
        }
      # MERGE STAND ALONE SURVIVAL TABLE WITH THE STRATA FILE
      if ( COMMPASS_SURVIVAL_ANALYSIS == T){
        if (STRATA_PT_ID_COL_HEADER == c("")){
          colnames(STRATA_FILE)[ grep( "public_id|PUBLIC_ID", colnames(STRATA_FILE) ) ] <- "public_id"
          CLINICAL_DATA_MATRIX <- merge(CLINICAL_DATA_MATRIX, STRATA_FILE[, c("publc_id", STRATA_VARIABLE)], by.x= "Patient_ID", by.y = "public_id")
        }else if(!is.null(STRATA_PT_ID_COL_HEADER)){
          CLINICAL_DATA_MATRIX <- merge(CLINICAL_DATA_MATRIX, STRATA_FILE[, c (STRATA_PT_ID_COL_HEADER, STRATA_VARIABLE)], by.x= "Patient_ID", by.y = STRATA_PT_ID_COL_HEADER)
        }
        # MERGE CLINICAL FILE TABLE WITH THE STRATA FILE
        }else if ( COMMPASS_SURVIVAL_ANALYSIS == F){
          CLINICAL_DATA_MATRIX <- merge(CLINICAL_DATA_MATRIX, STRATA_FILE[, c( STRATA_PT_ID_COL_HEADER, STRATA_VARIABLE)], by.x = "Patient_ID", by.y = STRATA_PT_ID_COL_HEADER)
        }
        # ALLOWS USER TO ENTER AN EXISTING DATA FRAME AND STRATIFY FROM THAT AFTER MERGING WITH CLINICAL FILE.
        }else if(STRATIFY_FROM_EXISTING_DATA_FRAME == T){
          YOUR_DATA_FRAME <- subset(YOUR_DATA_FRAME, select = c(STRATA_VARIABLE, STRATA_PT_ID_COL_HEADER))
          CLINICAL_DATA_MATRIX <- merge(CLINICAL_DATA_MATRIX, YOUR_DATA_FRAME, by.x = "Patient_ID", by.y = STRATA_PT_ID_COL_HEADER)
        }
    # IF NOT STRAFIYING BUT WANT THE SURVIAL LINES
    # FUNCTION ADDS A COLUMN OF ALL 1'S 
    }else if (NEW_STRATA == F){
      CLINICAL_DATA_MATRIX$STRATA <- rep( 1, nrow( CLINICAL_DATA_MATRIX ))
    }
    
 
# SUBSET SECTION
  if (SUBSET == T){
    idx <- which(colnames(CLINICAL_DATA_MATRIX) == SUBSET_HEADER)
    CLINICAL_DATA_MATRIX <- CLINICAL_DATA_MATRIX[CLINICAL_DATA_MATRIX[,idx] == SUBSET_VALUE,]
  }
  

  # rename STRATA_VARIABLE column to pass into the survfit function
  colnames(CLINICAL_DATA_MATRIX)[colnames(CLINICAL_DATA_MATRIX) == STRATA_VARIABLE] <- "STRATA"
  
  
  #survfit function
  CLINICAL_DATA_MATRIX$OS <- with(CLINICAL_DATA_MATRIX, Surv(Time_To_Censored_OS, OS_Censor_Flag == 1))
  if (COMMPASS_SURVIVAL_ANALYSIS==T |(RENAME_PFS_CENSOR_HEADER != c("") | RENAME_PFS_TIME_HEADER !=c(""))){
  CLINICAL_DATA_MATRIX$PFS <- with(CLINICAL_DATA_MATRIX, Surv(Time_To_Censored_PFS, PFS_Censor_Flag == 1))
  }
  
  
  # create the survival objects used to plot kaplan-meyer curves
  OS <- survfit(OS ~ STRATA, data = CLINICAL_DATA_MATRIX, conf.type = "log-log")
  if (COMMPASS_SURVIVAL_ANALYSIS== T | (RENAME_PFS_CENSOR_HEADER != c("") | RENAME_PFS_TIME_HEADER !=c(""))){
  PFS <- survfit(PFS ~ STRATA, data = CLINICAL_DATA_MATRIX, conf.type = "log-log")
  }
  
  
  # calculate the number of groups
  NUMBER_OF_GROUPS <-  length(levels(as.factor(CLINICAL_DATA_MATRIX$STRATA)))
  
  
  # determine graphing logistics based on number of groups in strata category
  colors <- c("#5ab4ac","#d8b365", "#01665e", "#8c510a","#35978f","#c3a5cf","#543005","#c7eae5", "#762a83", "#5aae61", "#bf812d", "#00441b", "#9970ab", "#1b7837", "#4b0f55")
  lab_number <- 1:NUMBER_OF_GROUPS
  color <- colors[lab_number]
  
  
  # set PFS_OS User variable
  # set the survival analysis type between overall survival and progression-free survival
  # also set the plot title to reflect the survival type
  if (DEFAULT_TITLE == T){
    if (PFS_OS == c("PFS")){
      PFS_OS_user <- "PFS"
      PFS_OS <- PFS
      TITLE <- c("Progression-free Survival")
    } else if (PFS_OS == c("OS")) {
      PFS_OS_user <- "OS"
      PFS_OS <- OS
      TITLE <- c("Overall Survival")
    }
  }else if(DEFAULT_TITLE == F){
    TITLE <- CUSTOM_TITLE
    if (PFS_OS == c("PFS")){
      PFS_OS_user <- "PFS"
      PFS_OS <- PFS
    } else if (PFS_OS == c("OS")) {
      PFS_OS_user <- "OS"
      PFS_OS <- OS
    }
  }
  
  
  # adds or removes the median survival line
  if (MEDIAN_SURVIVAL_LINE == T){
    line_type <- c("hv")
  } else if (MEDIAN_SURVIVAL_LINE == F) {
    line_type <- c("none")
  }
  
  
  # hard codes PFS tables to false if they are not entered in
  if ( COMMPASS_SURVIVAL_ANALYSIS == F &&( RENAME_PFS_CENSOR_HEADER == c("") || RENAME_PFS_TIME_HEADER == c("") ) ) {
    PFS_COXPH_TABLE <-  F
    PFS_TABLE <-  F
  }
  

  # forces p-value to only show if there are two survival lines and no custom annotations
  # allows the user to add custom plot annotations for HR, median times, etc.
  if (PLOT_ANNOTATIONS_CALL == T){
    p_value <- F
    plot_text <- PLOT_ANNOTATIONS_TEXT
  } else if (PLOT_ANNOTATIONS_CALL == F && NUMBER_OF_GROUPS == 2 && P_VALUE == T) {
    p_value <- T
    plot_text <- c("")
  }else{ 
    plot_text <- c("")
    p_value <- P_VALUE
  }
  
  
  # create the actual survival plot 
  if (DEFAULT_LEGEND_LABELS == T && NUMBER_OF_GROUPS != 1){
    p_survival <- ggsurvplot(PFS_OS,
                             data = CLINICAL_DATA_MATRIX,
                             log = PFS_OS,
                             log.rank.weights = c("survdiff"),
                             pval = p_value,
                             pval.method.size = 3,
                             conf.int = CONFIDENCE_INTERVAL,
                             censor = CENSOR_TICK,
                             surv.median.line = line_type,
                             risk.table = RISK_TABLE,
                             risk.table.title = "",
                             risk.table.fontsize = RISK_TABLE_FONT,
                             risk.table.height = RISK_TABLE_HEIGHT,
                             risk.table.y.text = T,
                             break.time.by = 250,
                             risk.table.pos = c("out"),
                             palette = color,
                             title = TITLE,
                             xlab = X_LABEL,
                             ylim = c(0, 1.0),
                             ylab = Y_LABEL,
                             font.main = c(30, "plain", "black"),
                             pval.size = 5,
                             font.x = c(20, "plain", "black"),
                             font.y = c(20, "plain", "black"),
                             font.legend = c(15, "plain"),
                             font.tickslab = c(15, "plain", "black"),
                             legend.labs = as.character(lab_number),
                             legend.title = "",
                             ggtheme = theme(plot.title = element_text(hjust = 0.5),legend.justification="center"))
    
    p_survival$plot <- p_survival$plot + ggplot2::annotate("text", x = PLOT_ANNOTATIONS_X, y = PLOT_ANNOTATIONS_Y, label = plot_text, size = 5)
    p_survival
  
  } else if (DEFAULT_LEGEND_LABELS == F && NUMBER_OF_GROUPS != 1){
    p_survival <- ggsurvplot(PFS_OS,
                             data = CLINICAL_DATA_MATRIX,
                             log = (PFS_OS),
                             log.rank.weights = c("survdiff"),
                             pval = p_value,
                             pval.method.size = 3,
                             conf.int = CONFIDENCE_INTERVAL,
                             censor = CENSOR_TICK,
                             surv.median.line = line_type,
                             risk.table = RISK_TABLE,
                             risk.table.title = "",
                             risk.table.fontsize = RISK_TABLE_FONT,
                             risk.table.height = RISK_TABLE_HEIGHT,
                             risk.table.y.text = T,
                             break.time.by = 250,
                             risk.table.pos = c("out"),
                             palette = color,
                             title = TITLE,
                             xlab = X_LABEL,
                             ylim = c(0, 1.0),
                             ylab = Y_LABEL,
                             font.main = c(25, "plain", "black"),
                             pval.size = 5,
                             font.x = c(20, "plain", "black"),
                             font.y = c(20, "plain", "black"),
                             font.legend = c(15, "plain"),
                             font.tickslab = c(15, "plain", "black"),
                             legend.labs = CUSTOM_LEGEND_LABELS,
                             legend.title = "",
                             ggtheme = theme(plot.title = element_text(hjust = 0.5),legend.justification="center"))
    
    p_survival$plot <- p_survival$plot + ggplot2::annotate("text", x = PLOT_ANNOTATIONS_X, y = PLOT_ANNOTATIONS_Y, label = plot_text, size = 5) + guides(colour = guide_legend(nrow = 1))
    
  } else if (NUMBER_OF_GROUPS == 1){
    p_survival <- ggsurvplot(PFS_OS,
                             data = CLINICAL_DATA_MATRIX,
                             log = (PFS_OS),
                             log.rank.weights = c("survdiff"),
                             pval = p_value,
                             pval.method.size = 3,
                             conf.int = CONFIDENCE_INTERVAL,
                             censor = CENSOR_TICK,
                             surv.median.line = line_type,
                             risk.table = RISK_TABLE,
                             risk.table.title = "",
                             risk.table.fontsize = RISK_TABLE_FONT,
                             risk.table.height = RISK_TABLE_HEIGHT,
                             risk.table.y.text = F,
                             break.time.by = 250,
                             risk.table.pos = c("out"),
                             palette = color,
                             title = TITLE,
                             xlab = X_LABEL,
                             ylim = c(0, 1.0),
                             ylab = Y_LABEL,
                             font.main = c(25, "plain", "black"),
                             pval.size = 5,
                             font.x = c(20, "plain", "black"),
                             font.y = c(20, "plain", "black"),
                             font.legend = c(15, "plain"),
                             font.tickslab = c(15, "plain", "black"),
                             legend.labs = lab_number[-1],
                             legend.title = "",
                             ggtheme = theme(plot.title = element_text(hjust = 0.5)))
    
    p_survival$plot <- p_survival$plot + ggplot2::annotate("text", x = PLOT_ANNOTATIONS_X, y = PLOT_ANNOTATIONS_Y, label = plot_text, size = PLOT_ANNOTATIONS_SIZE) + guides(colour = guide_legend(nrow = 1))
  }

  
# OS Table
  if(OS_TABLE == T && NUMBER_OF_GROUPS != 1){
    CLINICAL_DATA_MATRIX$OS <- with(CLINICAL_DATA_MATRIX, Surv(Time_To_Censored_OS, OS_Censor_Flag == 1))
    OS <- survfit(OS ~ STRATA, data = CLINICAL_DATA_MATRIX, conf.type = "log-log")
    print("~~~~~~~~~~~~~~~~~~~ Overall Survival ~~~~~~~~~~~~~~~~~~~")
    print(OS)
  } else if (OS_TABLE == T | NUMBER_OF_GROUPS == 1){
    CLINICAL_DATA_MATRIX$OS <- with(CLINICAL_DATA_MATRIX, Surv(Time_To_Censored_OS, OS_Censor_Flag == 1))
    OS <- survfit(OS ~ 1, data = CLINICAL_DATA_MATRIX, conf.type = "log-log")
    print("~~~~~~~~~~~~~~~~~~~ Overall Survival ~~~~~~~~~~~~~~~~~~~")
    print(OS)
  }    else if (OS_TABLE == F | NUMBER_OF_GROUPS == 1){
    print("~~~~~~~~~~~~~~~~~~~ Overall Survival Not Calculated~~~~~~~~~~~~~~~~~~~")
  }

    
# PFS table
  if(PFS_TABLE == T && NUMBER_OF_GROUPS != 1){
    CLINICAL_DATA_MATRIX$PFS <- with(CLINICAL_DATA_MATRIX, Surv(Time_To_Censored_PFS, PFS_Censor_Flag == 1))
    PFS <- survfit(PFS ~ STRATA, data = CLINICAL_DATA_MATRIX, conf.type = "log-log")
    print("~~~~~~~~~~~~~~~~~~~ Progression-free Survival ~~~~~~~~~~~~~~~~~~~")
    print(PFS)
  }else if (PFS_TABLE == T | NUMBER_OF_GROUPS == 1){
    CLINICAL_DATA_MATRIX$PFS <- with(CLINICAL_DATA_MATRIX, Surv(Time_To_Censored_PFS, PFS_Censor_Flag == 1))
    PFS <- survfit(PFS ~ 1, data = CLINICAL_DATA_MATRIX, conf.type = "log-log")
    print("~~~~~~~~~~~~~~~~~~~ Progression-free Survival ~~~~~~~~~~~~~~~~~~~")
    print(PFS)
  }else if (PFS_TABLE == F | NUMBER_OF_GROUPS == 1){
    print("~~~~~~~~~~~ Progression-free Survival Not Calculated ~~~~~~~~~~~")
  }
  
  
# OS COXPH Table
  if(OS_COXPH_TABLE == T && NUMBER_OF_GROUPS != 1){
    OS_coxph <- summary(coxph(formula = Surv(CLINICAL_DATA_MATRIX$Time_To_Censored_OS, OS_Censor_Flag) ~ STRATA, data = CLINICAL_DATA_MATRIX))
    print("~~~~~~~~~~~~~~~~~~~ Cox proportional hazards regression model for Overall Survival ~~~~~~~~~~~~~~~~~~~")
    print(OS_coxph)
    OS_COX_PROPORTIONAL_HAZARD <- coxph(formula = Surv(CLINICAL_DATA_MATRIX$Time_To_Censored_OS, OS_Censor_Flag) ~ STRATA, data = CLINICAL_DATA_MATRIX)
  }else if (OS_COXPH_TABLE == F | NUMBER_OF_GROUPS == 1){
    print("~~~~~~~~~~~~~~~~~~~ Cox proportional hazards regression model for Overall Survival Not Calculated~~~~~~~~~~~~~~~~~~~")
  }
  
# PFS COXPH Table
  if(PFS_COXPH_TABLE == T && NUMBER_OF_GROUPS != 1){
    PFS_coxph <- summary(coxph(formula = Surv(CLINICAL_DATA_MATRIX$Time_To_Censored_PFS, PFS_Censor_Flag) ~ STRATA, data = CLINICAL_DATA_MATRIX))
    print("~~~~~~~~~~~~~~~~~~~ Cox proportional hazards regression model for Progression-free Survival ~~~~~~~~~~~~~~~~~~~")
    print(PFS_coxph)
    PFS_COX_PROPORTIONAL_HAZARD <- coxph(formula = Surv(CLINICAL_DATA_MATRIX$Time_To_Censored_PFS, PFS_Censor_Flag) ~ STRATA, data = CLINICAL_DATA_MATRIX)
  } else if (PFS_COXPH_TABLE == F | NUMBER_OF_GROUPS == 1){
    print("~~~~~~~~~~~ Cox proportional hazards regression model for Progression-free Survival Not Calculated ~~~~~~~~~~~")
  }
  
  COX_PROPORTIONAL_HAZARD <- coxph
  
  
  # pvalue matrix for pairwise comparison between strata variable
  if(PVAL_MATRIX==T){
    pval_matrix_OS <- pairwise_survdiff(OS ~ STRATA, data = CLINICAL_DATA_MATRIX, p.adjust.method = "none")
    print("~~~~~~~~~~~~~~~~~~~ Overall Survival p-values ~~~~~~~~~~~~~~~~~~~")
    print(pval_matrix_OS)
    if (COMMPASS_SURVIVAL_ANALYSIS==T | ( RENAME_PFS_CENSOR_HEADER != c("") | RENAME_PFS_TIME_HEADER !=c(""))){
    pval_matrix_PFS <- pairwise_survdiff(PFS ~ STRATA, data = CLINICAL_DATA_MATRIX, p.adjust.method = "none")
    print("~~~~~~~~~~~~~~~~~~~ Progression-free Survival p-values ~~~~~~~~~~~~~~~~~~~")
    print(pval_matrix_PFS)
    }
    
    if(PRINT_PVAL_MATRIX == T){
      # theme for pvalue table
      mytheme <- gridExtra::ttheme_default(
        core = list(fg_params=list(cex = PVAL_FONT_SIZE)),
        colhead = list(fg_params=list(cex = PVAL_HEADER_SIZE)),
        rowhead = list(fg_params=list(cex = PVAL_HEADER_SIZE)))
      
      if(PFS_OS_user == "PFS"){
        pval_df_PFS <- data.frame(pval_matrix_PFS$p.value)
        names(pval_df_PFS) <- substring(names(pval_df_PFS), 2)
        pval_df_PFS <- round(pval_df_PFS, 3)
        pval_plot <- tableGrob(pval_df_PFS, theme = mytheme)
        
      }else if(PFS_OS_user == "OS"){
        pval_df_OS <- data.frame(pval_matrix_OS$p.value)
        names(pval_df_OS) <- substring(names(pval_df_OS), 2)
        pval_df_OS <- round(pval_df_OS, 3)
        pval_plot <- tableGrob(pval_df_OS, theme = mytheme)
      }
    }
  
    #save pvalue matrix to Rdata
    #to reload: load("pval_matrix.RData")    
    if(SAVE_PVAL_MATRIX==T){
      save(pval_matrix, file = "pval_matrix.RData")
    }
  }
  
  # save plot
  if (SAVE_PLOT == T){
    png(filename = paste(FILE_NAME,".png", sep = ""), res = 300, width = 12, height = 10, units = "in")
    p_survival
    print(p_survival)
    dev.off()
  }
  
  
  # plot combined survival, risk table, and pvalue table
  if(RISK_TABLE == T & PRINT_PVAL_MATRIX == T){
    grid.arrange(p_survival$plot, p_survival$table, pval_plot) 
  }else if(RISK_TABLE == T & PRINT_PVAL_MATRIX == F){
    grid.arrange(p_survival$plot, p_survival$table) 
  }else if(RISK_TABLE == F & PRINT_PVAL_MATRIX == T){
    grid.arrange(p_survival$plot, pval_plot)
  }else{
    grid.arrange(p_survival$plot)
  }
  
  #print p_survial    
  p_survival   
  }
  

  
Survival_Function(COMMPASS_SURVIVAL_ANALYSIS = F, # must enter T or F
                  IA_RELEASE = 14, # enter CoMMpass release number if COMMPASS_SURVIVAL_ANALYSIS=T
                  
                  
                  YOUR_CLINICAL_FILE_PATH = c("/Volumes/MMRF/commpass/clinical_flat_files/IA14/CoMMpass_IA14_FlatFiles/MMRF_CoMMpass_IA14_PER_PATIENT.csv"), # must enter for non-CoMMpass 
                  RENAME_PATIENT_ID_HEADER = c("public_id"), # insert your patient id variable
                  RENAME_OS_TIME_HEADER = c("ttt2line"), # insert your overall survival time variable
                  RENAME_OS_CENSOR_HEADER = c("censt2line"), # insert your overall survival flag variable
                  RENAME_PFS_TIME_HEADER = c("ttt2linw"), # insert your progression free survival time variable
                  RENAME_PFS_CENSOR_HEADER = c("censet2line"), # insert your progression free survival flag variable
                  
                  # strata is how the survival function will be split up
                  # example: stratifying by sex would compare males to females
                  NEW_STRATA = T, # T or F
                  STRATA_VARIABLE = c("D_PT_race"), #insert the variable you want to stratify by 
                  STRATA_PT_ID_COL_HEADER= c("PUBLIC_ID"), # insert the column header name of our patient id variable. 
                  NEW_STRATA_FILE = T, # F is stratifying from stand alone survival(coMMpass), or clinical file path entered abolve (non-coMMpass)
                  STRATA_FILE_PATH = c("/Volumes/MMRF/commpass/clinical_flat_files/IA14/CoMMpass_IA14_FlatFiles/MMRF_CoMMpass_IA14_PER_PATIENT.csv"), # path you wish to stratify from  
                  STRATIFY_FROM_EXISTING_DATA_FRAME = F, # enter T or F
                  YOUR_DATA_FRAME = riss3, # insert name of data frame if stratifying from an existing data frame
                  
                  # subset is which group you want to analyze. 
                  # example: subsetting by sex would compare only females or only males
                  SUBSET = F,  # T or F
                  SUBSET_HEADER = c("sctflag"), # must be filled in if subset = T 
                  SUBSET_VALUE = 1, # insert a number if subset = T


                  PFS_OS= c("PFS"),   # Set the type of survival you want plotted: PFS = Progression-free survival, OS = Overall survival

                  
                  DEFAULT_LEGEND_LABELS = T, # T or F
                  CUSTOM_LEGEND_LABELS = c("Maintenance", "No Maintenance"), # insert custom legend label(s) if default legend labels = F
                  DEFAULT_TITLE = T, # T or F
                  CUSTOM_TITLE = c(""), # insert custom title if default title = F
                  X_LABEL = c("Days"), #  depends on the data, can be days, months, years, etc 
                  Y_LABEL = c("Survival Probability"), # specify Y axis label
                  
                  
                  RISK_TABLE = T, # T or F
                  RISK_TABLE_HEIGHT = 0.25, # DEFAULT = 0.25
                  RISK_TABLE_FONT = 5, # DEFAULT = 5
                  
                  
                  PVAL_MATRIX = T, # T or F
                  PRINT_PVAL_MATRIX = T, # T or F
                  SAVE_PVAL_MATRIX = F, # T or F
                  
                  
                  P_VALUE = T, # T or F
                  PVAL_HEADER_SIZE = 1.0, # DEFAULT = 1
                  PVAL_FONT_SIZE = 0.75, # DEFAULT = 0.75
                  
                  
                  CONFIDENCE_INTERVAL = T, # T or F
                  CENSOR_TICK = F, #  T or F
                  MEDIAN_SURVIVAL_LINE = T,  # T or F


                  PLOT_ANNOTATIONS_CALL = F, # T or F
                  PLOT_ANNOTATIONS_TEXT = c(""), # insert text to annotate
                  PLOT_ANNOTATIONS_X = 250,  # DEFAULT = 250
                  PLOT_ANNOTATIONS_Y = 0.15, # DEFAULT = .15
                  PLOT_ANNOTATIONS_SIZE = 5, # DEFAULT= 5
                  OS_TABLE = T, # T or F
                  PFS_TABLE = T, # T or F
                  OS_COXPH_TABLE = T,  # T or F
                  PFS_COXPH_TABLE = T,  # T or F
                  SAVE_PLOT  = T, # T or F
                  FILE_NAME = c("plot") ) # name the file if SAVE_PLOT= T 
 

  
