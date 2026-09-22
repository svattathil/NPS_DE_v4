library(data.table)


rm(list = ls())



###### GOAL ######



##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######


### Files that exist
indir <- "2_Pipeline/02-Prepare_analysis_data/"
infiles <- sapply(cohorts, function(acohort) { paste0(indir, "phenos_cleaned_SVs_", acohort, ".txt") },
                  simplify = FALSE)


### Files to be created
outfiles <- list(out = paste0(indir, "tabulate_hall_del.txt"))


###### MAIN ######
### Read in data
dats <- sapply(infiles, fread, simplify = FALSE)


### Prepare data
tabs <- sapply(dats, function(x) {
   addmargins(table(x[, .(hall, del)], useNA = "ifany"))
    }, simplify = FALSE)


dts <- sapply(tabs, function(x) {
    dt <- data.table(as.data.frame.matrix(x), keep.rownames = "hall")
    setnames(dt, c("0", "1"), c("del control", "del case"))
    dt[hall == 0, hall := "hall control"]
    dt[hall == 1, hall := "hall case"]
    return(dt)
}, simplify = FALSE)




###### FINISH ######
write("hall v. del counts", file = outfiles$out)

for(acohort in cohorts) {
    write("", file = outfiles$out, append = TRUE)
    write(acohort, file = outfiles$out, append = TRUE)
    write.table(dts[[acohort]], file = outfiles$out, row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
}
