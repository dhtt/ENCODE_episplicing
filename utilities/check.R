library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(stringr)
library(stats)
library(NMF)

setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/file_transfer/res")
folder = list.files("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/file_transfer/res", 
                    full.names = FALSE, pattern = ".csv")
get_all_pairs.exp <- function(all_pairs.exp, thres_p = 0.01){
  pair.exp_list = vector("list", length(all_pairs.exp))
  for (i in 1:length(all_pairs.exp)){
    print(paste("Pair: ",i, sep=''))
    pair.exp = all_pairs.exp[[i]]
    pair.exp_list[[i]] = fread(pair.exp, header=TRUE, stringsAsFactors = FALSE, quote=FALSE)
    print(head(pair.exp_list[[i]]))
  }
  all_res_pair = vector("list", length(all_pairs.exp))
  for (i in 1:length(pair.exp_list)){
    print(paste("Pair: ", i, sep=''))
    exp_table = pair.exp_list[[i]]
    exp_res = exp_table %>%
      filter(padj <= thres_p) %>%
      group_by(groupID) %>%
      mutate(DEU = paste(featureID, collapse = ', ')) %>%
      dplyr::select(groupID, DEU) %>%
      unique()
    all_res_pair[[i]] = exp_res
  }
  names(all_res_pair) = sapply(folder, function(x) str_split(x, ".csv")[[1]][1] )
  return(all_res_pair)
}
all_DEU_res1 = get_all_pairs.exp(folder, 0.01)
all_DEU_res2 = get_all_pairs.exp(folder, 0.05)
all_DEU_res3 = get_all_pairs.exp(folder, 0.1)

length(intersect(all_DEU_res3$EH1_Esmallintestine_res$groupID, all_DEU_res3$RH1_Rsmallintestine_res$groupID))/
  length(union(all_DEU_res3$EH1_Esmallintestine_res$groupID, all_DEU_res3$RH1_Rsmallintestine_res$groupID)) * 100

length(intersect(all_DEU_res3$EH1_Espleen_res$groupID, all_DEU_res3$RH1_Rspleen_res$groupID))/
  length(union(all_DEU_res3$EH1_Espleen_res$groupID, all_DEU_res3$RH1_Rspleen_res$groupID)) * 100

length(intersect(all_DEU_res3$Espleen_Esmallintestine_res$groupID, all_DEU_res3$Rspleen_Rsmallintestine_res$groupID))/
  length(union(all_DEU_res3$Espleen_Esmallintestine_res$groupID, all_DEU_res3$Rspleen_Rsmallintestine_res$groupID))*100


all_DEU_res_en = get_all_pairs.exp(folder, 0.001)
all_DEU_res = get_all_pairs.exp(folder, 0.5)
length(all_DEU_res)
lapply(all_DEU_res, length)
temp = as.data.frame(lapply(all_DEU_res3, function(x) lapply(x, length)))
temp = temp[, seq(from = 1, to = length(temp), by = 2)]
temp

length(intersect(all_DEU_res$EH1_Esmallintestine_res$groupID, all_DEU_res_en$RH1_Rsmallintestine_res$groupID))/
  length(union(all_DEU_res$EH1_Esmallintestine_res$groupID, all_DEU_res_en$RH1_Rsmallintestine_res$groupID)) * 100

length(intersect(all_DEU_res$EH1_Espleen_res$groupID, all_DEU_res_en$RH1_Rspleen_res$groupID))/
  length(union(all_DEU_res$EH1_Espleen_res$groupID, all_DEU_res_en$RH1_Rspleen_res$groupID)) * 100

length(intersect(all_DEU_res$Espleen_Esmallintestine_res$groupID, all_DEU_res_en$Rspleen_Rsmallintestine_res$groupID))/
  length(union(all_DEU_res$Espleen_Esmallintestine_res$groupID, all_DEU_res_en$Rspleen_Rsmallintestine_res$groupID))*100


#==============================
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/file_transfer/res/other")
folder = list.files("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/file_transfer/res/other", 
                    full.names = FALSE, pattern = ".csv")
all_DEU_res2 = get_all_pairs.exp(folder)
length(all_DEU_res2)
lapply(all_DEU_res2, length)

temp2 = as.data.frame(lapply(all_DEU_res2, function(x) lapply(x, length)))
temp2 = temp2[, seq(from = 1, to = length(temp2), by = 2)]

