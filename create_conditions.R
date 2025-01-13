library(dplyr)
creat_conditions<- function(data){
c_names<- colnames(data)%>% .[. !="Geneid"]
conditions<- substr(c_names, 1, (nchar(c_names)-2) )
condition_list <- split(seq_along(conditions), conditions)
return(condition_list)
}
