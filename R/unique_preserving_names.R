unique_preserving_names<-function(vector) tapply(vector, names(vector), FUN = head, 1)  #https://stackoverflow.com/questions/42714163/how-do-i-get-unique-element-from-a-vector-keeping-its-name

#this returns named vector of unique values (unique within each name only)
