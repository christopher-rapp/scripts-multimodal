Error in FUN(X[[i]], ...) : object 'max.iterations' not found 
Traceback: 
 [1] "2: FUN(X[[i]], ...)"                                                    "1: lapply(data.ls, function(z) {"                                      
 [3] "       date.time.c <- lubridate::as_datetime(as.numeric(first(z[, "     "           time.index, with = F])))"                                   
 [5] "       {"                                                               "           bins.ix = str_which(colnames(z), \"\\\\d+(\\\\.\\\\d+)?$\")"
 [7] "           bins.nm = as.numeric(colnames(z)[bins.ix])"                  "           tmp.df <- z %>% select(!all_of(bins.ix))"                   
 [9] "           {"                                                           "               tmp1 <- as.numeric(bins.nm[-c(length(bins.nm))])/1000"  
[11] "               tmp2 <- as.numeric(bins.nm[-1])/1000"                    "               dlogDp = log10(tmp2) - log10(tmp1)"                     
[13] "           }"                                                           "           select.c <- bins.ix[1:(length(bins.ix) - 1)]"               
[15] "       }"                                                               "       {"                                                              
[17] "           dNdlogDp.mat <- z[, select.c, with = F]"                     "           dNdlogDp.mat <- sapply(dNdlogDp.mat, as.numeric)"           
[19] "           dN.mat = t(t(dNdlogDp.mat) * dlogDp)"                        "           dN.mat <- dN.mat[complete.cases(dN.mat), ]"                 
[21] "           if (is.null(nrow(dNdlogDp.mat))) {"                          "     ... at #137"                                                      

