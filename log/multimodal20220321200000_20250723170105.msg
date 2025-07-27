Error in FUN(X[[i]], ...) : object 'lower.limit' not found 
Traceback: 
 [1] "4: seq(lower.limit, upper.limit, by = 0.01) at #438"                    "3: data.frame(x = seq(lower.limit, upper.limit, by = 0.01)) at #438"   
 [3] "2: FUN(X[[i]], ...)"                                                    "1: lapply(data.ls, function(z) {"                                      
 [5] "       date.time.c <- lubridate::as_datetime(as.numeric(first(z[, "     "           time.index, with = F])))"                                   
 [7] "       {"                                                               "           bins.ix = str_which(colnames(z), \"\\\\d+(\\\\.\\\\d+)?$\")"
 [9] "           bins.nm = as.numeric(colnames(z)[bins.ix])"                  "           tmp.df <- z %>% select(!all_of(bins.ix))"                   
[11] "           {"                                                           "               tmp1 <- as.numeric(bins.nm[-c(length(bins.nm))])/1000"  
[13] "               tmp2 <- as.numeric(bins.nm[-1])/1000"                    "               dlogDp = log10(tmp2) - log10(tmp1)"                     
[15] "           }"                                                           "           select.c <- bins.ix[1:(length(bins.ix) - 1)]"               
[17] "       }"                                                               "       {"                                                              
[19] "           dNdlogDp.mat <- z[, select.c, with = F]"                     "           dNdlogDp.mat <- sapply(dNdlogDp.mat, as.numeric)"           
[21] "           dN.mat = t(t(dNdlogDp.mat) * dlogDp)"                        "           dN.mat <- dN.mat[complete.cases(dN.mat), ]"                 
[23] "           if (is.null(nrow(dNdlogDp.mat))) {"                          "     ... at #137"                                                      

