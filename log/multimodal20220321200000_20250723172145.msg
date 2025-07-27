Error in validate_rounding_nunit(.Call(C_parse_unit, as.character(unit))) : Rounding with minute > 60 is not supported. Use aseconds for arbitrary units. 
Traceback: 
[1] "5: stop(sprintf(\"Rounding with %s > %d is not supported. Use aseconds for arbitrary units.\", " "       nunit$unit, trunc_multi_limits[[nunit$unit]]))"                                          
[3] "4: validate_rounding_nunit(.Call(C_parse_unit, as.character(unit)))"                             "3: timechange:::parse_units(unit)"                                                              
[5] "2: .normalize_multi_week_unit(unit)"                                                             "1: round_date(data[[time.index]], paste(frequency, \" mins\")) at #123"                         

