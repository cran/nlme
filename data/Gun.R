### $Id: Gun.R,v 1.1 2000/03/17 22:21:21 saikat Exp $
### The Gun data from the standard S version 3 distribution but as a groupedData
###   object.  Given as Table 7.8 in Hicks (1993), "Fundamental Concepts in the
###   Design of Experiments".
"Gun" <-
  structure(list(
   rounds = c(20.2, 14.2, 22, 14.1, 23.1, 14.1, 26.2, 
     18, 22.6, 14, 22.9, 12.2, 23.8, 12.5, 22.9, 13.7, 21.8, 12.7, 
     24.1, 16.2, 23.5, 16.1, 22.9, 16.1, 26.9, 19.1, 24.6, 18.1, 23.7, 
     13.8, 24.9, 15.4, 25, 16, 23.5, 15.1),
   Method = structure(factor(c(1, 
     2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 
     1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2), levels=1:2),
     class = "factor", .Label = c("M1", "M2")),
   Team = structure(ordered(c(1, 
     1, 4, 4, 7, 7, 3, 3, 5, 5, 9, 9, 2, 2, 6, 6, 8, 8, 1, 1, 4, 4, 
     7, 7, 3, 3, 5, 5, 9, 9, 2, 2, 6, 6, 8, 8), levels=1:9),
     class = c("ordered", "factor"),
     .Label = c("T1S", "T3S", "T2S", "T1A", "T2A", "T3A", "T1H", "T3H", "T2H")),
   Physique = structure(ordered(c(1, 
     1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 
     3, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3), levels=1:3),
     class = c("ordered", "factor"), 
     .Label = c("Slight", "Average", "Heavy"))),
formula = rounds ~ Method | Team, 
row.names = as.character(1:36),
class = c("nffGroupedData", "nfGroupedData", "groupedData", "data.frame"),
labels = list(y = "Number of rounds fired"), 
units = list(y = "(rounds/min)"), 
FUN = function (x) mean(x, na.rm = TRUE),
order.groups = TRUE)
