### $Id: Machines.R,v 1.3 2000/07/03 18:22:44 bates Exp $
### Productivity scores of different workers on three different machines
### Milliken and Johnson (1992), section 23.1, p. 285.
"Machines" <-
  structure(list
  (Worker = structure(ordered(c(4, 4, 4, 2, 2, 2, 
     5, 5, 5, 3, 3, 3, 6, 6, 6, 1, 1, 1, 4, 4, 4, 2, 2, 2, 5, 5, 5, 
     3, 3, 3, 6, 6, 6, 1, 1, 1, 4, 4, 4, 2, 2, 2, 5, 5, 5, 3, 3, 3, 
     6, 6, 6, 1, 1, 1), levels=1:6), class = c("ordered", "factor"),
     .Label = c("6", "2", "4", "1", "3", "5")),
   Machine = structure(factor(c(1, 
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 
     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 
     3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3), levels=1:3), class = "factor",
     .Label = c("A", "B", "C")),
   score = c(52, 52.8, 
     53.1, 51.8, 52.8, 53.1, 60, 60.2, 58.4, 51.1, 52.3, 50.3, 50.9, 
     51.8, 51.4, 46.4, 44.8, 49.2, 62.1, 62.6, 64, 59.7, 60, 59, 68.6, 
     65.8, 69.7, 63.2, 62.8, 62.2, 64.8, 65, 65.4, 43.7, 44.2, 43, 
     67.5, 67.2, 66.9, 61.5, 61.7, 62.3, 70.8, 70.6, 71, 64.1, 66.2, 
     64, 72.1, 72, 71.1, 62, 61.4, 60.5)),
class = c("nffGroupedData", "nfGroupedData", "groupedData", "data.frame"),
row.names = as.character(1:54),
formula = score ~ Machine | Worker,
labels = list(y = "Productivity score"),
inner = ~ Machine)
