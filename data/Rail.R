### $Id: Rail.R,v 1.1 2000/03/17 22:21:21 saikat Exp $
### Travel times of ultrasonic head waves in several different railway
###  rails. 
### Devore (1995), example 10.10
"Rail" <-
  structure(list(
  Rail = structure(ordered(c(3, 3, 3, 1, 1, 1, 5, 
    5, 5, 6, 6, 6, 2, 2, 2, 4, 4, 4), levels=1:6),
    class = c("ordered", "factor"),
    .Label = c("2", "5", "1", "6", "3", "4")),
travel = c(55, 53, 54, 26, 37, 32, 78, 91, 85, 92, 100, 96, 49,
  51, 50, 80, 85, 83)),
class = c("nffGroupedData", "nfGroupedData", "groupedData", "data.frame"),
row.names = as.character(1:18),
labels = list(y = "Zero-force travel time"),
units = list(y = "(nanoseconds)"),
formula = travel ~ 1 | Rail,
order.groups = TRUE)
