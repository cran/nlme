### $Id: ergoStool.R,v 1.1 2000/03/17 22:21:21 saikat Exp $
### Effort required to arise from a stool. 9 different subjects and 4 different
###  types of stools.
### Devore (1995), exercise 11.9
"ergoStool" <-
  structure(list(
effort = c(12, 15, 12, 10, 10, 14, 13, 12, 7, 
    14, 13, 9, 7, 11, 10, 9, 8, 11, 8, 7, 9, 11, 11, 10, 8, 12, 12, 
    11, 7, 11, 8, 7, 9, 13, 10, 8),
Type = structure(c(1, 2, 3, 4, 
  1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 
  2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
  .Label = c("T1", "T2", "T3", "T4"), class = "factor"),
Subject = structure(c(8, 8, 8, 8, 9, 9, 
  9, 9, 6, 6, 6, 6, 3, 3, 3, 3, 2, 2, 2, 2, 5, 5, 5, 5, 7, 7, 7, 
  7, 1, 1, 1, 1, 4, 4, 4, 4),
  .Label = c("8", "5", "4", "9", "6", "3", "7", "1", "2"),
  class = c("ordered", "factor"))),
row.names = as.character(1:36),
class = c("nffGroupedData", "nfGroupedData", "groupedData", "data.frame"),
formula = effort ~ Type | Subject,
labels = list(x = "Type of stool", y = "Effort required to arise"),
units = list(y = "(Borg scale)"),
FUN = function (x) mean(x, na.rm = TRUE),
order.groups = TRUE)
