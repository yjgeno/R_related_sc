read.table(): a general R base function that can be used to read a file in table format. The data will be imported as a data frame.
read.csv(): for reading “comma separated value” files (“.csv”).
read.csv2(): variant used in countries that use a comma “,” as decimal point and a semicolon “;” as field separators.
read.delim(): for reading “tab-separated value” files (“.txt”). By default, point (“.”) is used as decimal points.
read.delim2(): for reading “tab-separated value” files (“.txt”). By default, comma (“,”) is used as decimal points.

examples:
read.csv(file, header = TRUE, sep = ",", dec = ".", ...)
read.delim(file, header = TRUE, sep = "\t", dec = ".", ...)
