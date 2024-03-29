homstart_data <- function(dir = NULL, load = TRUE, tdir = NULL) {
  if (is.null(tdir)) {
 tdir <- tempfile()
 dir.create(tdir)
 on.exit(unlink(tdir))
  }
  tdir <- path.expand(tdir)
  stopifnot(file.exists(tdir))
  if (is.null(dir)) 
 dir <- tdir
  dir <- path.expand(dir)
  stopifnot(file.exists(dir))
  rain <- TRUE
  elevation <- TRUE
  tempmin <- tempmax <- TRUE
  stations <- NULL
  all_stations <- stations
  include_rain <- rain
  include_tempmin <- tempmin
  include_tempmax <- tempmax
  include_elevation <- elevation
  savelist <- c("homstart", "rain", "tempmin", "tempmax", "levels")[c(TRUE, include_rain, include_tempmin, include_tempmax, TRUE)]
  if (!file.exists(file.path(dir, "homstart0.rda"))) {
 stations <- "id;name;long;lat;elevation;nied;tx;tn\ns100000;Hieflau;14,75;47,6;779;0;1;0\ns104000;Muerzzuschlag;15,6858333333333;47,6030555555556;758;0;0;1\ns105100;Reichenau an der Rax;15,8369444444444;47,6997222222222;485;0;0;1\ns111000;Bregenz;9,75;47,5;436;1;0;0\ns111100;Feldkirch;9,6;47,2666666666667;439;0;1;1\ns112000;Schoppernau;10,0186111111111;47,3111111111111;835;0;0;0\ns113000;Schroecken;10,0852777777778;47,2630555555556;1260;0;1;1\ns114000;Holzgau;10,3491666666667;47,2625;1100;0;2;2\ns115000;Reutte;10,75;47,5;870;0;1;1\ns118000;Innsbruck airport;11,3552777777778;47,2588888888889;579;0;1;0\ns118010;Innsbruck university;11,385;47,2605555555556;577;0;1;1\ns119000;Jenbach;11,75;47,3833333333333;530;0;1;1\ns123200;Zell am See;12,7833333333333;47,3;751;0;2;1\ns131100;Seckau;14,7833333333333;47,2833333333333;874;0;0;1\ns133000;Bruck an der Mur;15,2666666666667;47,4166666666667;489;0;0;0\ns140000;Kollerschlag;13,84;48,6063888888889;725;2;0;1\ns144000;Landeck;10,5666666666667;47,1333333333333;818;1;1;1\ns148100;Patscherkofel;11,4622222222222;47,2094444444445;2247;1;2;2\ns150000;Mayrhofen;11,85;47,15;643;1;1;1\ns153100;Mooserboden;12,7166666666667;47,15;2036;2;0;2\ns154000;Rauris;13;47,2166666666667;941;2;1;1\ns154010;Sonnblick;12,9575;47,0541666666667;3105;2;2;2\ns155000;Bad Gastein;13,1333333333333;47,1166666666667;1089;2;2;2\ns157100;Tamsweg;13,81;47,1247222222222;1022;2;1;1\ns159100;Stolzalpe;14,2;47,1166666666667;1299;0;1;1\ns160000;Freistadt;14,5;48,5166666666667;548;2;1;1\ns161000;Zeltweg;14,7833333333333;47,2;670;0;0;0\ns163000;Lobming;15,1783333333333;47,0430555555556;414;0;1;1\ns164000;Graz airport;15,4477777777778;46,9947222222222;337;0;1;1\ns164020;Graz university;15,4477777777778;47,0797222222222;366;0;1;1\ns164210;Schoeckl;15,4663888888889;47,1986111111111;1436;1;2;2\ns165000;Gleisdorf;15,7108333333333;47,1133333333333;375;0;0;1\ns167100;Woerterberg;16,0983333333333;47,2272222222222;400;0;1;2\ns170000;Galtuer;10,1833333333333;46,9666666666667;1587;1;2;2\ns173000;Obergurgl;11,0272222222222;46,8675;1938;0;2;2\ns177000;St. Jakob im Def.;12,3544444444444;46,9172222222222;1388;0;1;0\ns179000;Lienz;12,7833333333333;46,8166666666667;659;0;1;1\ns181000;Kolbnitz;13,3122222222222;46,8730555555556;603;1;1;1\ns188000;Preitenegg;14,9166666666667;46,9333333333333;1060;0;1;1\ns191000;Stift Zwettl;15,2;48,6166666666667;505;0;1;0\ns192010;Bad Gleichenberg;15,9019444444444;46,8666666666667;300;0;0;1\ns198000;Reisach;13,1541666666667;46,6483333333333;646;0;1;1\ns200200;Villacher Alpe;13,6733333333333;46,6036111111111;2140;0;2;2\ns201000;Kanzelhoehe;13,9066666666667;46,6780555555556;1526;1;1;1\ns202100;Klagenfurt;14,3333333333333;46,65;447;0;2;1\ns204000;St. Michael ob Bleiburg;14,75;46,5833333333333;500;0;0;0\ns211000;Loibl ;14,25;46,4444444444444;1098;0;1;1\ns240000;Laa an der Thaya;16,3852777777778;48,7261111111111;185;0;1;1\ns241000;Oberleis;16,37;48,5594444444444;420;2;1;0\ns260000;Hohenau;16,9044444444444;48,6172222222222;155;2;1;1\ns290000;Reichersberg;13,37;48,3361111111111;350;0;0;1\ns341000;Pabneukirchen;14,8194444444444;48,3188888888889;595;0;1;1\ns380000;Krems;15,6166666666667;48,4166666666667;190;2;1;1\ns500010;Hoersching;14,1911111111111;48,2411111111111;298;0;1;0\ns501000;Kremsmuenster;14,1322222222222;48,0552777777778;383;0;1;1\ns560000;St. Poelten;15,6166666666667;48,2;272;2;1;1\ns580400;Wien - Mariabrunn;16,2333333333333;48,2;227;2;2;2\ns590100;Wien - Hohe Warte;16,3577777777778;48,25;198;0;1;1\ns599000;Schwechat;16,5708333333333;48,1108333333333;184;0;1;1\ns630000;Salzburg airport;13,0016666666667;47,8013888888889;430;1;1;1\ns651000;Mondsee;13,3688888888889;47,8477777777778;482;0;1;2\ns661000;Feuerkogel;13,7183333333333;47,8177777777778;1618;2;2;2\ns691000;Grossraming;14,5177777777778;47,8933333333333;379;0;1;1\ns720000;Mariazell;15,3166666666667;47,7666666666667;865;0;1;2\ns770000;Eisenstadt;16,5166666666667;47,85;184;0;1;1\ns900000;Retz;15,95;48,7666666666667;256;0;2;2\ns901000;Kufstein;12,1638888888889;47,5741666666667;492;0;0;0\ns961000;Bad Ischl;13,6316666666667;47,7166666666667;469;0;1;1\ns962000;Krippenstein;13,7;47,5166666666667;2050;1;2;2\ns964000;Bad Aussee;13,7827777777778;47,6111111111111;660;2;1;2\ns981000;Irdning - Gumpenstein;14,1;47,5;710;0;1;1"
 cat(stations, file = file.path(tdir, "stations.csv"))
 homstart <- read.csv2(file.path(tdir, "stations.csv"))
 levels <- rownames(homstart)
 unlink(file.path(tdir, "stations.csv"))
 if (!file.exists(file.path(tdir, "HOMSTART-Daten.zip"))) {
   url <- "http://www.zamg.ac.at/docs/forschung/klimatologie/HOMSTART-Daten.zip"
   cat("downloading homstart data.\n")
   download.file(url, destfile = file.path(tdir, "HOMSTART-Daten.zip"))
 }
 if (!file.exists(file.path(tdir, "DATEN"))) {
   url <- "http://www.zamg.ac.at/docs/forschung/klimatologie/HOMSTART-Daten.zip"
   bsname <- basename(url)
   if (!file.exists(file.path(tdir, bsname))) 
  download.file(url, file.path(tdir, bsname))
   unzip(file.path(tdir, bsname), exdir = tdir)
 }
 stations <- if (is.null(all_stations)) {
   homstart$id
 }
 else {
   if (all(all_stations %in% homstart$id)) 
  all_stations <- which(homstart$id %in% all_stations)
   if (all(all_stations %in% homstart$name)) 
  all_stations <- which(homstart$name %in% all_stations)
   homstart$id[all_stations]
 }
 readhomstart <- function(f, which = 2, prefix = "s") {
   x <- zoo::read.zoo(f, header = FALSE, skip = 3, format = "%Y%m%d", na.string = "NA", drop = FALSE)
   colnames(x) <- paste(prefix, strsplit(f, "_", fixed = TRUE)[[1]][2], sep = "")
   x
 }
 cat("Extracting data.\n")
 wd <- getwd()
 setwd(tdir)
 on.exit(setwd(wd), add = TRUE)
 stations1 <- stations2 <- stations3 <- NULL
 if (include_rain) {
   rain <- do.call("merge", lapply(Sys.glob("DATEN/Niederschlag/*.txt"), readhomstart))
   rain <- rain[, names(rain) %in% stations, drop = FALSE]
   stations1 <- names(rain)
 }
 if (include_tempmin) {
   tempmin <- do.call("merge", lapply(Sys.glob("DATEN/Tmin/*.txt"), readhomstart))
   tempmin <- tempmin[, names(tempmin) %in% stations, drop = FALSE]
   stations2 <- names(tempmin)
 }
 if (include_tempmax) {
   tempmax <- do.call("merge", lapply(Sys.glob("DATEN/Tmax/*.txt"), readhomstart))
   tempmax <- tempmax[, names(tempmax) %in% stations, drop = FALSE]
   stations2 <- names(tempmax)
 }
 stations <- if (!is.null(stations2)) 
   intersect(stations1, stations2)
 else stations1
 unlink(file.path(tdir, "DATEN"))
 homstart <- homstart[homstart$id %in% stations, ]
 save(list = savelist, file = file.path(tdir, "homstart0.rda"), compress = "bzip2")
 remove(list = savelist)
  }
  if (!file.exists(file.path(dir, "homstart.rda"))) {
 load(file.path(tdir, "homstart0.rda"))
 harmonic <- function(x, order = 2, freq = 365) {
   x <- x/freq
   order <- round(order)
   stopifnot(order <= freq/2)
   rval <- outer(2 * pi * x, 1:order)
   rval <- cbind(apply(rval, 2, cos), apply(rval, 2, sin))
   colnames(rval) <- if (order == 1) {
  c("cos", "sin")
   }
   else {
  c(paste("cos", 1:order, sep = ""), paste("sin", 1:order, sep = ""))
   }
   if ((2 * order) == freq) 
  rval <- rval[, -(2 * order)]
   return(rval)
 }
 yday365 <- function(x) {
   x <- as.POSIXlt(x)
   mdays <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
   cumsum(c(0, mdays))[1 + x$mon] + x$mday
 }
 tspp <- function(x, breaks = c(0.01, 0.99, 4.99), labels = c("none", "low", "medium", "high"), origin = 1970, order = 2) {
   x0 <- x
   x <- x[, 1]
   tx <- time(x)
   stopifnot(inherits(tx, "Date"))
   is_weekday <- function(z) {
  if (inherits(z, "zoo")) 
  z <- time(z)
  z <- as.POSIXlt(z)$wday
  z > 0 & z < 6
   }
   is_weekend <- function(z) {
  if (inherits(z, "zoo")) 
  z <- time(z)
  z <- as.POSIXlt(z)$wday
  z < 1 | z > 5
   }
   is_Feb29 <- function(z) {
  if (inherits(z, "zoo")) 
  z <- time(z)
  z <- as.POSIXlt(z)
  z$mon == 1 & z$mday == 29
   }
   x <- x[i <- !is_Feb29(tx)]
   tx <- time(x)
   rval <- data.frame(raw = as.numeric(x), cens = ifelse(zoo::coredata(x) >= 0, zoo::coredata(x), 0), bin = factor(zoo::coredata(x) > 0, levels = c(FALSE, TRUE), labels = c("no", "yes")), cat = cut(zoo::coredata(x), breaks = c(-Inf, breaks, Inf), labels = labels), trend = 1900 - origin + as.POSIXlt(tx)$year + (yday365(tx) - 1)/365, month = factor(as.POSIXlt(tx)$mon + 1, levels = 1:12, labels = month.abb), year = factor(1900 + as.POSIXlt(tx)$year), day = yday365(tx), weekend = factor(is_weekend(tx), 
  labels = c("no", "yes")))
   if (include_tempmin) 
  rval$tempmin <- as.numeric(x0[i, "tempmin"])
   if (include_tempmax) 
  rval$tempmax <- as.numeric(x0[i, "tempmax"])
   rval$harmon <- harmonic(rval$day)
   return(rval)
 }
 k <- 1
 K <- length(names(rain))
 cat("Create station number", 1, "out of", K, "stations.")
 for (i in names(rain)) {
   if (!(i %in% homstart$id)) 
  next
   if (k > 1) 
  cat("\rCreate station number", k, "out of", K, "stations.")
   tm <- rain[, i]
   if (include_tempmin) 
  tm <- cbind(tm, tempmin = if (i %in% names(tempmin)) 
  tempmin[, i]
  else rep(NA, length(rain[, i])))
   if (include_tempmax) 
  tm <- cbind(tm, tempmax = if (i %in% names(tempmax)) 
  tempmax[, i]
  else rep(NA, length(rain[, i])))
   dat <- tspp(tm)
   dat <- data.frame(raw = dat$raw, cens = dat$cens, bin = as.integer(dat$bin) - 1, cat = as.integer(dat$cat) - 1, trend = dat$trend, month = dat$month, year = dat$year, day = dat$day, lon = homstart$long[homstart$id == i], lat = homstart$lat[homstart$id == i], id = as.integer(rep(rownames(homstart)[homstart$id == i], length = nrow(dat))), as.data.frame(dat$harmon), weekend = as.integer(dat$weekend) - 1)
   if (include_elevation) 
  dat$elevation <- homstart$elevation[homstart$id == i]
   for (ch in c(".", ":", "_", "-")) names(dat) <- gsub(ch, "", names(dat), fixed = TRUE)
   append <- if (k < 2) 
  FALSE
   else TRUE
   write.table(dat, file = file.path(tdir, "homstart_rain.raw"), quote = FALSE, row.names = FALSE, col.names = !append, append = append)
   k <- k + 1
 }
 cat("\n")
 cat("Creating final object homstart.rda\n")
 homstart <- read.table(file.path(tdir, "homstart_rain.raw"), header = TRUE)
 homstart$id = factor(as.integer(homstart$id), levels = as.integer(levels), labels = levels)
 homstart$weekend = factor(as.integer(homstart$weekend), levels = c(0, 1), labels = c("no", "yes"))
 homstart$bin = factor(as.integer(homstart$bin), levels = c(0, 1), labels = c("no", "yes"))
 homstart$cat = factor(as.integer(homstart$cat), levels = c(0, 1, 2, 3), labels = c("none", "low", "medium", "high"))
 homstart$raw[homstart$raw < 0] <- 0
 save(homstart, file = file.path(dir, "homstart.rda"), compress = "bzip2")
 remove(list = savelist)
  }
  if (load) {
 if ("homstart" %in% ls(envir = .GlobalEnv)) 
   remove("homstart", envir = .GlobalEnv)
 load(file.path(dir, "homstart.rda"), envir = .GlobalEnv)
  }
}
