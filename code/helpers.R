#Function to aggregate data to user-specified timestep
#arguments:
##data: data frame
##dcol: character, name of data column containing datetimes
##vcol: character, name of data column containing data to be aggregated
##sd: datetime object (e.g., POSIXct) or character string interpretable as datetime, start datetime
##ed: datetime object (e.g., POSIXct) or character string interpretable as datetime, end datetime
##by: numeric, time bin size, in seconds
##func: character, aggregation function name (e.g., "mean", "sum", "min")
#Value: data frame with start and end datetimes for each bin, aggregated data, and count of observations

tagg = function(data, dcol, vcol, sd = min(data[,dcol], na.rm = TRUE), 
                ed = max(data[,dcol], na.rm = TRUE), by = 12*3600, 
                func = "mean"){
  
  if(length(na.omit(data[,dcol])) != nrow(data)){
    data = data[!is.na(data[,dcol]),]
    warning("Missing datetime values")
  }
  sd = as.POSIXct(sd)
  ed = as.POSIXct(ed)
  dts = seq(sd, ed, by = by)
  
  aggd = data.frame("sd" = dts, 
                    "ed" = dts + by - 1,
                    "val" = numeric(length(dts)), "nobs" = numeric(length(dts)))
  
  for(i in 1:length(dts)){
    data.sub = data[data[,dcol] >= sd & data[,dcol] <= ed,]
    aggd$val[i] = do.call(func, list(data.sub[,vcol], na.rm = TRUE))
    aggd$nobs[i] = length(na.omit(data.sub[,vcol]))
  }
  
  return(aggd)
}
