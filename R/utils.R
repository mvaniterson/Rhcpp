#convert seconds to minutes, hours and days
.sec2time <- function(x) {
    time <- if(x < 60)
        paste0(round(x, 0), "s.")
    else if(x < 60*60)
        paste0(x%/%60, "m:", .sec2time(x%%60))
    else if(x < 60*60*24)
        paste0(x%/%(60*60), "h:", .sec2time(x%%(60*60)))
    else if(x < 60*60*24*7)
        paste0(x%/%(60*60*24), "d:", .sec2time(x%%(60*60*24)))
    time
}
