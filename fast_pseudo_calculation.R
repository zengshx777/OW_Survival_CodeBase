# By Dayne Batten
rmst_on_summary <- function(df, tmax) {
  not_dead <- mapply(function(eligible, event) {
    if (eligible == 0) {
      return(1)
    } else {
      return((eligible - event) / eligible)
    }
  },
  df$eligible[df$time < tmax],
  df$event[df$time < tmax])
  time_grid = c(df$time[df$time < tmax], tmax)
  time_diff = c(time_grid[1], diff(time_grid))
  return(sum(c(1, cumprod(not_dead)) * time_diff))
}

# Surv function
surv_on_summary <- function(df, tmax) {
  not_dead <- mapply(function(eligible, event) {
    if (eligible == 0) {
      return(1)
    } else {
      return((eligible - event) / eligible)
    }
  },
  df$eligible[df$time < tmax],
  df$event[df$time < tmax])
  return(prod(not_dead))
}


fast.pseudo <- function (time, event, tmax, type = "mean")
{
  if (any(is.na(time)))
    stop("Missing values in time vector.")
  if (any(time < 0))
    stop("Times must be nonnegative.")
  if (any(is.na(event)))
    stop("Missing values in event vector.")
  if (any(event != 0 & event != 1))
    stop("Event must be a binary variable (alive = 0, dead = 1).")
  if (missing(tmax) && (type == "surv"))
    stop("Need time point for survival function.")
  if (missing(tmax))
    tmax <- max(time[event == 1])
  if (length(tmax) > 1)
    stop("More than one value specified for tmax.")
  if (is.na(tmax))
    stop("Missing tmax value.")
  n <- length(time)
  df <- data.frame(id = 1:n, time, event)
  df$time[df$time >= tmax] <- tmax
  time_list <- sort(unique(time[event == 1]))
  
  summary <- do.call(rbind,
                     lapply(time_list,
                            function(x)
                              data.frame(
                                time = x,
                                eligible = sum(x <= df$time),
                                event = sum(df$event[df$time == x])
                              )))
  
  if (sum(summary$event[summary$time < tmax]) == 0)
    stop('No events occured before time tmax.')
  
  if (type == "mean") {
    main_res <- rmst_on_summary(summary, tmax)
  } else if (type == "surv") {
    main_res <- surv_on_summary(summary, tmax)
  } else{
    stop("Not supported type")
  }
  
  results <- do.call(rbind,
                     mapply(function(x, y) {
                       temp <- summary
                       temp$eligible[temp$time <= x] <-
                         temp$eligible[temp$time <= x] - 1
                       
                       if (y == 1) {
                         temp$event[temp$time == x] <- temp$event[temp$time == x] - 1
                       }
                       
                       if (type == "mean") {
                         this_res <- rmst_on_summary(temp, tmax)
                       } else if (type == "surv") {
                         this_res <- surv_on_summary(temp, tmax)
                       } else{
                         stop("Not supported type")
                       }
                       return(data.frame(
                         time = x,
                         event = y,
                         pseudo = (n * main_res) - ((n - 1) * this_res)
                       ))
                     },
                     df$time, df$event,
                     SIMPLIFY = FALSE))
  return(results$pseudo)
}