
fn <- function(p, mod, data, parname, yobs,pred_col,pred=FALSE) {
  p <- lapply(p,exp) %>% rlang::set_names(parname)
  mod <- param(mod,p)
  out <- mrgsim_q(mod,data,output="df")
  if(pred) return(out)
  sum((yobs - out[[pred_col]])^2,na.rm=TRUE)
}

fit_newuoa <- function(mod, dv_col, pred_col, theta,...) {
  fit_nloptr(mod,dv_col,pred_col,theta,...)
}


fit_nloptr <- function(mod, dv_col, pred_col, theta, algorithm = "NLOPT_LN_NEWUOA", 
                       xtol_rel = 1E-4) {
  opts <- list(algorithm = algorithm, xtol_rel = xtol_rel)
  fit <- nloptr(theta, fn, mod = mod, data = data, 
                       parname = names(theta), yobs = data[[dv_col]], 
                       pred_col = pred_col, opts = opts, pred = FALSE)
  data[["PRED"]] <- fn(fit$solution,mod,data,names(theta),yobs,pred_col,pred=TRUE)[[pred_col]]
  list(par = fit$solution, full_data = data, fit = fit, data = filter(data, evid==0))
}


