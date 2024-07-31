

#' Function for creating a input list for emmeans cov.reduce argument
make_cov_reduce_list <- function(fixedterms,dat){
  cov_reduce_list <- list()
  for(termm in fixedterms){
    if(!is.factor(dat[[termm]])){
      cov_reduce_list[[termm]] <- function(x){ c(mean(x) + sd(x), mean(x))}   
    }
  }
  return(cov_reduce_list)
}


#' Function for converting nebula model to emmeans object
#' 
prepare_mod_emmeans <- function(mod, term, fitmethod, fixed_formula, dat, cov_reduce_list){
  ## nebula (neg binomial or poisson mixed models) 
  ## requires special output coercion to emmeans format 
  if(fitmethod=="nebula::nebula"){
    neb_summ <- nebula_summary(mod)
    bhat <- neb_summ$msumm$est
    names(bhat) <- neb_summ$msumm$term
    tt <- terms(fixed_formula)
    allvar <- as.character(attr(tt, "variables"))[-1]
    offsetv <- allvar[attr(tt,"offset")] 
    offsetuse <- 1
    if(length(offsetv) > 0){
      offsetcol <- all.vars(as.formula(paste0("~",offsetv)))  
      offsetuse <- log(mean(dat[[offsetcol]]))
    }
    
    mod <- 
      emmeans::qdrg(formula=fixed_formula
                    ,data=dat
                    ,vcov = neb_summ$neb_cov
                    ,coef=bhat
                    ,link='log'
                    ,offset=offsetuse
                    ,cov.reduce = cov_reduce_list
      )
    return(mod) 
  } else {
    return(mod) 
  }
}
