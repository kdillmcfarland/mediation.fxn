"
Mediation analysis and plotting

#################

Kim Dill-McFarland & Max Segnitz
University of Washington
Copyright (C) 2020

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
    This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. <http://www.gnu.org/licenses/>

Runs mediation analysis of 1 mediator on 1 independent and dependent variable

REQUIRED
  dat = data frame containing all variables in model
  dv.family = character string specifying family for models with outcome DV. 
              Default is 'gaussian'. See ?family for options
  mediator.family = character string specifying family for models with
                    outcome mediator. Default is 'gaussian'. 
                    See ?family for options

If modeling all variables against each other:
  iv = character string of independent variable(s) in dat. Must specify
       interaction terms as 'a:b'
  dv = character string of dependent variable(s) in dat
  mediator = character string of mediator variable(s) in dat
  
If modeling specific pairs within variables
(interaction terms not supported)
  var.dat = data frame with columns 'iv', 'dv', 'mediator' specifying
            character strings for variables to run in models. Each row is
            assessed together
  
OPTIONAL
  intercept = logical if should include an intercept in the model. Default is TRUE
  coVar = character string of co-variates to include in model
  randVar = character string of random effects to include in model
  boot = logical. If FALSE, quasi-Bayesian approximation  used for confidence
         intervals. If TRUE, nonparametric bootstrap used. Default is FALSE
  plot = logical if should output a flowchart plot of the mediation results.
         Default is FALSE
  outdir = character string for output directory. Default is working directory
  prefix = character string to name outputs saved to disk. Default is 'mediation'
"

mediation.fxn <- function(dat, iv=NULL, dv=NULL, mediator=NULL,
                          var.dat=NULL,
                          intercept=TRUE,
                          dv.family = "gaussian", 
                          mediator.family = "gaussian",
                          coVar=NULL, randVar=NULL, plot=FALSE,
                          boot=FALSE, outdir=NULL, prefix="mediation", ...){
  #### Setup ####
  #load packages
  require(plyr, quietly = TRUE)
  require(tidyverse, quietly = TRUE)
  require(broom, quietly = TRUE) 
  require(broom.mixed, quietly = TRUE)
  require(lme4, quietly = TRUE)
  library(car, quietly = TRUE)
  require(mediation, quietly = TRUE)
  #silence warnings
  options(warn=-1)
  
  #Blanks lists to hold results
  model.result <- list()
  mediation.result <- list()
  
  #Assign some variables to global environment for mediate
  assign("dv.family", dv.family, envir = .GlobalEnv)
  assign("mediator.family", mediator.family, envir = .GlobalEnv)
  
  #Force boot=TRUE if using binomial models
  if(dv.family=="binomial" & boot==FALSE){
    boot=TRUE
    message("Binomial DV models require nonparametric bootstrap. Setting boot to TRUE.")
  }
  
  #### Loops for iv, dv, mediator vectors ####
  if(is.null(var.dat)){
    for(DV in dv){
      for(MED in mediator){
        for(IV in iv){
          print(IV)
          #### Filter complete cases ####
          #Include covariates and random effects if supplied
          dat.sub <- dat %>% 
            dplyr::select(1, all_of(c(IV, DV, MED, coVar, randVar))) %>% 
            drop_na()
          
          #### Make model formulae ####
          # Note that reformulate( ) works to run models but these models fail in mediate( )
          # due to variable name issues. Thus, complex paste( ) statements must be used
          
          # WITHOUT random effects
          if(is.null(randVar)){
            m.totEffect <- paste(DV, "~", paste(c(IV, coVar), collapse=" + "), sep=" ")
            m.mediator <- paste(MED, "~", paste(c(IV, coVar), collapse=" + "), sep=" ")
            m.dv <- paste(DV, "~", paste(c(IV, MED, coVar), collapse=" + "), sep=" ")
          } else {
            # WITH random effects
            m.totEffect <- paste(DV, "~", paste(c(IV, coVar), collapse=" + "), 
                                 "+ (1|", randVar, ")", sep=" ")
            m.mediator <- paste(MED, "~", paste(c(IV, coVar), collapse=" + "), 
                                "+ (1|", randVar, ")", sep=" ")
            m.dv <- paste(DV, "~", paste(c(IV, MED, coVar), collapse=" + "), 
                          "+ (1|", randVar, ")", sep=" ")
          }
          
          # Remove intercept if necessary
          if(intercept == FALSE){
            m.totEffect <- gsub("~", "~ 0 +", m.totEffect)
            m.mediator <- gsub("~", "~ 0 +", m.mediator)
            m.dv <- gsub("~", "~ 0 +", m.dv)
          }
          
          #### Run models ####
          # WITHOUT random effects
          if(is.null(randVar)){ 
            # Total effect of independent var on dependent var
            fit.totEffect <- glm(m.totEffect, data=dat.sub, 
                                 family = dv.family)
            # Effect of independent var on mediator
            fit.mediator <- glm(m.mediator, data=dat.sub, 
                                family = mediator.family)
            # Effect mediator on dependent var
            fit.dv <- glm(m.dv, data=dat.sub, 
                          family = dv.family)
            
          } else{
            # WITH random effects
            # Total effect of independent var on dependent var
            fit.totEffect <- glmer(m.totEffect, data=dat.sub, 
                                   family = dv.family)
            # Effect of independent var on mediator
            fit.mediator <- glmer(m.mediator, data=dat.sub,
                                  family = mediator.family)
            # Effect mediator on dependent var
            fit.dv <- glmer(m.dv, data=dat.sub, 
                            family = dv.family)
          }
          
          #### MEDIATION ####
          #Assign models to global environment for mediate
          assign("m.mediator", m.mediator, envir = .GlobalEnv)
          assign("m.dv", m.dv, envir = .GlobalEnv)
        
          tryCatch({
          mediation <- mediate(model.m = fit.mediator,
                               model.y = fit.dv, 
                               treat = IV, 
                               mediator = MED,
                               covariates = coVar,
                               boot = boot)
          
          #### Extract results ####
          # WITHOUT random effects
          if(is.null(randVar)){
            result <- tidy(fit.totEffect) %>% 
              bind_rows(tidy(fit.mediator)) %>% 
              bind_rows(tidy(fit.dv)) %>% 
              mutate(model = c(rep("DV~IV", nrow(tidy(fit.totEffect))),
                               rep("MED~IV", nrow(tidy(fit.mediator))),
                               rep("DV~IV+MED", nrow(tidy(fit.dv))))) %>% 
              dplyr::select(model, everything())
          } else{
            # WITH random effects
            m.result <- tidy(fit.totEffect) %>% 
              bind_rows(tidy(fit.mediator)) %>% 
              bind_rows(tidy(fit.dv)) %>% 
              mutate(model = c(rep("DV~IV", nrow(tidy(fit.totEffect))),
                               rep("MED~IV", nrow(tidy(fit.mediator))),
                               rep("DV~IV+MED", nrow(tidy(fit.dv)))))
            #Estimate p-values
            p.result <- tidy(Anova(fit.totEffect)) %>% 
              bind_rows(tidy(Anova(fit.mediator))) %>% 
              bind_rows(tidy(Anova(fit.dv))) %>% 
              mutate(model = c(rep("DV~IV", nrow(tidy(Anova(fit.totEffect)))),
                               rep("MED~IV", nrow(tidy(Anova(fit.mediator)))),
                               rep("DV~IV+MED", nrow(tidy(Anova(fit.dv)))))) %>% 
              mutate(effect="fixed") %>% 
              rename(statistic.Anova = statistic)
            
            result <- full_join(m.result, p.result) %>% 
              dplyr::select(model, everything())
          }
          
          #Mediation analysis
          mediation.summ <- summary(mediation)
          
          result.mediation <- data.frame(
            #ACME = average causal mediation effects
            #ADE = average direct effects
            model="mediation",
            term = c("ACME","ADE","Total Effect","Prop. Mediated"),
            estimate = c(mediation.summ$d0, mediation.summ$z0,
                         mediation.summ$tau.coef, mediation.summ$n0),
            
            CI95_lower =c(mediation.summ$d0.ci[1], mediation.summ$z0.ci[1],
                          mediation.summ$tau.ci[1], mediation.summ$n0.ci[1]),
            
            CI95_upper = c(mediation.summ$d0.ci[2], mediation.summ$z0.ci[2],
                           mediation.summ$tau.ci[2], mediation.summ$n0.ci[2]),
            
            p.value = c(mediation.summ$d0.p, mediation.summ$z0.p,
                        mediation.summ$tau.p, mediation.summ$n0.p))
          
          #### Save results to list ####
          name <- paste(DV,IV,MED, sep=".")
          model.result[[name]] <- result
          mediation.result[[name]] <- result.mediation
          
          #### Plot ####
          if(plot){
            basic.plot(result=result, DV=DV, IV=IV, MED=MED,
                       result.mediation=result.mediation,
                       outdir=outdir, prefix=prefix, name=name)
          }
          }, error=function(e){
            print(paste("Could not complete mediation for ", 
                        DV, " ~ ", IV, " + ", MED))
          })
        }}}
  } else{
  #### Loop for var.dat pairs ####
    for(index in 1:nrow(var.dat)){
      IV <- var.dat$iv[index]
      DV <- var.dat$dv[index]
      MED <- var.dat$mediator[index]
      print(IV)
      
    #### Filter complete cases ####
    #Include covariates and random effects if supplied
    dat.sub <- dat %>% 
      dplyr::select(1, all_of(c(IV, DV, MED, coVar, randVar))) %>% 
      drop_na()
          
      assign("dat.sub",dat.sub, envir = .GlobalEnv)
    #### Make model formulae ####
    # Note that reformulate( ) works to run models but these models fail in mediate( )
    # due to variable name issues. Thus, complex paste( ) statements must be used
          
    # WITHOUT random effects
    if(is.null(randVar)){
      m.totEffect <- paste(DV, "~", paste(c(IV, coVar), collapse=" + "), sep=" ")
      m.mediator <- paste(MED, "~", paste(c(IV, coVar), collapse=" + "), sep=" ")
      m.dv <- paste(DV, "~", paste(c(IV, MED, coVar), collapse=" + "), sep=" ")
    } else {
      # WITH random effects
      m.totEffect <- paste(DV, "~", paste(c(IV, coVar), collapse=" + "), 
                                 "+ (1|", randVar, ")", sep=" ")
      m.mediator <- paste(MED, "~", paste(c(IV, coVar), collapse=" + "), 
                                "+ (1|", randVar, ")", sep=" ")
      m.dv <- paste(DV, "~", paste(c(IV, MED, coVar), collapse=" + "), 
                          "+ (1|", randVar, ")", sep=" ")
    }
          
    # Remove intercept if necessary
    if(intercept == FALSE){
      m.totEffect <- gsub("~", "~ 0 +", m.totEffect)
      m.mediator <- gsub("~", "~ 0 +", m.mediator)
      m.dv <- gsub("~", "~ 0 +", m.dv)
    }
          
    #### Run models ####
    # WITHOUT random effects
    if(is.null(randVar)){ 
      # Total effect of independent var on dependent var
      fit.totEffect <- glm(m.totEffect, data=dat.sub, 
                                 family = dv.family)
      # Effect of independent var on mediator
      fit.mediator <- glm(m.mediator, data=dat.sub, 
                                family = mediator.family)
      # Effect mediator on dependent var
      fit.dv <- glm(m.dv, data=dat.sub, 
                          family = dv.family)
            
    } else{
      # WITH random effects
      # Total effect of independent var on dependent var
      fit.totEffect <- glmer(m.totEffect, data=dat.sub, 
                                   family = dv.family)
      # Effect of independent var on mediator
      fit.mediator <- glmer(m.mediator, data=dat.sub,
                                  family = mediator.family)
      # Effect mediator on dependent var
      fit.dv <- glmer(m.dv, data=dat.sub, 
                            family = dv.family)
    }
          
    #### MEDIATION ####
    #Assign models to global environment for mediate
    assign("m.mediator", m.mediator, envir = .GlobalEnv)
    assign("m.dv", m.dv, envir = .GlobalEnv)
     
    tryCatch({     
    mediation <- mediate(model.m = fit.mediator,
                               model.y = fit.dv, 
                               treat = IV, 
                               mediator = MED,
                               covariates = coVar,
                               boot = boot)
    
    #### Extract results ####
    # WITHOUT random effects
    if(is.null(randVar)){
      result <- tidy(fit.totEffect) %>% 
        bind_rows(tidy(fit.mediator)) %>% 
        bind_rows(tidy(fit.dv)) %>% 
        mutate(model = c(rep("DV~IV", nrow(tidy(fit.totEffect))),
                         rep("MED~IV", nrow(tidy(fit.mediator))),
                         rep("DV~IV+MED", nrow(tidy(fit.dv))))) %>% 
        dplyr::select(model, everything())
    } else{
      # WITH random effects
      m.result <- tidy(fit.totEffect) %>% 
        bind_rows(tidy(fit.mediator)) %>% 
        bind_rows(tidy(fit.dv)) %>% 
        mutate(model = c(rep("DV~IV", nrow(tidy(fit.totEffect))),
                         rep("MED~IV", nrow(tidy(fit.mediator))),
                         rep("DV~IV+MED", nrow(tidy(fit.dv)))))
      #Estimate p-values
      p.result <- tidy(Anova(fit.totEffect)) %>% 
        bind_rows(tidy(Anova(fit.mediator))) %>% 
        bind_rows(tidy(Anova(fit.dv))) %>% 
        mutate(model = c(rep("DV~IV", nrow(tidy(Anova(fit.totEffect)))),
                         rep("MED~IV", nrow(tidy(Anova(fit.mediator)))),
                         rep("DV~IV+MED", nrow(tidy(Anova(fit.dv)))))) %>% 
        mutate(effect="fixed") %>% 
        rename(statistic.Anova = statistic)
      
      result <- full_join(m.result, p.result) %>% 
        dplyr::select(model, everything())
    }
    
    #Mediation analysis
    mediation.summ <- summary(mediation)
    
    result.mediation <- data.frame(
      #ACME = average causal mediation effects
      #ADE = average direct effects
      model="mediation",
      term = c("ACME","ADE","Total Effect","Prop. Mediated"),
      estimate = c(mediation.summ$d0, mediation.summ$z0,
                   mediation.summ$tau.coef, mediation.summ$n0),
      
      CI95_lower =c(mediation.summ$d0.ci[1], mediation.summ$z0.ci[1],
                    mediation.summ$tau.ci[1], mediation.summ$n0.ci[1]),
      
      CI95_upper = c(mediation.summ$d0.ci[2], mediation.summ$z0.ci[2],
                     mediation.summ$tau.ci[2], mediation.summ$n0.ci[2]),
      
      p.value = c(mediation.summ$d0.p, mediation.summ$z0.p,
                  mediation.summ$tau.p, mediation.summ$n0.p))
    
    #### Save results to list ####
    name <- paste(DV,IV,MED, sep=".")
    model.result[[name]] <- result
    mediation.result[[name]] <- result.mediation
    
    #### Plot ####
    if(plot){
      basic.plot(result=result, DV=DV, IV=IV, MED=MED,
                 result.mediation=result.mediation,
                 outdir=outdir, prefix=prefix, name=name)
    }
    
    }, error=function(e){
      print(paste("Could not complete mediation for ", 
                  DV, " ~ ", IV, " + ", MED))
    })
        }
  }
  
  #### Save tables to disk ####
  df.file1 <- paste(outdir, prefix, ".mediation.models.csv", sep="")
  df.file2 <- paste(outdir, prefix, ".mediation.summary.csv", sep="")
  
  plyr::ldply(model.result) %>% 
    dplyr::rename(DV.IV.MED=`.id`) %>% 
    write_csv(path=df.file1)
  
  plyr::ldply(mediation.result) %>% 
    dplyr::rename(DV.IV.MED=`.id`) %>% 
    write_csv(path=df.file2)
  
}

#### Basic triangle plot ####
basic.plot <- function(result=result, DV=DV, IV=IV, MED=MED,
                       result.mediation=result.mediation,
                       outdir=outdir, prefix=prefix, name=name, ...){
  require(DiagrammeR, quietly = TRUE)
  require(DiagrammeRsvg, quietly = TRUE)
  require(rsvg, quietly = TRUE)
  
  ##### Edge labels #####
  #create edge labels with model effect size and p-value symbol
  iv.med.lab <- result %>% 
    filter(model == "MED~IV" & term == IV) %>% 
    distinct(estimate, p.value) %>% 
    unlist(use.names=FALSE)
  iv.med.lab <- paste("IV effect = ", round(iv.med.lab[1], digits=3), 
                      " P = ", round(iv.med.lab[2], digits=3), sep="")
  
  med.dv.lab <- result %>% 
    filter(model == "DV~IV+MED" & term == MED) %>% 
    distinct(estimate, p.value) %>% 
    unlist(use.names=FALSE)
  med.dv.lab <- paste("MED effect controlled for IV = ", round(med.dv.lab[1], digits=3), 
                      " P = ", round(med.dv.lab[2], digits=3), sep="")
  
  iv.dv.lab1 <- result %>% 
    filter(model == "DV~IV" & term == IV) %>% 
    distinct(estimate, p.value) %>% 
    unlist(use.names=FALSE)
  iv.dv.lab1 <- paste("IV total effect = ", round(iv.dv.lab1[1], digits=3), 
                      " P = ", round(iv.dv.lab1[2], digits=3), sep="")
  
  iv.dv.lab2 <- result %>% 
    filter(model == "DV~IV+MED" & term == IV) %>% 
    distinct(estimate, p.value) %>% 
    unlist(use.names=FALSE)
  iv.dv.lab2 <- paste("IV effect controlled for MED = ", 
                      round(iv.dv.lab2[1], digits=3), 
                      " P = ", round(iv.dv.lab2[2], digits=3), sep="")
  
  iv.dv.lab <- paste(iv.dv.lab1, iv.dv.lab2, med.dv.lab, sep="\n")
  
  ##### Format df for diagram #####
  nodes <- create_node_df(n=3,
                          nodes = c("A","B","C"),
                          type="lower",
                          label = c(IV,MED,DV),
                          shape = "rectangle",
                          color = "black",
                          style = "filled", fillcolor="white",
                          fixedsize = FALSE,
                          fontcolor = "black")
  edges <- create_edge_df(from = c(1,1,2),
                          to = c(3,2,3),
                          color = "black",
                          label = c(iv.dv.lab, iv.med.lab, ""))
  
  graph <- create_graph(nodes_df = nodes,
                        edges_df = edges,
                        attr_theme = "lr")
  #### Save plot ####
  title.plot <- result.mediation %>% 
    distinct(term, estimate, p.value) %>% 
    mutate(p.value = ifelse(p.value < 0.001, "P<0.001", paste("P=",p.value, sep="")))
  title.plot <- paste(title.plot$term, round(title.plot$estimate, digits=3), 
                      title.plot$p.value, collapse = ", ")
  title.plot <- gsub(" NA$", "", title.plot)
  
  plot.dir <- paste(outdir, prefix, ".figs/", sep="")
  dir.create(plot.dir, showWarnings=FALSE, recursive=TRUE)
  plot.file <- paste(plot.dir, "/", name, ".png", sep="")
  
  export_graph(graph=graph, file_name=plot.file, title=title.plot,
               height=200, width=800)
}
  