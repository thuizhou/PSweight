#' Plot the distribution of propensity scores and balance statistics
#'
#' Summarize the SumStat x, generate histogram or density of estimated propensity scores and
#' plot the balance statistics under weighting versus no weighting.
#'
#' @param x a \code{SumStat} object obtained with \code{\link{SumStat}} function.
#' @param type a character indicating the type of plot to produce, including histogram of estimated propensity scores (\code{"hist"}), density of estimated propensity scores (\code{"density"}), and plot of balance statistics (\code{"balance"}).
#' @param weighted.var logical. Indicating whether weighted variance should be used in calculating the balance statistics. Default is \code{TRUE}.
#' @param threshold an optional numeric value indicating the balance threshold for the balance plot. Default is 0.1. Only valid when \code{type = "balance"}.
#' @param metric a character indicating the type of metric used in balance plot. Only \code{"ASD"} or \code{"PSD"} is allowed. If not specified, the default is \code{"ASD"}. See  \code{\link{summary.SumStat}} for additional details on balance metrics.
#' @param ... further arguments passed to or from other methods.
#'
#' @details For the balance plot, a vertical line at \code{threshold} is used to define balance on covariates.
#' The default value is \code{threshold = 0.1} following Austin and Stuart (2015). If more than 2 treatments
#' are considered, only density of the estimated generalized propensity scores will be produced, regardless of
#' whether \code{type = "density"} or \code{type = "hist"}.
#'
#' @return Plot of the indicated type.
#'
#' @references
#'
#' Austin, P.C. and Stuart, E.A. (2015). Moving towards best practice when using inverse probability of treatment weighting (IPTW) using the propensity score to estimate causal treatment effects in observational studies.
#' Statistics in Medicine, 34(28), 3661-3679.
#'
#' @export
#'
#' @examples
#' data("psdata")
#' ps.formula<-trt~cov1+cov2+cov3+cov4+cov5+cov6
#' msstat <- SumStat(ps.formula, trtgrp="2", data=psdata,
#'    weight=c("ATE","ATO","ATT"))
#'
#' plot(msstat, type="hist")
#' plot(msstat, type="balance", weighted.var=TRUE, threshold=0.1, metric="ASD")
#'
#' @import ggplot2
#' @importFrom  stats binomial coef cov formula glm lm model.matrix plogis poisson predict qnorm quantile sd
#' @importFrom  utils capture.output combn
#' @importFrom  graphics hist legend
#'
#'
plot.SumStat<-function(x, type="balance", weighted.var=TRUE, threshold=0.1, metric="ASD",...){
  #get object info
  m<-length(names(x$ps.weights)[-c(1,2)])+1
  zname<-names(x$ps.weights)[1]
  Z<-unlist(x$ps.weights[zname])
  ncate<-length(unique(Z))
  metric<-toupper(metric)

  #extract original treatment labels
  dic0<-rep(NA,ncate)
  for (i in 1:ncate){
    dic0[i]<-as.character(x$ps.weights[which(x$ps.weights$zindex==i)[1],1])
  }

  if (type=="balance"){
    if(!(metric %in% c("ASD","PSD"))){
      warning("Balance metric unrecognized, 'ASD' calculated instead.")
      metric="ASD"
    }
    wt_list<-c("unweighted",names(x$ps.weights)[-c(1,2)])
    p<-dim(x$unweighted.sumstat$ASD.unweighted.var)[1]
    wt.type<-rep(wt_list,each=p)
    var.name<-rep(rownames(x$unweighted.sumstat$ASD.unweighted.var),m)

    if(metric=="ASD"){
      ab.smd=c()
      for (wt in wt_list){
        wt_ex<-paste0(wt,".sumstat")
        if(weighted.var==TRUE){
          ab.smd<-c(ab.smd,apply(abs(x[[wt_ex]]$ASD.weighted.var),1,max))
        } else{
          ab.smd<-c(ab.smd,apply(abs(x[[wt_ex]]$ASD.unweighted.var),1,max))
        }
      }

      plot.df<-data.frame(wt.type,var.name,ab.smd)
      xl = min(plot.df$ab.smd, 0)
      xu = max(plot.df$ab.smd, 1)

      #SMD Plot
      pt<-ggplot(plot.df, aes(x=ab.smd, y=var.name))+theme_bw()+
        geom_point(aes(color=wt.type, shape=wt.type, stroke=1.5), size=4)+
        #scale_color_manual(values=rep("black", m))+
        scale_x_continuous(limits=c(xl,xu), breaks = c(0,threshold,1))+
        scale_shape_manual(values=1:m)+
        scale_y_discrete(limits = rev(levels(plot.df$var.name)))+
        xlab("Standardized Mean Differences") + ylab("")+
        geom_vline(xintercept=0, linetype="dashed")+
        theme(plot.title = element_text(hjust = 0.5, size=18),
              axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),
              axis.title.x=element_text(size=16), legend.text=element_text(size=16),
              legend.title=element_text(size=16))+
        geom_vline(xintercept = threshold, linetype="dotted",lwd=1.1)
      print(pt)
    }else if(metric=="PSD"){
      ab.smd=c()
      for (wt in wt_list){
        wt_ex<-paste0(wt,".sumstat")
        if(weighted.var==TRUE){
          ab.smd<-c(ab.smd,apply(abs(x[[wt_ex]]$PSD.weighted.var),1,max))
        }else{
          ab.smd<-c(ab.smd,apply(abs(x[[wt_ex]]$PSD.unweighted.var),1,max))
        }
      }
      #SMD forest plot
      plot.df<-data.frame(wt.type,var.name,ab.smd)
      xl = min(plot.df$ab.smd, 0)
      xu = max(plot.df$ab.smd, 1)

      #SMD Plot
      pt<-ggplot(plot.df, aes(x=ab.smd, y=var.name))+theme_bw(base_size=18)+
        geom_point(aes(color=wt.type, shape=wt.type, stroke=1.5), size=4)+
        #scale_color_manual(values=rep("black", m))+
        scale_x_continuous(limits=c(xl,xu), breaks = c(0,threshold,1))+
        scale_shape_manual(values=1:m)+
        scale_y_discrete(limits = rev(levels(plot.df$var.name)))+
        xlab("Standardized Mean Differences") + ylab("")+
        geom_vline(xintercept=0, linetype="dashed")+
        theme(plot.title=element_text(hjust=0.5, size=18),
              axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),
              axis.title.x=element_text(size=16), legend.text=element_text(size=16),
              legend.title=element_text(size=16))+
        geom_vline(xintercept = threshold, linetype="dotted",lwd=1.1)
      print(pt)
    }
  }

  if (type=="density"){

    #plot/check the density of prepensity score
    z<-x$ps.weights$zindex
    propensity<-x$propensity

    for(i in 1:ncate){
      proptmp<-propensity[,i]
      group<-as.character(Z)
      df<-data.frame(proptmp,group,zindex=z)
      pt<-ggplot(df,aes(x=proptmp,linetype=group, colour=group))+
        geom_density(size=1.2,show.legend=FALSE)+
        stat_density(aes(x=proptmp, linetype=group, colour=group),
                     geom="line",position="identity",size=0)+
        scale_color_brewer(palette="Set1")+
        # scale_linetype_manual(name=zname, values=c(1:ncate))+
        xlab(paste("Propensity score for group",dic0[i]))+
        theme_bw(base_size = 18)+
        theme(axis.line = element_line(),
              #panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank())+
        guides(colour = guide_legend(override.aes=list(size=1.0)))

      #pt
      print(pt)
      message(paste("Propensity score for group",dic0[i]))
      if (i<ncate) {
        invisible(readline(prompt="Press [enter] to continue"))
      }
    }
  }
  if (type=="hist"){
    if (ncate>2) {
      #plot/check the density of prepensity score
      z<-x$ps.weights$zindex
      propensity<-x$propensity

      warning("Histogram only available for binary treatment. Density plot provided instead.")

      for(i in 1:ncate){
        proptmp<-propensity[,i]
        group<-as.character(Z)
        df<-data.frame(proptmp,group,zindex=z)
        pt<-ggplot(df,aes(x=proptmp,linetype=group, colour=group))+
          geom_density(size=1.2,show.legend=FALSE)+
          stat_density(aes(x=proptmp, linetype=group, colour=group),
                       geom="line",position="identity",size=0)+
          scale_color_brewer(palette="Set1")+
          # scale_linetype_manual(name=zname, values=c(1:ncate))+
          xlab(paste("Propensity score for group",dic0[i]))+
          theme_bw(base_size = 18)+
          theme(axis.line = element_line(),
                #panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank())+
          guides(colour = guide_legend(override.aes=list(size=1.0)))

        #pt
        print(pt)
        message(paste("Propensity score for group",dic0[i]))
        if (i<ncate) {
          invisible(readline(prompt="Press [enter] to continue"))
        }
      }
    }
    else{
      e<-x$propensity[,2]
      z<-x$ps.weights$zindex
      hist(e[z==1],breaks=50,col="gray77",border="gray77",
           main=NULL,ylab=NULL,xlim=c(max(min(e)-0.2,0),min(max(e)+0.2,1)),xlab="Estimated propensity score",
           freq=F,cex.lab = 1.5, cex.axis = 1.5 ,cex.main = 2, cex = 2)
      hist(e[z==2],breaks=50,add=TRUE,freq=F)
      legend("topright",title="group",legend=dic0,col=c("gray77","black"),
             lty=1,lwd=1.5,bty='n',cex=1.5)
    }
  }
}




