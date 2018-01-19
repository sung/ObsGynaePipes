library(data.table)
library(ggpubr)
source("~/Pipelines/config/graphic.R")
time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p')
		########################
		## Mapping efficiency ##
        ## Low-PSG7 sample    ##
		########################
		#https://www.r-bloggers.com/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
		# Scatter plot colored by Library
		dt.map<-fread("~/results/RNA-Seq/Plasma.2017/Meta/PE75.PE150.read.cnt.txt",stringsAsFactors=TRUE); dt.map$Type<-factor(dt.map$Type,levels=rev(levels(dt.map$Type)))
		dt.foo<-dt.map[Library %in% c("SLX-9342","SLX-9345") & Type=="PE75",.(Library,Seq=Yield/10^6,Eff=Mapping_Efficiency)]
		p.scatter<- ggscatter(dt.foo, x="Seq",y="Eff",color = "Library", palette = "jco",size = 6, alpha = 0.5)

		# Box plot of the x variable
		xbp <- ggboxplot(dt.foo$Seq, width = 0.2, fill = "lightgray") + rotate() + theme_transparent()
		# Box plot of the y variable
		ybp <- ggboxplot(dt.foo$Eff, width = 0.2, fill = "lightgray") + theme_transparent()
		# Create the external graphical objects
		# called a "grop" in Grid terminology
		xbp_grob <- ggplotGrob(xbp)
		ybp_grob <- ggplotGrob(ybp)
		# Place box plots inside the scatter plot
		#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		xmin <- min(dt.foo$Seq); xmax <- max(dt.foo$Seq)
		ymin <- min(dt.foo$Eff); ymax <- max(dt.foo$Eff)
		yoffset <- (1/15)*ymax; xoffset <- (1/15)*xmax
		# Insert xbp_grob inside the scatter plot
		p1 <- p.scatter + 
			annotation_custom(grob = xbp_grob, xmin = xmin, xmax = xmax, ymin = ymin-yoffset, ymax = ymin+yoffset) +
			annotation_custom(grob = ybp_grob, xmin = xmin-xoffset, xmax = xmin+xoffset, ymin = ymin, ymax = ymax) # Insert ybp_grob inside the scatter plot


		# Marginal density plot of x (top panel) and y (right panel)
		xplot <- ggdensity(dt.foo, "Seq", fill = "Library", palette = "jco")
		yplot <- ggdensity(dt.foo, "Eff", fill = "Library", palette = "jco") + rotate()
		# Cleaning the plots
		yplot <- yplot + clean_theme() 
		xplot <- xplot + clean_theme()
		# Arranging the plot
		p2<- ggarrange(xplot, NULL, p.scatter, yplot, ncol = 2, nrow = 2,  align = "hv", widths = c(2, 1), heights = c(1, 2), common.legend = TRUE)
