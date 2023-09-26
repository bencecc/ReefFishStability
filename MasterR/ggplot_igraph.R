# Utility function to convert a network to a data frame for plotting with ggplot
# or to directly plot the network. The first option is appropriate to combine
# data for different metacommunity networks and to plot the with faceting. 
# The funciton also return network statistics (degree and cloeseness centrality
# that can be used to fill network nodes or set their size.
###############################################################################

ggplot_igraph <- function(plot.graph, MPA, net.name, mpa.extent, net.extent, conn.dec=3, plot=F, ...) {
	
	# compute centrality metrics (i.e. degree and closeness centrality)
	cc.tmp <- igraph::closeness(plot.graph, normalize=T)
	cc <- cc.tmp/sum(cc.tmp) # rescale
	dc.tmp <- igraph::degree(plot.graph, normalize=F, mode="all")
	dc <- dc.tmp/sum(dc.tmp) # rescale
	# minumum spanning treee
	# inet <- mst(inet.tmp)
	
	# set graphical parameters
	pal <- fish(2, option = "Cirrhilabrus_solorensis", begin=0.05, end=0.95, direction=1)
	V(plot.graph)$MPA <- MPA
	V(plot.graph)$color <- ifelse(V(plot.graph)$MPA=="No", pal[1], pal[2])
	V(plot.graph)$node.id <- as.numeric(vertex_attr(plot.graph)$name)
	V(plot.graph)$mpa.extent <- mpa.extent 
	V(plot.graph)$net.extent <- net.extent 
	V(plot.graph)$cc <- cc
	V(plot.graph)$dc <- dc
	
	my.layout <- as.data.frame(layout_with_kk(plot.graph, dim = 2))
	my.layout$color <- V(plot.graph)$color
	my.layout$node.id <- V(plot.graph)$node.id
	my.layout$mpa.extent <- V(plot.graph)$mpa.extent
	my.layout$net.extent <- V(plot.graph)$net.extent
	my.layout$cc <- V(plot.graph)$cc
	my.layout$dc <- V(plot.graph)$dc
	edge.info <- get.data.frame(plot.graph) 
	edge.info$from.x <- my.layout$V1[match(edge.info$from, my.layout$node.id)]  
	edge.info$from.y <- my.layout$V2[match(edge.info$from, my.layout$node.id)]
	edge.info$to.x <- my.layout$V1[match(edge.info$to, my.layout$node.id)]  
	edge.info$to.y <- my.layout$V2[match(edge.info$to, my.layout$node.id)]
	
	if(plot) {
		out <- ggplot(data=g3, aes(x=V1, y=V2))+
			geom_segment(data=edge.info, aes(x=from.x,xend = to.x, y=from.y,yend = to.y),
					linewidth=0.8, colour="darkgrey")+
			geom_point(aes(color=color), alpha=1, size=1.5)+
			labs(x="", y="",
					subtitle=paste(net.name, " (GC = ", round(mean(cc), conn.dec), ")", sep="")) +
			theme_bw() +
			facet_wrap(.~ID) +
			theme(
					legend.position="none",
					strip.background=element_rect(colour="grey90",
							fill="grey90"),
					plot.subtitle = element_text(size = 10),					
					panel.spacing.x=unit(0.1,"line"),
					panel.border=element_rect(colour="grey60", linewidth=0.2),
					panel.background = element_blank(),
					panel.grid.major = element_blank(), 
					panel.grid.minor = element_blank(),
					axis.text.x = element_blank(),
					axis.text.y = element_blank(),
					axis.ticks.x = element_blank(),
					axis.ticks.y = element_blank(),
					plot.margin = margin(0, 0.1, 0.1, 0.1, "cm")
					)
			
		} else {
			
			out <- cbind(ID=net.name, MPA=MPA, my.layout, rbind(edge.info, NA))
			
		}
#	clos.cent <- data.frame(ID=net.name, METRIC_TYPE="Closeness", centrality=cc)
#	deg.cent <- data.frame(ID=net.name, METRIC_TYPE="Degree", centrality=dc)
#			
#	out <- list(plot.net, clos.cent, deg.cent)
#	names(out) <- rep(c(paste(net.name, "Plot", sep="_")), 3)
	return(out)
}
