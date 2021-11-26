phsc.plot.transmission.network<- function(df, di, point.size=10, point.size.couple=point.size*1.4, edge.gap=0.04, edge.size=0.4, curvature= -0.2, arrow=arrow(length=unit(0.02, "npc"), type="open"), curv.shift=0.08, label.size=3, node.label='ID', node.shape=NA_character_, node.fill=NA_character_, node.shape.values=c('M' = 15, 'F' = 17), node.fill.values=c('F'='hotpink2', 'M'='steelblue2') , threshold.linked=NA_real_)
{	
	#point.size=10; point.size.couple=14; edge.gap=0.04; edge.size=0.4; curvature= -0.2; arrow=arrow(length=unit(0.04, "npc"), type="open"); curv.shift=0.08; label.size=3
	#node.label='ID'; threshold.linked=0.6; node.shape=NA_character_; node.fill=NA_character_; node.shape.values=NA_integer_; node.fill.values=NA_character_
	#node.shape='IN_COUPLE'; node.fill='SEX'
	#node.shape.values=c('not in long-term\nrelationship'=18,'in long-term\nrelationship'=16)

  # point.size=10; point.size.couple=point.size*1.4; edge.gap=0.04; edge.size=0.2; curvature= -0.2;
  # arrow=arrow(length=unit(0.04, "npc"), type="open")
  # curv.shift=0.08; label.size=3; node.label='ID'
  # node.shape=NA_character_; node.fill=NA_character__; threshold.linked= NA_real_

	if(is.na(node.label))
	{
		node.label<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.label, NA_character_)
	}
	if(is.na(node.shape))
	{
		node.shape<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.shape, 'NA')
	}
	if(is.na(node.fill))
	{
		node.fill<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.fill, 'NA')
	}
	if(any(is.na(node.fill.values)))
	{
		z						<- unique(di[[node.fill]])
		node.fill.values		<- heat.colors(length(z))
		names(node.fill.values)	<- z
	}
	if(any(is.na(node.shape.values)))
	{
		z						<- unique(di[[node.shape]])
		node.shape.values		<- seq_along(z)
		names(node.shape.values)<- z
	}
	setnames(di, c(node.label, node.shape, node.fill), c('NODE_LABEL','NODE_SHAPE','NODE_FILL'))
	tmp	<- c('NODE_LABEL','NODE_SHAPE','NODE_FILL')[which(c(node.label, node.shape, node.fill)=='ID')]
	if(length(tmp))
		set(di, NULL, 'ID', di[[tmp]])
	di	<- subset(di, select=c(ID, NODE_LABEL, NODE_SHAPE, NODE_FILL))
	
	layout	<- as.data.table(ggnet2(network(unique(subset(df, select=c(ID1,ID2))), directed=FALSE, matrix.type="edgelist"))$data[,c("label", "x", "y")])
	setnames(layout, c('label','x','y'), c('ID1','ID1_X','ID1_Y'))
	df		<- merge(df, layout, by='ID1')
	setnames(layout, c('ID1','ID1_X','ID1_Y'), c('ID2','ID2_X','ID2_Y'))
	df		<- merge(df, layout, by='ID2')
	setnames(layout, c('ID2','ID2_X','ID2_Y'),  c('ID','X','Y'))	
	layout	<- merge(layout,di, by='ID')	

	df[, EDGETEXT_X:= (ID1_X+ID2_X)/2]
	df[, EDGETEXT_Y:= (ID1_Y+ID2_Y)/2]
	#
	#	calculate score for linked
	if(is.na(threshold.linked))
	{
		df	<- merge(df,df[, 	{
							z<- rep('edge_col_1', length(TYPE))
							z[which.max(POSTERIOR_SCORE)]	<- 'edge_col_2'
							list(EDGE_COL=z, TYPE=TYPE)	
						}, by=c('ID1','ID2')], by=c('ID1','ID2','TYPE'))		
	}
	if(!is.na(threshold.linked))
	{
		tmp	<- subset(df, TYPE!='not close/disconnected')[, list( EDGE_COL=as.character(factor(sum(POSTERIOR_SCORE)>=threshold.linked, levels=c(TRUE, FALSE), labels=c('edge_col_2','edge_col_1'))) ), by=c('ID1','ID2')]
		df	<- merge(df, tmp, by=c('ID1','ID2'))		
	}	
	#	for edges, move the start and end points on the line between X and Y
	#	define unit gradient
	df[, MX:= (ID2_X - ID1_X)]	
	df[, MY:= (ID2_Y - ID1_Y)]
	tmp		<- df[, sqrt(MX*MX+MY*MY)]
	set(df, NULL, 'MX', df[, MX/tmp])
	set(df, NULL, 'MY', df[, MY/tmp])	
	set(df, NULL, 'ID1_X', df[, ID1_X + MX*edge.gap])
	set(df, NULL, 'ID1_Y', df[, ID1_Y + MY*edge.gap])
	set(df, NULL, 'ID2_X', df[, ID2_X - MX*edge.gap])
	set(df, NULL, 'ID2_Y', df[, ID2_Y - MY*edge.gap])	
	#	label could just be move on the tangent vector to the line
	#	define unit tangent
	df[, TX:= -MY]
	df[, TY:= MX]
	tmp		<- df[, which(TYPE=='12')]
	set(df, tmp, 'EDGETEXT_X', df[tmp, EDGETEXT_X + TX*curv.shift])
	set(df, tmp, 'EDGETEXT_Y', df[tmp, EDGETEXT_Y + TY*curv.shift])
	tmp		<- df[, which(TYPE=='21')]
	set(df, tmp, 'EDGETEXT_X', df[tmp, EDGETEXT_X - TX*curv.shift])
	set(df, tmp, 'EDGETEXT_Y', df[tmp, EDGETEXT_Y - TY*curv.shift])
	#	

	p		<- ggplot() +			
			geom_point(data=layout, aes(x=X, y=Y, colour=NODE_FILL, pch=NODE_SHAPE), size=point.size) +
			geom_segment(data=subset(df, EDGE_COL=='edge_col_1' & TYPE=='ambiguous' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), lineend="butt") +
			geom_curve(data=subset(df, EDGE_COL=='edge_col_1' & TYPE=='12' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +
			geom_curve(data=subset(df, EDGE_COL=='edge_col_1' & TYPE=='21' & KEFF>0), aes(x=ID2_X, xend=ID1_X, y=ID2_Y, yend=ID1_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +
			geom_segment(data=subset(df, EDGE_COL=='edge_col_2' & TYPE=='ambiguous' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), lineend="butt") +
			geom_curve(data=subset(df, EDGE_COL=='edge_col_2' & TYPE=='12' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +
			geom_curve(data=subset(df, EDGE_COL=='edge_col_2' & TYPE=='21' & KEFF>0), aes(x=ID2_X, xend=ID1_X, y=ID2_Y, yend=ID1_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +									
			scale_colour_manual(values=c(node.fill.values, 'edge_col_1'='red', 'edge_col_2'='blue','NA'='green')) +
			scale_shape_manual(values=c(node.shape.values, 'NA'=21)) +
			scale_fill_manual(values=c(node.fill.values, 'NA'='grey50')) +
			scale_size_identity() +
			geom_text(data=subset(df, TYPE!='not close/disconnected' & KEFF>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=paste0(round(100*POSTERIOR_SCORE,d=1),'%')), size=label.size) +
			geom_text(data=layout, aes(x=X, y=Y, label=NODE_LABEL)) +
			theme_void() +
			guides(colour='none', fill='none',size='none', pch='none') 
	layout		<- subset(layout, select=c(ID,X,Y))
	setnames(layout, c('ID','X','Y'), c('label','x','y'))	
	# p$layout	<- layout
	return(p)

  # ggsave(paste0(gsub('-', '', Sys.Date()), '_network.png'), p)
}

plot_age_source_recipient <- function(data, title, lab, outdir){
  
  data <- data[!is.na(age_infection.SOURCE) & !is.na(age_infection.RECIPIENT)]
  data[, `Cohort round recipient` := cohort_round.RECIPIENT]
  data[, `Cohort round source` := cohort_round.SOURCE]
  data[, `Community recipient` := comm.RECIPIENT]
  data[, `Community source` := comm.SOURCE]
  data[, `Date infection recipient` := (date_first_positive.RECIPIENT - 1) <= as.Date('2017-01-01')]
  data[, `Date infection source` := (date_first_positive.SOURCE - 1) <= as.Date('2017-01-01')]
  
  # all pairs
  p <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point() + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) 
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_AllPairs_', lab, '.png')), w = 4, h = 4)
  
  # by cohort round
  p1 <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point(aes(col = `Cohort round source`)) + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')

  p2 <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point(aes(col = `Cohort round recipient`)) + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')

  p <- ggarrange(p1, p2, ncol = 2)
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_CohortRound_', lab, '.png')), w = 9, h = 7)
  
  # by age infection round
  p <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point(aes(col = `Date infection recipient`)) + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient',
         col = 'Date infection recipient before 2017') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_DateInfectionRecipient_', lab, '.png')), w = 5, h = 5)
  
  
  # by community
  p1 <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point(aes(col = `Community source`)) + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')

  p2 <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point(aes(col =`Community recipient`)) + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')
  
  p <- ggarrange(p1, p2, ncol = 2)
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_CommunityRecipient_', lab, '.png')), w = 9, h = 7)
  
}

plot_hist_age_infection <- function(pairs, outdir){
  
  pairs[, Sex := sex.SOURCE]
  p1 <- ggplot(pairs, aes(x = age_infection.SOURCE)) + 
    geom_histogram(bins = 30) + 
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection source') +
    scale_x_continuous(limits = range(c(pairs$age_infection.SOURCE, pairs$age_infection.RECIPIENT)))
  
  pairs[, Sex := sex.RECIPIENT]
  p2 <- ggplot(pairs, aes(x = age_infection.RECIPIENT)) + 
    geom_histogram(bins = 30) +     
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection recipient')  +
    scale_x_continuous(limits = range(c(pairs$age_infection.SOURCE, pairs$age_infection.RECIPIENT)))
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  
  file = file.path(outdir, paste0('hist_age_infection_', lab, '.png'))
  cat('saving', file)
  ggsave(p, file = file, w = 6, h = 6)
  
  return(p)
}

plot_time_infection <- function(pairs, outdir){
  
  pairs[, `Round source` := cohort_round.SOURCE]
  p1 <- ggplot(pairs, aes(x = date_first_positive.SOURCE)) + 
    geom_histogram(bins = 30) +
    facet_wrap(~`Round source`, nrow = 2) +
    theme_bw() + 
    labs(x = 'Date of infection source') + 
    geom_vline(xintercept = as.Date('2017-01-01'), linetype = 'dashed')
  
  pairs[, `Round recipient` := cohort_round.RECIPIENT]
  p2 <- ggplot(pairs, aes(x = date_first_positive.RECIPIENT)) + 
    geom_histogram(bins = 30) +
    facet_wrap(~`Round recipient`, nrow = 2) +
    theme_bw() + 
    labs(x = 'Date of infection recipient') + 
    geom_vline(xintercept = as.Date('2017-01-01'), linetype = 'dashed')

  p <- ggarrange(p1, p2, ncol = 2)
  
  file = file.path(outdir, paste0('hist_date_infection_', lab, '.png'))
  cat('saving', file)
  ggsave(p, file = file, w = 6, h = 6)
}

plot_hist_age_infection_diff_threshold <- function(pairs, outdir){
  
  chain <- keep.likely.transmission.pairs(as.data.table(dchain), 0.5) # Think I need to add chains.env prior to dchain
  pairs <- pairs.get.meta.data(chain, meta_data)
  pairs$threshold = '0.5'
  
  chain <- keep.likely.transmission.pairs(as.data.table(dchain), 0.6)
  pairs2 <- pairs.get.meta.data(chain, meta_data)
  pairs2$threshold = '0.6'
  
  pairs = rbind(pairs, pairs2)
  
  if(!include.mrc){
    cat('Keep only pairs in RCCS')
    pairs <- pairs[cohort.RECIPIENT == 'RCCS' & cohort.SOURCE == 'RCCS']
  }
  if(include.only.heterosexual.pairs){
    cat('Keep only heterosexual pairs')
    pairs <- pairs[(sex.RECIPIENT == 'M' & sex.SOURCE == 'F') | (sex.RECIPIENT == 'F' & sex.SOURCE == 'M')]
  }
  
  pairs[, Sex := sex.SOURCE]
  p1 <- ggplot(pairs, aes(x = age_infection.SOURCE)) + 
    geom_density(aes( group = threshold, fill = threshold), alpha = 0.5) + 
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection source') +
    scale_x_continuous(limits = range(c(pairs$age_infection.SOURCE, pairs$age_infection.RECIPIENT)))
  
  pairs[, Sex := sex.SOURCE]
  p2 <- ggplot(pairs, aes(x = age_infection.RECIPIENT)) + 
    geom_density(aes( group = threshold, fill = threshold), alpha = 0.5) + 
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection recipient')  +
    scale_x_continuous(limits = range(c(pairs$age_infection.SOURCE, pairs$age_infection.RECIPIENT)))
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  
  file = file.path(outdir, paste0('hist_age_infection_source_', gsub('(.+)_threshold.*', '\\1', lab), '.png'))
  cat('saving', file)
  ggsave(p, file = file, w = 6, h = 6)
  
  return(p)
}

plot_age_source <- function(pairs, outdir){
  
  tmp <- copy(pairs)
  tmp[, date_infection_before_2017.RECIPIENT := (date_first_positive.RECIPIENT - 1) < as.Date('2017-01-01')]
  tmp[, age_infection.SOURCE := floor(age_infection.SOURCE)]
  tmp[, age_infection.RECIPIENT := floor(age_infection.RECIPIENT)]
  tmp <- merge(tmp, df_age, by = c('age_infection.RECIPIENT', 'age_infection.SOURCE'))
  
  ps <- c(0.5, 0.2, 0.8)
  p_labs <- c('M','CL','CU')
  
  quantile(subset(tmp, sex.SOURCE == 'F' & sex.RECIPIENT == 'M' & age_infection_reduced.RECIPIENT == 46)$age_infection.SOURCE, prob=ps)
  
  tmp = tmp[, list(q= quantile(age_infection.SOURCE, prob=ps, na.rm = T), q_label=p_labs), 
            by=c('sex.SOURCE', 'sex.RECIPIENT', 'date_infection_before_2017.RECIPIENT', 'age_infection_reduced.RECIPIENT')]	
  tmp = dcast(tmp, sex.SOURCE + sex.RECIPIENT + date_infection_before_2017.RECIPIENT + age_infection_reduced.RECIPIENT ~ q_label, value.var = "q")
  
  subset(tmp, sex.SOURCE == 'F' & sex.RECIPIENT == 'M' & age_infection_reduced.RECIPIENT == 46)
  
  # FM
  tmp1 <- subset(tmp, sex.SOURCE == 'F' & sex.RECIPIENT == 'M') 
  p1 <- ggplot(tmp1, aes(x = age_infection_reduced.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_2017.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_2017.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Age at infection male recipient', y = 'Age at infection female source', col = 'Date infection of recipient before 2017') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 

  # MF
  tmp1 <- subset(tmp, sex.SOURCE == 'M' & sex.RECIPIENT == 'F') 
  p2 <- ggplot(tmp1, aes(x = age_infection_reduced.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_2017.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_2017.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Age at infection female recipient', y = 'Age at infection male source', col = 'Date infection of recipient before 2017') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 

  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  ggsave(p, filename = file.path(outdir, paste0('AgeInfectionSource_', lab, '.png')), w = 8, h = 5)
  
}

