# used to merge and plot the number of publications regarding machine learning 
# in different types of omics

# the raw csv files were downloaded from pubmed.ncbi.nlm.nih.gov as the date of August 7, 2021
# search url is https://pubmed.ncbi.nlm.nih.gov/?term=machine+learning%2C+[omics_type]&filter=years.2006-2021
# where the omics_type refer to the name of the omics which is given in each csc file
# Note: we restrict the data between years 2006-2021; but some omics studies on
# have publications in later than 2006, and only those years considered.

require("tidyverse")

# set working dir as source file
# for Rstudio IDE
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# for otherwise
# setwd(getSrcDirectory()[1])
omics_type = c("Genomics", "Transcriptomics", "Proteomics", "Metabolomics","Single_cell",
                "Imaging_omics", "Multiomics")
#loading data from csv files and merge them
allDf = NULL
for (omics in omics_type) {
  
  #omics = "Imaging_omics"
  inFile = paste0(omics,"_PubMed_Timeline_Results_by_Year.csv")
  inDf = read.csv(inFile, skip = 1) # skip the first line of the file
  colnames(inDf)[2] <- omics
  if(is.null(allDf)){
    allDf = inDf
  }else{
    allDf = allDf %>% full_join(inDf, by = "Year")
  }
}
rm(inDf) # delete the unused variables

# plot the Circular barplot

# convert the data from wide to long format
data <- allDf %>% pivot_longer(!Year, names_to = "group", values_to = "count")
data$value = log(data$count) # convert the count to log2, so that scale down the diversity of different omics


# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(match(group, omics_type)) # group the table by omics type in the same order in manuscript sections
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id)) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))%>% 
  arrange(match(group, omics_type))

base_angle <- label_data %>% group_by(group) %>% 
  summarize(angle = median(angle))
base_data$angle <- ifelse(base_angle$angle > 0, base_angle$angle-90, base_angle$angle-60)

angle = -360 * base_data$title/nrow(data)     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
base_data$angle <- ifelse(angle < -90 , angle+180, ifelse(angle < -180 , angle-180, ifelse(angle < -270 , angle+360, angle)))
base_data$angle <- ifelse(angle < -270 , angle+360, ifelse(angle < -180 , angle-180, ifelse(angle < -90 , angle+180, angle)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

scales = c(2, 4, 6, 8)
# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = scales[4], xend = start, yend = scales[4]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = scales[3], xend = start, yend = scales[3]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = scales[2], xend = start, yend = scales[2]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = scales[1], xend = start, yend = scales[1]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = scales, label = scales , color="grey", size=2.5 , angle=0, fontface="bold", hjust=1) +
  annotate("text", x = max(data$id)*0.985, y = (max(scales)+min(scales))/2, label = "log2(count)" , color="grey", size=2.5 , angle=95, fontface="bold", hjust=0.5) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-min(scales)*4,max(scales)+min(scales)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  )+
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+min(scales)/5, label=count, hjust=hjust), 
            color="black", alpha=0.6, size=2, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_text(data=label_data, aes(x=id, y=-min(scales)/2, label=Year, hjust=hjust), 
            color="black", alpha=0.6, size=1.5, angle= label_data$angle, inherit.aes = FALSE ) +
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -min(scales)*0.7, xend = end, yend = -min(scales)*0.7), colour = "black", alpha=1, size=0.6 , inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x = title, y = -min(scales)*1.1 , label=group), colour = "black", alpha=1, size=2.5, angle = base_data$angle, inherit.aes = FALSE)

pdf("pltos.pdf", 6,6)
p + annotate("text", x=0, y=-min(scales)*4, label= "Number of publications\n for Machine learning\n in omics", size = 3, fontface = "bold")
dev.off()

