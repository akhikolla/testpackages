\name{view1}
\alias{view1}
\title{View1 data}
\usage{view1}
\description{
	 One of the two views simulated to mimic datasets from a real study, in which genes are characterized with expression patterns. 
}
\format{
	999x12 matrix with rows representing genes and columns representing various sources, such as different tissues or time points.
}
\details{
	We mimicked datasets from a real study in which genes are characterized with expression patterns and simulated datasets with implanted block structures that give both clusters of subjects and variables. Two views of data for 1000 subjects were created. There were 12 variables in view 1, and 15 variables in view 2. The data matrix of each view is created by randomly setting 0 or 1 to each entry with varying probability that is determined according to prefixed block structures projected in the data. More specifically, we start from a data matrix filled with all 0. Then we reset data entries inside and outside the blocks to 1 with probability 0.9 and 0.1, respectively. For simplifying the process and easy presentation, we had subjects in the two datasets well aligned and indexed from 1 to 1000; and variables were also indexed using consecutive number starting from 1. View 1 was designed to have two blocks. The first block consists of subjects from 1 to 400 and variables from 1 to 3. The second includes the 200 subjects indexed from 481 to 680 and variables from 4 to 6. Three blocks were included in view 2. The first block contains subjects from 1 to 240 and the first three variables. The second block consists of subjects from 241 to 480 and variable 4, 5 and 6. The last block includes 320 subjects indexed from 481 to 800 and variables from 7 to 9. By comparing blocks of the two views, it is obvious that there are three consistent blocks (i.e., containing same subjects) among the two views. Block 1 consists of 240 subjects and associates with variables from 1 to 3 in both view. There are 200 subjects in block 2. The associated variables are 4, 5 and 6 in view 1 and 7, 8 and 9 in view 2. Block 3 consists of 160 subjects and associates with variable from 1 to 3 in view 1 and from 4 to 6 in view 2.   
}