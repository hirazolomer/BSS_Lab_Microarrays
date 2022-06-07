
# sensor systems lab 7 DNA microarrays

# install Bioconductor package
if(!requireNamespace("BiocManager",quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("affy")

# load affy library
library("affy")

# load phenodata.txt
phenodata = read.AnnotatedDataFrame(file='C:\\Users\\Rives\\OneDrive\\StudiumRives\\MSC\\Semester2\\Biomedical Sensor Systems Lab\\NR7 DNA Biosensors & Microarrays\\R tutorial\\GSE41328\\phenodata.txt')

# verify the read function
class(phenodata)

# access information from file
varMetadata(phenodata)
pData(phenodata)

# review sample names
sampleNames(phenodata)

# import .CEL files
celfiles = ReadAffy(filenames=sampleNames(phenodata),
                    celfile.path='C:\\Users\\Rives\\OneDrive\\StudiumRives\\MSC\\Semester2\\Biomedical Sensor Systems Lab\\NR7 DNA Biosensors & Microarrays\\R tutorial\\GSE41328',
                    phenoData=phenodata)


phenoData(celfiles)
pData(celfiles)
varLabels(celfiles)


image(celfiles[,1])
image(celfiles[,2])
image(celfiles[,3])
image(celfiles[,4])
image(celfiles[,5])
image(celfiles[,6])
image(celfiles[,7])
image(celfiles[,8])

# pm() intensities of perfect match probes
# head() prints the top row of a matrix
head(pm(celfiles))
dim(pm(celfiles))

# produce a histogram
hist(pm(celfiles[,1]))
hist(log2(pm(celfiles[,1])))

# access probe set identifiers
gn = geneNames(celfiles)
head(gn)

# length of vector
length(gn)

# extract number of corresponding identifiers, fir example "1552277_a_at"
p = probeset(celfiles, "1552277_a_at")
p
p[["1552277_a_at"]]
p[[1]]

# returns slot names
slotNames(p[["1552277_a_at"]])

# access PM intensities in p
p[["1552277_a_at"]]@pm

# PM and MM for "1007_s_at"
q = probeset(celfiles, "1007_s_at")
q[["1007_s_at"]]@pm

# visualize intensities for "1552277_a_at"
barplot.ProbeSet( p[["1552277_a_at"]], col.pm="yellow", col.mm="blue")

# visualize intensities for "1007_s_at"
barplot.ProbeSet( q[["1007_s_at"]], col.pm="yellow", col.mm="blue")



# examine a few more probes


# visualize the distribution of the data
hist(celfiles, which="both", col=1:8, lty=1, lwd=1)
# add a legend
legend("topright", rownames(pData(celfiles)), col=1:8, lty=1, lwd=1)
# create a boxplot
boxplot(celfiles, which="both", col="red", las=2)
