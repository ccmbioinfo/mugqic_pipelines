
group <- c()
samps <- c()
count <- 0
for (num in contrast) {
        count <- count + 1
	if (num == 1) {
                group <- c(group, "CTRL")
		samps <- c(samps, decoder$Sample[count])
                }
        if (num == 2) {
                group <- c(group, "CASE")
		samps <- c(samps, decoder$Sample[count])
                }
}

countFiles <- paste0("jctseq/rawCts/", samps, "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz");

jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
     sample.names = samps,
     condition=factor(group),
     flat.gff.file = gff.file,
     nCores = 1,
     analysis.type = "junctionsAndExons"
     );

writeCompleteResults(jscs=jscs,
        outfile.prefix = paste0("./", folder, "/results/"),
        save.jscs = TRUE
        );

buildAllPlots(jscs=jscs,
        outfile.prefix = paste0("./", folder, "/plots/"),
        use.plotting.device = "png",
        FDR.threshold = 0.01
        );
