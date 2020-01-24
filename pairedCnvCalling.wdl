version 1.0

# Copyright (c) 2020 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import "sample.wdl" as sampleWorkflow

workflow PairedCnvCalling {
    input {
        String caseSampleName
        File caseBam
        File caseBamIndex
        String controlSampleName
        File controlBam
        File controlBamIndex
        File? PON
        File? annotatedIntervals
        File preprocessedIntervals
        File commonVariantSites
        String outputDir = "."
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai

        Map[String, String] dockerImages = {
            "gatk": "quay.io/biocontainers/gatk4:4.1.2.0--1"
        }
    }

    call sampleWorkflow.Sample as controlSample {
        input:
            preprocessedIntervals = preprocessedIntervals,
            PON = PON,
            annotatedIntervals = annotatedIntervals,
            inputBam = controlBam,
            inputBamIndex = controlBamIndex,
            commonVariantSites = commonVariantSites,
            sampleName = controlSampleName,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            outputDir = outputDir + "/" + controlSampleName + "/",
            dockerImages = dockerImages
    }

    call sampleWorkflow.Sample as caseSample {
        input:
            preprocessedIntervals = preprocessedIntervals,
            PON = PON,
            annotatedIntervals = annotatedIntervals,
            matchedNormalAllelicCounts = controlSample.allelicCounts,
            inputBam = caseBam,
            inputBamIndex = caseBamIndex,
            commonVariantSites = commonVariantSites,
            sampleName = caseSampleName,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            outputDir = outputDir + "/" + caseSampleName + "/",
            dockerImages = dockerImages
    }

    output {
        File caseAllelicCounts = caseSample.allelicCounts
        File caseReadCounts = caseSample.readCounts
        File caseStandardizedCopyRatios = caseSample.standardizedCopyRatios
        File caseDenoisedCopyRatios = caseSample.denoisedCopyRatios
        File caseHetrozygousAllelicCounts = caseSample.hetrozygousAllelicCounts
        File caseNormalHetrozygousAllelicCounts = caseSample.normalHetrozygousAllelicCounts
        File caseCopyRatioSegments = caseSample.copyRatioSegments
        File caseCopyRatioCBS = caseSample.copyRatioCBS
        File caseAlleleFractionCBS = caseSample.alleleFractionCBS
        File caseUnsmoothedModeledSegments = caseSample.unsmoothedModeledSegments
        File caseUnsmoothedCopyRatioParameters = caseSample.unsmoothedCopyRatioParameters
        File caseUnsmoothedAlleleFractionParameters = caseSample.unsmoothedAlleleFractionParameters
        File caseModeledSegments = caseSample.modeledSegments
        File caseCopyRatioParameters = caseSample.copyRatioParameters
        File caseAlleleFractionParameters = caseSample.alleleFractionParameters
        File caseCalledSegments = caseSample.calledSegments
        File caseCalledSegmentsIgv = caseSample.calledSegmentsIgv
        File caseDenoisedCopyRatiosPlot = caseSample.denoisedCopyRatiosPlot
        File caseDenoisedCopyRatiosLimitedPlot = caseSample.denoisedCopyRatiosLimitedPlot
        File caseStandardizedMedianAbsoluteDeviation = caseSample.standardizedMedianAbsoluteDeviation
        File caseDenoisedMedianAbsoluteDeviation = caseSample.denoisedMedianAbsoluteDeviation
        File caseDeltaMedianAbsoluteDeviation = caseSample.deltaMedianAbsoluteDeviation
        File caseDeltaScaledMedianAbsoluteDeviation = caseSample.deltaScaledMedianAbsoluteDeviation
        File caseModeledSegmentsPlot = caseSample.modeledSegmentsPlot

        File controlAllelicCounts = controlSample.allelicCounts
        File controlReadCounts = controlSample.readCounts
        File controlStandardizedCopyRatios = controlSample.standardizedCopyRatios
        File controlDenoisedCopyRatios = controlSample.denoisedCopyRatios
        File controlHetrozygousAllelicCounts = controlSample.hetrozygousAllelicCounts
        File controlNormalHetrozygousAllelicCounts = controlSample.normalHetrozygousAllelicCounts
        File controlCopyRatioSegments = controlSample.copyRatioSegments
        File controlCopyRatioCBS = controlSample.copyRatioCBS
        File controlAlleleFractionCBS = controlSample.alleleFractionCBS
        File controlUnsmoothedModeledSegments = controlSample.unsmoothedModeledSegments
        File controlUnsmoothedCopyRatioParameters = controlSample.unsmoothedCopyRatioParameters
        File controlUnsmoothedAlleleFractionParameters = controlSample.unsmoothedAlleleFractionParameters
        File controlModeledSegments = controlSample.modeledSegments
        File controlCopyRatioParameters = controlSample.copyRatioParameters
        File controlAlleleFractionParameters = controlSample.alleleFractionParameters
        File controlCalledSegments = controlSample.calledSegments
        File controlCalledSegmentsIgv = controlSample.calledSegmentsIgv
        File controlDenoisedCopyRatiosPlot = controlSample.denoisedCopyRatiosPlot
        File controlDenoisedCopyRatiosLimitedPlot = controlSample.denoisedCopyRatiosLimitedPlot
        File controlStandardizedMedianAbsoluteDeviation = controlSample.standardizedMedianAbsoluteDeviation
        File controlDenoisedMedianAbsoluteDeviation = controlSample.denoisedMedianAbsoluteDeviation
        File controlDeltaMedianAbsoluteDeviation = controlSample.deltaMedianAbsoluteDeviation
        File controlDeltaScaledMedianAbsoluteDeviation = controlSample.deltaScaledMedianAbsoluteDeviation
        File controlModeledSegmentsPlot = controlSample.modeledSegmentsPlot
    }

    parameter_meta {
        caseSampleName: {description: "The name of the case sample.", category: "required"}
        caseBam: {description: "The BAM file for the case sample.", category: "required"}
        caseBamIndex: {description: "The index for the case sample's BAM file.", category: "required"}
        controlSampleName: {description: "The name of the control sample.", category: "required"}
        controlBam: {description: "The BAM file for the control sample.", category: "required"}
        controlBamIndex: {description: "The index for the control sample's BAM file.", category: "required"}
        preprocessedIntervals: {description: "Intervals to operate on. Should be produced by gatk PreprocessIntervals (eg. using CNV-PON.wdl).",
                                category: "required"}
        PON: {description: "A read counts panel of normals as generated by gatk CreateReadCountPanelOfNormals (eg. using CNV-PON.wdl).",
                            category: "common"}
        annotatedIntervals: {description: "An annotated set of intervals as generated by AnnotateIntervals (eg. using CNV-PON.wdl). Will be ignored if PON is provided.",
                             category: "common"}
        matchedNormalAllelicCounts: {description: "The allelicCounts as generate by CollectAllelicCounts for a matched normal.", category: "common"}
        commonVariantSites: {description: "Interval list of common variamt sites (to retrieve the allelic counts for).", category: "required"}
        sampleName: {description: "The name of the sample, used for file naming.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        outputDir: {description: "The directory the output should be written to.", category: "common"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.",
                       category: "advanced"}
    }
}