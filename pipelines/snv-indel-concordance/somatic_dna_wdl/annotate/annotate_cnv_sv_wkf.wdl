version 1.0

import "annotate.wdl" as annotate
import "../wdl_structs.wdl"

workflow AnnotateCnvSv {
    input {
        String tumor
        String normal
        String pairName
        Array[String] listOfChroms

        # cnv
        File bicseq2
        File cytoBand
        File dgv
        File thousandG
        File cosmicUniqueBed
        File cancerCensusBed
        File ensemblUniqueBed

        # sv
        File filteredMantaSV
        IndexedVcf gridssVcf
        String vepGenomeBuild

        # gap,DGV,1000G,PON,COSMIC
        File gap
        File dgvBedpe
        File thousandGVcf
        File svPon
        File cosmicBedPe

        # output
        File svBedPePath = "~{pairName}.sv.annotated.v7.somatic.final.bedpe"
        File svHighConfidenceBedPePath = "~{pairName}.sv.annotated.v7.somatic.high_confidence.final.bedpe"

        File svSupplementalBedPePath = "~{pairName}.sv.annotated.v7.somatic.supplemental.bedpe"
        File svHighConfidenceSupplementalBedPePath = "~{pairName}.sv.annotated.v7.somatic.high_confidence.supplemental.bedpe"

        Int annotateBicSeq2CnvMem = 36
    }

    call annotate.annotateBicSeq2Cnv {
        input:
            pairName=pairName,
            listOfChroms=listOfChroms,
            tumor=tumor,
            normal=normal,
            bicseq2=bicseq2,
            cytoBand=cytoBand,
            dgv=dgv,
            thousandG=thousandG,
            cosmicUniqueBed=cosmicUniqueBed,
            cancerCensusBed=cancerCensusBed,
            ensemblUniqueBed=ensemblUniqueBed,
            memoryGb = annotateBicSeq2CnvMem
    }

    call annotate.mergeSv {
        input:
            pairName=pairName,
            listOfChroms=listOfChroms,
            vepGenomeBuild=vepGenomeBuild,
            tumor=tumor,
            normal=normal,
            filteredMantaSV=filteredMantaSV,
            gridssVcf=gridssVcf

    }

    call annotate.annotateSv as annotateSvFinal {
        input:
            pairName=pairName,
            tumor=tumor,
            normal=normal,
            gap=gap,
            dgvBedpe=dgvBedpe,
            thousandGVcf=thousandGVcf,
            svPon=svPon,
            cosmicBedPe=cosmicBedPe,
            svMergedBedPe=mergeSv.svMergedFinalBedPe

    }

    call annotate.annotateSv as annotateSvSupplemental {
        input:
            pairName=pairName,
            tumor=tumor,
            normal=normal,
            gap=gap,
            dgvBedpe=dgvBedpe,
            thousandGVcf=thousandGVcf,
            svPon=svPon,
            cosmicBedPe=cosmicBedPe,
            svMergedBedPe=mergeSv.svMergedSupplementalBedPe

    }

    call annotate.annotateGenesSv {
        input:
            pairName=pairName,
            tumor=tumor,
            normal=normal,
            ensemblUniqueBed=ensemblUniqueBed,
            cancerCensusBed=cancerCensusBed,
            svMergedAnnotatedFinalBedPe=annotateSvFinal.svMergedAnnotatedBedPe
     }

    call annotate.annotateGenesSvSupplemental {
        input:
            pairName=pairName,
            tumor=tumor,
            normal=normal,
            ensemblUniqueBed=ensemblUniqueBed,
            cancerCensusBed=cancerCensusBed,
            svMergedAnnotatedSupplementalBedPe=annotateSvSupplemental.svMergedAnnotatedBedPe

    }

    call annotate.annotateWithCnvSv as annotateWithCnvSvFinal {
        input:
            pairName=pairName,
            tumor=tumor,
            normal=normal,
            cnvAnnotatedFinalBed=annotateBicSeq2Cnv.cnvAnnotatedFinalBed,
            svGeneAnnotatedBedPe=annotateGenesSv.svGeneAnnotatedFinalBedPe

    }

    call annotate.annotateWithCnvSv as annotateWithCnvSvSupplemental {
        input:
            pairName=pairName,
            tumor=tumor,
            normal=normal,
            cnvAnnotatedFinalBed=annotateBicSeq2Cnv.cnvAnnotatedFinalBed,
            svGeneAnnotatedBedPe=annotateGenesSvSupplemental.svGeneAnnotatedSupplementalBedPe

    }

    call annotate.filterBedPe as filterBedPeFinal {
        input:
            pairName=pairName,
            svCnvAnnotatedBedPe=annotateWithCnvSvFinal.svCnvAnnotatedBedPe,
            svBedPePath=svBedPePath,
            svHighConfidenceBedPePath=svHighConfidenceBedPePath

    }

    call annotate.filterBedPe as filterBedPeSupplemental {
        input:
            pairName=pairName,
            svCnvAnnotatedBedPe=annotateWithCnvSvSupplemental.svCnvAnnotatedBedPe,
            svBedPePath=svSupplementalBedPePath,
            svHighConfidenceBedPePath=svHighConfidenceSupplementalBedPePath

    }

    output {
        File svFinalBedPe = "~{filterBedPeFinal.svBedPe}"
        File svHighConfidenceFinalBedPe = "~{filterBedPeFinal.svHighConfidenceBedPe}"

        File svSupplementalBedPe = "~{filterBedPeSupplemental.svBedPe}"
        File svHighConfidenceSupplementalBedPe = "~{filterBedPeSupplemental.svHighConfidenceBedPe}"

        File cnvAnnotatedFinalBed  = "~{annotateBicSeq2Cnv.cnvAnnotatedFinalBed}"
        File cnvAnnotatedSupplementalBed  = "~{annotateBicSeq2Cnv.cnvAnnotatedSupplementalBed}"
    }
}
