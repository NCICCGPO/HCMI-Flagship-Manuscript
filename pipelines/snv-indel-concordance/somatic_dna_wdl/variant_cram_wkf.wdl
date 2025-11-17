version 1.0
import "wdl_structs.wdl"
import "tasks/utils.wdl" as utils


workflow CollectVariantSupport {
    input {
        IndexedReference referenceFa
        # -> Patient-specific inputs
        PairInfo pairInfo
        File mainVcf
        File cnvAnnotatedFinalBed
        File svFinalBedPe
        # needed for work on cloud without localizing
        File? serviceAccountKey
        String? gcpProject
    }
    # -> get insert size estimate if not provided
    Int tumorDiskSize = ceil(size(pairInfo.tumorFinalBam.bam, "GB")) + 30
    Int normalDiskSize = ceil(size(pairInfo.normalFinalBam.bam, "GB")) + 30
    # create resources to replace full BAMs for most review
    call utils.FineGrainedCov as tumorFineGrainedCov {
        input:
            finalBam = pairInfo.tumorFinalBam,
            sampleId = pairInfo.tumorId,
            diskSize = tumorDiskSize + 50
    }
    call utils.FineGrainedCov as normalFineGrainedCov {
        input:
            finalBam = pairInfo.normalFinalBam,
            sampleId = pairInfo.normalId,
            diskSize = normalDiskSize + 50
    }
    call utils.MakeVariantCram as tumorMakeVariantCram {
        input:
            finalBam = pairInfo.tumorFinalBam,
            sampleId = pairInfo.tumorId,
            mainVcf = mainVcf,
            cnvAnnotatedFinalBed = cnvAnnotatedFinalBed,
            svFinalBedPe = svFinalBedPe,
            referenceFa = referenceFa,
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject,
            diskSize = 50
    }
    call utils.MakeVariantCram as normalMakeVariantCram {
        input:
            finalBam = pairInfo.normalFinalBam,
            sampleId = pairInfo.normalId,
            mainVcf = mainVcf,
            cnvAnnotatedFinalBed = cnvAnnotatedFinalBed,
            svFinalBedPe = svFinalBedPe,
            referenceFa = referenceFa,
            serviceAccountKey = serviceAccountKey,
            gcpProject = gcpProject,
            diskSize = 50
    }
    output {
        File tumorDistMosdepth = tumorFineGrainedCov.distMosdepthRegion
        File tumorSummaryMosdepth = tumorFineGrainedCov.summaryMosdepth
        File tumorRegionsMosdepth = tumorFineGrainedCov.regionsMosdepth
        File tumorRegionsMosdepthIndex = tumorFineGrainedCov.regionsMosdepthIndex
        File normalDistMosdepth = normalFineGrainedCov.distMosdepthRegion
        File normalSummaryMosdepth = normalFineGrainedCov.summaryMosdepth
        File normalRegionsMosdepth = normalFineGrainedCov.regionsMosdepth
        File normalRegionsMosdepthIndex = normalFineGrainedCov.regionsMosdepthIndex
        Cram tumorVariantCram = tumorMakeVariantCram.variantCram
        Cram normalVariantCram = normalMakeVariantCram.variantCram
    }
}