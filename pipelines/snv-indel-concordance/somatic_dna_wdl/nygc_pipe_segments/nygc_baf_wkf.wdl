version 1.0
import "../wdl_structs.wdl"
import "../tasks/utils.wdl" as utils
import "../pre_process/qc.wdl" as qc
import "../pre_process/conpair_wkf.wdl" as conpair
import "../baf/baf_wkf.wdl" as baf
import "../alignment_analysis/msi_wkf.wdl" as msi
import "../alignment_analysis/telomeasures_wkf.wdl" as telomeasures

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2024) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    James Roche (jroche@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================


workflow NygcBaf {
    input {
        Boolean local = false
        IndexedReference referenceFa
        # patient-specific inputs
        PairInfo pairInfo
        # unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf
        File haplotypecallerAnnotatedVcf
        # need to create docker with gcloud
        File? serviceAccountKey
        String? gcpProject
        File? alleleCountsTxtPreexist
    }
    if (!defined(alleleCountsTxtPreexist)) {
        call baf.Baf {
            input:
                referenceFa = referenceFa,
                local = local,
                pairName = pairInfo.analysisPairId,
                sampleId = pairInfo.normalId,
                tumorFinalBam = pairInfo.tumorFinalBam,
                normalFinalBam = pairInfo.normalFinalBam,
                germlineVcf = haplotypecallerAnnotatedVcf,
                serviceAccountKey = serviceAccountKey,
                gcpProject = gcpProject
        }
    }
    output {
        File alleleCountsTxt = select_first([Baf.alleleCountsTxt, alleleCountsTxtPreexist])
    }
}
