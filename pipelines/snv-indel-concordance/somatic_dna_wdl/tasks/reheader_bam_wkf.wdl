version 1.0

import "../wdl_structs.wdl"
import "../alignment_analysis/alignment_analysis.wdl" as alignmentAnalysis
import "../tasks/utils.wdl" as utils

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2022) by the New York
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

workflow Reheader {
    input {
        
        Bam finalBam
        String sampleId
        IndexedReference referenceFa
    
    }

    # using small disk size because the file is not localized (on servers that support this)
    Int basicDiskSize = 4
    
    call alignmentAnalysis.GetSampleNameFile {
        input:
            finalBam = finalBam.bam,
            finalBai = finalBam.bamIndex,
            fileId = sampleId,
            referenceFa = referenceFa,
            diskSize = basicDiskSize
    }
    
    call utils.CompareStrings {
        input:
            fileInputString = GetSampleNameFile.bamSampleId,
            inputString = sampleId,
            fileId = sampleId
    }
    
    if (size(CompareStrings.resultFile) == 0) {
        Int renameDiskSize = (ceil( size(finalBam.bam, "GB") )  * 2 ) + 4
        
        call alignmentAnalysis.UpdateBamSampleName {
            input:
                finalBam = finalBam,
                referenceFa = referenceFa,
                sampleId = sampleId,
                outputPrefix = sampleId,
                diskSize = renameDiskSize
        }
    }
    
    Bam sampleBam = select_first([UpdateBamSampleName.reheaderBam, finalBam])
    
    output {
        Bam sampleBamMatched = sampleBam
    }
}
