version 1.0

import "wdl_structs.wdl"

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2021) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Minita Shah (mshah@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================

task BedToIntervalList {
    input {
        Int memoryGb = 7
        Int diskSize = 10

        ## Inputs 
        File inputBed
        File? dict
        String outputIntervalListPath = sub(basename(inputBed), ".bed$", ".interval_list")

    }

    command {
            /gatk/gatk \
            --java-options "-Xmx7G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            BedToIntervalList \
            -I ~{inputBed} \
            -O ~{outputIntervalListPath}  \
            -SD ~{dict}
    }

    output {
        File outputIntervalList = "~{outputIntervalListPath}"
    }

    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "us.gcr.io/broad-gatk/gatk:4.1.8.0"
        runtime_minutes: "90"
    }
}

task IntervalListToBed {
    input {
        Int memoryGb = 7
        Int diskSize = 10

        ## Inputs 
        File inputIntervalList
        String outputBedPath = sub(basename(inputIntervalList), ".interval_list$", ".bed")

    }

    command {
            /gatk/gatk \
            --java-options "-Xmx7G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            IntervalListToBed \
            -I ~{inputIntervalList} \
            -O ~{outputBedPath}
    }

    output {
        File outputBed = "~{outputBedPath}"
    }

    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "us.gcr.io/broad-gatk/gatk:4.1.8.0"
        runtime_minutes: "90"
    }
}

workflow MatchIntervalListDict {
    input {
        Array[File]+ inputIntervalLists
        IndexedReference referenceFa
        String name
        
    }
    
    if ( size(referenceFa.dict) > 0 ) {
        scatter(i in range(length(inputIntervalLists))) {
            call IntervalListToBed {
                input:
                    inputIntervalList = inputIntervalLists[i]
            }
            
            String outputIntervalListBasename = sub(basename(inputIntervalLists[i]), ".interval_list$", "")
            
            call BedToIntervalList {
                input:
                    inputBed = IntervalListToBed.outputBed,
                    dict = referenceFa.dict,
                    outputIntervalListPath = "~{outputIntervalListBasename}.~{name}.~{i}.interval_list"
            }
        }
     }
    
    output {
        Array[File]? outputIntervalLists = BedToIntervalList.outputIntervalList
    }
}
