version 1.0

import "../../../QC/wdl/tasks/extract_reads_toGZ.wdl" as extractReadsToGZ_t
import "filter_hifi_adapter.wdl" as adapter_t
import "seqkit_filter_fastq.wdl" as seqkit_filter_wf

workflow runTrioHifiasm{
    input {
        File paternalYak
        File maternalYak
        Array[File] childReadsHiFi
        Array[File] childReadsONT=[]
        Int? homCov
        Float genomeSizeGb = 3.1
        Int minOntReadLength=50000
        String childID
        String? hifiasmExtraOptions
        File? inputBinFilesTarGz
        File? referenceFasta
        Boolean filterAdapters
        Int? min_ont_qscore
        Array[Float] offsetMem = [40, 20, 20]
        Array[Float] memCovRatios = [4.7, 3.8, 3.6]
        String excludeStringReadExtraction=""
        Int threadCount
        Int preemptible
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "humanpangenomics/hifiasm@sha256:50483d0757e5ce1eaf68f33da8ac8087f8b40f1d1310f60100945fc3d56721c5"
        String zones = "us-west2-a"
    }

    scatter (readFile in childReadsHiFi) {
        call extractReadsToGZ_t.extractReadstoGZ as childReadsHiFiExtractedGz {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                excludeString=excludeStringReadExtraction,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
        if (filterAdapters){
            call adapter_t.cutadapt as filterAdapterHiFi {
                input:
                    readFastqGz = childReadsHiFiExtractedGz.extractedRead,
                    diskSizeGB = fileExtractionDiskSizeGB
            }
        }
        File hifiProcessedFastqGz = select_first([filterAdapterHiFi.filteredReadFastqGz, childReadsHiFiExtractedGz.extractedRead])
    }

    Float childReadHiFiSize = size(hifiProcessedFastqGz, 'GB')

    # If homCov is not given, estimate the HiFi read depth (total bases / genome size)
    # to size the memory of step 1. The --hom-cov value for the later steps is taken
    # from hifiasm's own estimate in step 1 (peak_hom).
    if (!defined(homCov)) {
        scatter (fastqGz in hifiProcessedFastqGz) {
            call countBases {
                input:
                    readFastqGz = fastqGz,
                    dockerImage = dockerImage
            }
        }
        call estimateDepth {
            input:
                baseCounts = countBases.bases,
                genomeSizeGb = genomeSizeGb,
                dockerImage = dockerImage
        }
    }
    Int step1Cov = select_first([homCov, estimateDepth.depth])


    # if ONT reads are provided
    if (length(childReadsONT) != 0){
        scatter (readFile in childReadsONT) {
            call extractReadsToGZ_t.extractReadstoGZ as childReadsOntExtractedGz {
                input:
                    readFile=readFile,
                    referenceFasta=referenceFasta,
                    memSizeGB=4,
                    threadCount=4,
                    excludeString=excludeStringReadExtraction,
                    diskSizeGB=fileExtractionDiskSizeGB,
                    dockerImage=dockerImage
            }

            # filter ONT reads for quality
            if (defined(min_ont_qscore)) {
                call seqkit_filter_wf.filter_fastq as qscore_filter_ont {
                    input:
                        input_fastq  = childReadsOntExtractedGz.extractedRead,
                        min_size     = minOntReadLength,
                        min_q        = min_ont_qscore
                }
            }
            File ont_reads = select_first([qscore_filter_ont.filteredFastq, childReadsOntExtractedGz.extractedRead])
         }
         ## get size of all ultralong reads
         Float ul_size_all = size(ont_reads, 'GB')
    }
    # if no ONT data is provided then it would be zero
    Float readULSize = select_first([ul_size_all, 0])

    call trioHifiasm as hifiasmStep1{
        input:
            childReadsHiFi=hifiProcessedFastqGz,
            homCov = homCov,
            childID=childID,
            extraOptions="--bin-only --telo-m CCCTAA",
            memSizeGB=ceil(memCovRatios[0] * step1Cov + offsetMem[0]),
            threadCount=threadCount,
            diskSizeGB= floor((childReadHiFiSize + readULSize) * 2.5) + 64,
            preemptible=preemptible,
            dockerImage=dockerImage,
            zones = zones
    }
    # homCov given by the user, otherwise hifiasm's estimate from step 1
    Int homCovResolved = select_first([homCov, hifiasmStep1.peakHom])

    # If ONT (ultralong) reads are provided, run step2 and step3.
    # If not, skip them and run step4 directly using step1's bin files.
    if (length(childReadsONT) != 0) {
        call trioHifiasm as hifiasmStep2{
            input:
                childReadsUL=ont_reads,
                paternalYak=paternalYak,
                maternalYak=maternalYak,
                homCov = homCovResolved,
                minOntReadLength = minOntReadLength,
                childID=childID,
                extraOptions="--bin-only --telo-m CCCTAA",
                inputBinFilesTarGz=hifiasmStep1.outputBinFiles,
                memSizeGB=ceil(memCovRatios[1] * homCovResolved + offsetMem[1]),
                threadCount=threadCount,
                diskSizeGB= floor((childReadHiFiSize + readULSize) * 2.5) + 128,
                preemptible=preemptible,
                dockerImage=dockerImage,
                zones = zones
        }
        call trioHifiasm as hifiasmStep3{
            input:
                childReadsUL=ont_reads,
                homCov = homCovResolved,
                minOntReadLength = minOntReadLength,
                childID=childID,
                extraOptions="--telo-m CCCTAA",
                inputBinFilesTarGz=hifiasmStep2.outputBinFiles,
                memSizeGB=ceil(memCovRatios[2] * homCovResolved + offsetMem[2]),
                threadCount=threadCount,
                diskSizeGB= floor((childReadHiFiSize + readULSize) * 2.5) + 128,
                preemptible=preemptible,
                dockerImage=dockerImage,
                zones = zones
        }
    }

    # ONT reads are NOT provided -> skip step2/3 and run step4
    if (length(childReadsONT) == 0) {
        call trioHifiasm as hifiasmStep4{
            input:
                paternalYak=paternalYak,
                maternalYak=maternalYak,
                homCov = homCovResolved,
                childID=childID,
                extraOptions="--telo-m CCCTAA",
                inputBinFilesTarGz=hifiasmStep1.outputBinFiles,
                memSizeGB=ceil(memCovRatios[1] * homCovResolved + offsetMem[1]),
                threadCount=threadCount,
                diskSizeGB= floor((childReadHiFiSize + readULSize) * 2.5) + 128,
                preemptible=preemptible,
                dockerImage=dockerImage,
                zones = zones
        }
    }

    output {
        File outputPaternalGfa = select_first([hifiasmStep3.outputPaternalGfa, hifiasmStep4.outputPaternalGfa])
        File outputMaternalGfa = select_first([hifiasmStep3.outputMaternalGfa, hifiasmStep4.outputMaternalGfa])
        File outputPaternalContigGfa = select_first([hifiasmStep3.outputPaternalContigGfa, hifiasmStep4.outputPaternalContigGfa])
        File outputMaternalContigGfa = select_first([hifiasmStep3.outputMaternalContigGfa, hifiasmStep4.outputMaternalContigGfa])
        File outputRawUnitigGfa = select_first([hifiasmStep3.outputRawUnitigGfa, hifiasmStep4.outputRawUnitigGfa])
        File outputBinFiles = select_first([hifiasmStep3.outputBinFiles, hifiasmStep4.outputBinFiles])
        Int homCovUsed = homCovResolved
    }
}

# hifiasm steps
# 1st: pass HiFi; no need to pass UL (extraOptions="--bin-only")
# 2nd: pass UL and both yak files (extraOptions="--bin-only")
# 3rd: pass UL no need to pass yak files
task trioHifiasm {
    input{
        File? paternalYak
        File? maternalYak
        Array[File]? childReadsHiFi
        Array[File]? childReadsUL
        Int? homCov
        Int minOntReadLength = 50000
        String childID
        String? extraOptions
        File? inputBinFilesTarGz
        File? referenceFasta
        # runtime configurations
        Int memSizeGB
        Int threadCount
        Int diskSizeGB
        Int preemptible
        String dockerImage
        String zones
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## if bin files are given we have to extract them in the directory where hifiasm is being run
        ## this enables hifiasm to skip the time-consuming process of finding read overlaps
        if [ ! -v ~{inputBinFilesTarGz} ]; then
            tar -xzf ~{inputBinFilesTarGz} --strip-components 1
            rm -rf ~{inputBinFilesTarGz} || true
        fi

        # keep a copy of hifiasm's log (stderr) to extract its coverage estimate
        exec 2> >(tee -a hifiasm.log >&2)

        # make a fake fastq to use for step 2 and 3
        printf "@fake\nA\n+\nI\n" > fake.fq

        ## run trio hifiasm https://github.com/chhylp123/hifiasm
        # If ONT ultra long reads are provided
        if [[ -n "~{sep="" childReadsUL}" ]]; then
            if [[ -n "~{paternalYak}" ]]; then
                # Run step 2
                hifiasm \
                    ~{extraOptions} \
                    -o ~{childID} \
                    --ul ~{sep="," childReadsUL} \
                    --ul-cut ~{minOntReadLength} \
                    ~{if select_first([homCov, 0]) > 0 then "--hom-cov " + homCov else ""} \
                    -t~{threadCount} \
                    -1 ~{paternalYak} \
                    -2 ~{maternalYak} \
                    fake.fq
            else
                # Keep only necessary bin files for step 3
                mkdir kept_bin_files
                mv *.ec.bin *.ovlp.reverse.bin *.ovlp.source.bin *.hap1.phase.bin *.hap2.phase.bin *.ul.ovlp.bin kept_bin_files
                rm -rf *.bin
                mv kept_bin_files/* .
                rm -rf kept_bin_files

                # Run step 3
                hifiasm \
                    ~{extraOptions} \
                    -o ~{childID} \
                    --ul ~{sep="," childReadsUL} \
                    ~{if select_first([homCov, 0]) > 0 then "--hom-cov " + homCov else ""} \
                    --ul-cut ~{minOntReadLength} \
                    --dual-scaf \
                    -t~{threadCount} \
                    -3 ~{childID}.hap1.phase.bin \
                    -4 ~{childID}.hap2.phase.bin \
                    fake.fq
            fi
        else
            # If no UL reads are provided, run either step 1 (HiFi-only)
            # or step 4 (yak + dual-scaf; uses bins from step 1).
            if [[ -n "~{paternalYak}" ]]; then
                # Run step 4
                hifiasm \
                    ~{extraOptions} \
                    -o ~{childID} \
                    ~{if select_first([homCov, 0]) > 0 then "--hom-cov " + homCov else ""} \
                    --dual-scaf \
                    -t~{threadCount} \
                    -1 ~{paternalYak} \
                    -2 ~{maternalYak} \
                    fake.fq
            else
                # Run step 1
                hifiasm \
                    ~{extraOptions} \
                    -o ~{childID} \
                    -t~{threadCount} \
                    ~{sep=" " childReadsHiFi}
            fi
        fi

        # hifiasm's estimate of the homozygous k-mer peak (last value reported); 0 if absent
        sed -n 's/^\[M::[a-z_]*\] peak_hom: \([0-9][0-9]*\);.*/\1/p' hifiasm.log | tail -n 1 > peak_hom.txt
        [ -s peak_hom.txt ] || echo 0 > peak_hom.txt

        #Move bin and gfa files to saparate folders and compress them
        mkdir ~{childID}.raw_unitig_gfa
        mkdir ~{childID}.pat.contig_gfa
        mkdir ~{childID}.mat.contig_gfa
        mkdir ~{childID}.binFiles

        # before hardlinking gfa files to the corresponding directory make sure they exist
        # first and second step of hifiasm does not output gfa files
        if [[ -n $(find . -maxdepth 1 -iname "*.hap1.p_ctg.*") ]];
        then
            ln ~{childID}.dip.r_utg.* ~{childID}.raw_unitig_gfa
            ln *.hap1.p_ctg.* ~{childID}.pat.contig_gfa
            ln *.hap2.p_ctg.* ~{childID}.mat.contig_gfa
        else
            # To avoid making a separete new task for steps 1 and 2 of hifiasm
            # we make empty gfa files since we cannot have optional outputs in wdl
            touch empty.hap1.p_ctg.gfa
            touch empty.hap2.p_ctg.gfa
        fi

        ln *.bin ~{childID}.binFiles


        # make archives
        tar -cf ~{childID}.raw_unitig_gfa.tar ~{childID}.raw_unitig_gfa
        tar -cf ~{childID}.pat.contig_gfa.tar ~{childID}.pat.contig_gfa
        tar -cf ~{childID}.mat.contig_gfa.tar ~{childID}.mat.contig_gfa
        tar -cf ~{childID}.binFiles.tar ~{childID}.binFiles

        # compress
        pigz -p~{threadCount} ~{childID}.raw_unitig_gfa.tar
        pigz -p~{threadCount} ~{childID}.pat.contig_gfa.tar
        pigz -p~{threadCount} ~{childID}.mat.contig_gfa.tar
        pigz -p~{threadCount} ~{childID}.binFiles.tar
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        cpuPlatform: "Intel Cascade Lake"
        zones : zones
    }

    output {
        File outputPaternalGfa = glob("*.hap1.p_ctg.gfa")[0]
        File outputMaternalGfa = glob("*.hap2.p_ctg.gfa")[0]
        File outputPaternalContigGfa = "~{childID}.pat.contig_gfa.tar.gz"
        File outputMaternalContigGfa = "~{childID}.mat.contig_gfa.tar.gz"
        File outputRawUnitigGfa = "~{childID}.raw_unitig_gfa.tar.gz"
        File outputBinFiles = "~{childID}.binFiles.tar.gz"
        Int peakHom = read_int("peak_hom.txt")
    }
}

task countBases {
    input {
        File readFastqGz
        Int memSizeGB = 4
        Int threadCount = 4
        String dockerImage
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        pigz -dc -p~{threadCount} ~{readFastqGz} | awk 'NR % 4 == 2 {s += length($0)} END {printf "%.0f\n", s}' > bases.txt
    >>>
    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        preemptible: 1
    }
    output {
        Float bases = read_float("bases.txt")
    }
}

task estimateDepth {
    input {
        Array[Float] baseCounts
        Float genomeSizeGb
        String dockerImage
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        awk -v g=~{genomeSizeGb} '{s += $1} END {d = s / (g * 1e9); printf "%d\n", (d == int(d)) ? d : int(d) + 1}' ~{write_lines(baseCounts)} > depth.txt
        echo "estimated HiFi depth: $(cat depth.txt)x" >&2
    >>>
    runtime {
        docker: dockerImage
        memory: "1 GB"
        cpu: 1
        preemptible: 1
    }
    output {
        Int depth = read_int("depth.txt")
    }
}
