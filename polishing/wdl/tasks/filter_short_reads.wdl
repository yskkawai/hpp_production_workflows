version 1.0

import "extract_reads.wdl" as extractReads_t


workflow FilterShortReads {
    input {
        Array[File] readFiles
        Int minReadLength
        File? referenceFasta
    }

    scatter (readFile in readFiles){
        # fastq.gz は filterShortReads が展開しながら直接読むので extractReads
        # （gunzip で非圧縮 fastq を丸ごと書き出す）を通さない。非圧縮 fastq が
        # 1ファイル100GB超になりディスクと時間を大きく消費するため。
        # bam/cram 等それ以外は従来どおり extractReads で fastq にする。
        String readBaseName = basename(readFile)
        Boolean isFastqGz = sub(readBaseName, "\\.(fastq|fq)\\.gz$", "") != readBaseName

        if (!isFastqGz) {
            call extractReads_t.extractReads as extractReads {
                input:
                    readFile=readFile,
                    referenceFasta=referenceFasta,
                    memSizeGB=4,
                    threadCount=4,
                    diskSizeGB=ceil(3 * size(readFile, "GB")) + 64,
                    dockerImage="mobinasri/bio_base:v0.2"
            }
        }
        File readFastq = select_first([extractReads.extractedRead, readFile])

        call filterShortReads{
            input:
                readFastq = readFastq,
                diskSizeGB = ceil(3 * size(readFastq, "GB")) + 64,
                minReadLength = minReadLength
        }
    }

    output {
        Array[File] longReadFastqGzArray = filterShortReads.longReadFastqGz
    }
}


# This task does not work properly if sequences have methylation tags
task filterShortReads {
    input{
        File readFastq
        Int minReadLength
        # runtime configurations
        Int memSizeGB=8
        Int threadCount=4
        Int diskSizeGB=512
        Int preemptible=1
        String dockerImage="mobinasri/bio_base:latest"
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        FILENAME=$(basename -- "~{readFastq}")
        # .gz は外してから拡張子を落とす（extractReads 経由の X.fastq と
        # 直接入力の X.fastq.gz で出力名 X.gt_*kb.fastq.gz を揃える）
        FILENAME="${FILENAME%.gz}"
        PREFIX="${FILENAME%.*}"
        minLenKb=$(echo ~{minReadLength} | awk '{printf "%.0f",$1/1e3}')
        # gzip 圧縮された入力は展開しながら読む（pigz -dc -f は非圧縮ならそのまま出す）
        # filter reads shorter than minReadLength
        pigz -dc -f ~{readFastq} | awk 'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0&&length(b)>~{minReadLength}{print a"\n"b"\n"c"\n"$0;}' | pigz -p8 - > ${PREFIX}.gt_${minLenKb}kb.fastq.gz
        OUTPUTSIZE=`du -s -BG *.fastq.gz | sed 's/G.*//'`
        echo $OUTPUTSIZE > outputsize
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: preemptible
    }

    output {
        File longReadFastqGz = glob("*.fastq.gz")[0]
        Int fileSizeGB = read_int("outputsize")
    }
}
