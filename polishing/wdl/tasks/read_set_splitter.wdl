version 1.0

workflow runReadSetSplitter{
    call readSetSplitter
}

task readSetSplitter {
    input {
        Array[File] readFastqs
        Int splitNumber
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=1024
        # pigz と split --filter が必要（tpesout/hpp_base には pigz が無い）
        String dockerImage="mobinasri/bio_base@sha256:8f724fbb8801e3a8e73dfa8bcabebd9188d769db5e3343cd6763218d266e1063"
        Int preemptible=2
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # 入力は gzip 圧縮のまま渡されることがあるので展開しながら読み
        # （pigz -dc -f は非圧縮ならそのまま出す）、分割結果も gzip 圧縮で書く。
        # 非圧縮の分割 fastq を書き出すとリード全量の数倍のディスクを消費するため。
        # 出力名は xaa.fq.gz 等（aligner が ${fileBasename%.*.*} で xaa.bam にする）。
        readsPerFile=$(pigz -dc -f ~{sep=" " readFastqs} | wc -l | awk '{print int($1/(4 * ~{splitNumber})) + 1}')
        linesPerFile=$(( readsPerFile * 4 ))
        pigz -dc -f ~{sep=" " readFastqs} | split -l $linesPerFile --filter='pigz -p~{threadCount} > $FILE.gz' --additional-suffix ".fq"
        for i in *.fq.gz;do du -s -BG $i | sed 's/G.*//';done > outputsize.txt
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        Array[File] splitReadFastqs = glob("*.fq.gz")
        Array[Float] readSizes = read_lines("outputsize.txt")
    }
}
