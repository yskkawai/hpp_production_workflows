version 1.0

task max {
    input {
        Array[Int] integers
    }

    command <<<
        max=~{integers[0]}
        for i in ~{sep=" " integers}; do
            if [[ "$i" -gt "$max" ]]; then
                max=$i
            fi
        done
        echo $max
    >>>

    output {
        Int value = read_int(stdout())
    }

    runtime {
        docker: "tpesout/hpp_base@sha256:2e21cc70c7e74e580ea4b9a6c406e6a801e10fbd916113ab53500f61d9eb3440"
    }
}

task sum {
    input {
        Array[Int] integers
    }

    command <<<
        echo $((~{sep="+" integers}))
    >>>

    output {
        Int value = read_int(stdout())
    }

    runtime {
        docker: "tpesout/hpp_base@sha256:2e21cc70c7e74e580ea4b9a6c406e6a801e10fbd916113ab53500f61d9eb3440"
        preemptible: 1
    }
}
