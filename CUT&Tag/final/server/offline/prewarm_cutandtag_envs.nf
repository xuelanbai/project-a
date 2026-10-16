#!/usr/bin/env nextflow

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import nextflow.util.CacheHelper

nextflow.enable.dsl = 2

params.pipeline_dir = null
params.cache_dir = null
params.install = false

def hasCondaPackage(Path meta, String name) {
    def matches = Files.newDirectoryStream(meta, "${name}-*.json")
    try {
        return matches.iterator().hasNext()
    }
    finally {
        matches.close()
    }
}

def isEnvironmentReady(Path prefix, String spec) {
    Path meta = prefix.resolve('conda-meta')
    if (!Files.isRegularFile(meta.resolve('history'))) {
        return false
    }
    for (String token : spec.split(/\s+/)) {
        String name = token.replaceFirst(/^.*::/, '').split('=')[0]
        if (!hasCondaPackage(meta, name)) {
            return false
        }
    }
    if (spec.startsWith('bioconda::seacr=')) {
        for (String name : ['SEACR_1.3.sh', 'R', 'bedtools']) {
            if (!Files.isExecutable(prefix.resolve('bin').resolve(name))) {
                return false
            }
        }
    }
    return true
}

def prewarm(String pipelineDir, String cacheDir, boolean install) {
    if (!pipelineDir || !cacheDir) {
        throw new IllegalArgumentException('Specify --pipeline_dir and --cache_dir')
    }

    Path modules = Paths.get(pipelineDir, 'modules')
    Path cache = Paths.get(cacheDir)
    if (!Files.isDirectory(modules) || !Files.isDirectory(cache)) {
        throw new IllegalArgumentException("Missing modules or cache directory: ${modules}, ${cache}")
    }

    // Verify that this Nextflow installation calculates the same Conda-cache
    // name observed in the failed 43984 run. Stop before any installation if not.
    String seacr = 'bioconda::seacr=1.3 conda-forge::r-base=4.0.2 bioconda::bedtools=2.30.0'
    String observed = 'e978bebf288708a606a2e98bd44d43b1'
    String calculated = CacheHelper.hasher(seacr).hash().toString()
    if (calculated != observed) {
        throw new IllegalStateException("Conda hash mismatch: expected ${observed}, got ${calculated}")
    }

    Set<String> specs = new TreeSet<String>()
    def paths = Files.walk(modules)
    try {
        paths.forEach { Path path ->
            if (Files.isRegularFile(path) && path.toString().endsWith('.nf')) {
                Files.readAllLines(path).each { String line ->
                    def match = line =~ /^\s*conda\s+["']([^"']+)["']\s*$/
                    if (match.matches()) {
                        specs.add(match.group(1))
                    }
                }
            }
        }
    }
    finally {
        paths.close()
    }

    if (specs.isEmpty()) {
        throw new IllegalStateException("No literal Conda specifications found under ${modules}")
    }

    int ready = 0
    int missing = 0
    int incomplete = 0
    Map<String,Path> toInstall = new LinkedHashMap<String,Path>()
    println "Found ${specs.size()} distinct Conda specifications; install=${install}"

    for (String spec : specs) {
        Path prefix = cache.resolve('env-' + CacheHelper.hasher(spec).hash().toString())

        if (Files.isDirectory(prefix) && isEnvironmentReady(prefix, spec)) {
            println "READY   ${prefix}  ${spec}"
            ready++
            continue
        }

        if (Files.exists(prefix)) {
            println "INCOMPLETE ${prefix}  ${spec}"
            incomplete++
            continue
        }

        println "MISSING ${prefix}  ${spec}"
        missing++
        toInstall.put(spec, prefix)
    }

    println "READY=${ready}; MISSING=${missing}; INCOMPLETE=${incomplete}; install=${install}"
    if (incomplete > 0 && install) {
        throw new IllegalStateException('Incomplete environments must be inspected and moved aside before installation')
    }
    if (!install) {
        return
    }

    for (Map.Entry<String,Path> entry : toInstall.entrySet()) {
        String spec = entry.key
        Path prefix = entry.value
        List<String> command = new ArrayList<String>([
            'mamba', 'create', '--yes', '--prefix', prefix.toString(),
            '-c', 'conda-forge', '-c', 'bioconda'
        ])
        command.addAll(spec.split(/\s+/).toList())
        int status = new ProcessBuilder(command).inheritIO().start().waitFor()
        if (status != 0 || !isEnvironmentReady(prefix, spec)) {
            throw new IllegalStateException("Conda creation failed or incomplete for ${prefix}; exit=${status}")
        }
    }

    println 'Installation finished; run check mode again to verify every environment.'
}

workflow {
    prewarm(
        params.pipeline_dir?.toString(),
        params.cache_dir?.toString(),
        params.install.toString().toBoolean()
    )
}
