nextflow.enable.dsl=2

params.manifest = null
params.vcf = null
params.reference = null
params.sample_id = 'PSEUDONYM_001'
params.arrayevidence_bin = null
params.python_bin = 'python3'
params.staging_dir = null

process BUILD_EVIDENCE {
    tag 'synthetic_build'

    input:
    path manifest
    path vcf
    path reference
    path reference_index

    output:
    path 'bundle'

    script:
    """
    ${params.arrayevidence_bin} build \
      --manifest ${manifest} \
      --vcf ${vcf} \
      --reference ${reference} \
      --sample-id ${params.sample_id} \
      --output bundle
    """
}


process VERIFY_EVIDENCE {
    stageInMode 'copy'
    tag 'verify_bundle'

    input:
    path bundle

    output:
    path 'verification.ok'

    script:
    """
    ${params.arrayevidence_bin} verify bundle
    touch verification.ok
    """
}

process QC_EVIDENCE {
    stageInMode 'copy'
    tag 'qc_bundle'

    input:
    path bundle
    path verification_marker

    output:
    path 'qc.json'

    script:
    """
    test -f verification.ok
    ${params.arrayevidence_bin} qc bundle > qc.json
    """
}


process STAGE_EVIDENCE {
    stageInMode 'copy'
    tag 'validated_stage'

    input:
    path bundle
    path verification_marker
    path qc_report

    output:
    path 'finalization.ok'

    script:
    def helper = new File(
        workflow.projectDir.parent.toString(),
        'scripts/finalize_staging.py'
    ).absolutePath

    """
    test -f verification.ok
    test -s qc.json

    ${params.python_bin} "${helper}" \
      bundle \
      qc.json \
      "${params.staging_dir}" \
      "${params.arrayevidence_bin}" \
      "${workflow.sessionId}"

    touch finalization.ok
    """
}

workflow {
    if (!params.manifest || !params.vcf ||
        !params.reference || !params.arrayevidence_bin ||
        !params.staging_dir) {
        error 'Missing required workflow parameters'
    }

    // STAGING_PREFLIGHT: reserve synchronously before process submission.
    def reservationHelper = new File(
        workflow.projectDir.parent.toString(),
        'scripts/reserve_staging.py'
    ).absolutePath

    def reservationCommand = [
        params.python_bin.toString(),
        reservationHelper,
        params.staging_dir.toString(),
        workflow.sessionId.toString(),
        workflow.resume.toString()
    ]

    def reservationProcess = new ProcessBuilder(reservationCommand)
        .redirectErrorStream(true)
        .start()

    def reservationOutput = reservationProcess.inputStream.getText('UTF-8')
    def reservationStatus = reservationProcess.waitFor()

    if (reservationStatus != 0) {
        error "STAGING_PREFLIGHT: reservation failed: ${reservationOutput.trim()}"
    }

    bundles = BUILD_EVIDENCE(
        file(params.manifest),
        file(params.vcf),
        file(params.reference),
        file("${params.reference}.fai")
    )
    verified = VERIFY_EVIDENCE(bundles)
    qc_reports = QC_EVIDENCE(bundles, verified)
    STAGE_EVIDENCE(bundles, verified, qc_reports)
}
